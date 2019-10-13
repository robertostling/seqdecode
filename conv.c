#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>


typedef struct {
    // constraint length
    int k;
    // number of polynomials
    int n_poly;
    // polynomials, there are n_poly of them. poly[0] == 1 indicates a
    // systematic code.
    uint32_t poly[8];
    // index of quick look-in bit, negative values indicate no quick look-in
    // currently all quick look-in is assumed to be
    //      ((poly[0] ^ poly[1]) >> qli_bit) & 1
    int qli_bit;
} conv_code;

/* Maximum number of message bits. Since they are stored as uint64_t in
 * the stack_item struct, this should be a multiple of 64. If zero-valued
 * flushing (tail bits) are used, they should also fit within this limit. */
#ifndef MAX_WORDS
#define MAX_WORDS   5
#endif
#define MAX_BITS    (64*MAX_WORDS)

int max_words = MAX_WORDS;

conv_code available_codes[2] = {
    // From Table 1 of Johannesson & Ståhl (1999).
    // http://portal.research.lu.se/portal/files/5311959/1058578.pdf
    {   32,     // K = 32
        2,      // rate 1/2
        { 0x00000001, 0xd8c6993b },
        0,      // systematic bit at position 0
    },
    // From fano.c by Phil Karn, KA9Q and Joe Taylor, K1JT.
    // Reference: https://ntrs.nasa.gov/search.jsp?R=19720026674
    {   32,     // K = 32
        2,      // rate 1/2
        { 0xf2d05351, 0xe4613c47 },
        -1,     // no quick look-in
    }
};

/* Basic parity computation method from
 * https://graphics.stanford.edu/~seander/bithacks.html#ParityParallel
static inline unsigned int parity_u32(uint32_t v) {
    v ^= v >> 16;
    v ^= v >> 8;
    v ^= v >> 4;
    return (0x6996U >> (v & 0xf)) & 1;
}
*/

static inline unsigned int encode_state(const conv_code *code, uint32_t s) {
    uint32_t x[code->n_poly];

    for (int i=0; i<code->n_poly; i++)
        x[i] = s & code->poly[i];

    for (int i=0; i<code->n_poly; i++)
        x[i] ^= x[i] >> 16;

    for (int i=0; i<code->n_poly; i++)
        x[i] ^= x[i] >> 8;

    for (int i=0; i<code->n_poly; i++)
        x[i] ^= x[i] >> 4;

    unsigned int r = 0;
    unsigned int pattern = 0x6996;
    unsigned int mask = 1;
    for (int i=0; i<code->n_poly; i++) {
        r |= (pattern >> (x[i] & 0xf)) & mask;
        mask <<= 1;
        pattern <<= 1;
    }
    return r;
}

typedef struct {
    uint64_t history[MAX_WORDS];
    float score;
    int depth;
} stack_item;

typedef stack_item mmheap_item_t;
typedef float mmheap_key_t;
static inline mmheap_key_t mmheap_get_key(const mmheap_item_t *item) {
    return item->score;
}
#include "mmheap.c"

#define MAKE_NAME(NAME) map_u64_u32 ## NAME
#define MAX_FIXED       0
#define INDEX_TYPE      size_t
#define KEY_TYPE        uint64_t
#define VALUE_TYPE      uint32_t
#define EMPTY_KEY       0xffffffffffffffffULL
#define HASH_KEY        hash_u64_u64
#include "hash.c"
#include "natmap.c"

/* Decoder for convolutional codes.
 *
 * encoded contains the soft input symbols, where encoded[i*code->n_poly + j]
 * is the probability of polynomial j having generated a 1 bit at position i.
 *
 * decoded is a buffer for decoded data, where the i:th bit can be accessed as
 * (decoded[i/64] >> (63 - i%64)) & 1
 * that is, most significant bit first order within uint64s.
 *
 * n_bits is the number of original data bits, excluding the all-zero flushing
 * bits. Note that n_bits+n_flush must be at most MAX_BITS, and for 
 *
 * n_flush is the number of zero bits at the end used for flushing, typically
 * this would be equal to the constraint length CODE_K.
 *
 * tailbiting is 0 if tail bits (all 0) are used, or non-zero if the code is
 * tail-biting.
 */
static int decode_block(
        const conv_code *code,
        const float *encoded,
        uint64_t *decoded,
        size_t n_bits,
        size_t n_flush,
        long int max_iter_per_bit,
        int log2_stack_size,
        int tailbiting)
{
    const long int max_iter = n_bits*max_iter_per_bit;
    stack_item item;
    const size_t stack_size = 1 << log2_stack_size;

    if (tailbiting && code->qli_bit < 0) {
        // Tail-biting codes need to be systematic or have quick look-in
        return 1;
    }

    // encoded_joint[(i << code->n_poly) + x] contains the negative
    // log-probability of the i:th codeword being x. This is computed as the
    // log of the product of each individual bit probability.
    float encoded_joint[(n_bits+n_flush) << code->n_poly];
    for (int i=0; i<n_bits+n_flush; i++) {
        for (int x=0; x<(1 << code->n_poly); x++) {
            float p = 1.0f;
            for (int j=0; j<code->n_poly; j++) {
                const float p1 = encoded[i*code->n_poly + j];
                p *= ((x >> j) & 1)? p1: 1.0f-p1;
            }
            const float log_p = logf((p < 1e-10f)? 1e-10f: p);
            encoded_joint[(i << code->n_poly) + x] = -log_p;
        }
    }

    // This will be set to the depth of decoding. In case of tail-biting codes
    // this is several times the length of the data, otherwise it is
    // n_bits + n_flush.
    int decode_depth;

    mmheap *h = mmheap_create(stack_size);

    memset(&item, 0, sizeof(item));
    item.score = 0.0f;

    if (tailbiting) {
        // Tail-biting code, which means we are not afforded the luxury of a
        // known initial state. We now need to do the following:
        //
        //  1. find the subsequence of the data (using the systematic bits or
        //     quick look-in) with the least uncertainty
        //  2. start decoding from there, and keep decoding until the data has
        //     been processed at least 4 times, so that we have time to move
        //     away from errors in the initial estimate

        // Temporary heap used for constructing the list of possible initial
        // states.
        mmheap *h_init = mmheap_create(stack_size);
        // Keys in this map indicate which states have been explored.
        // Values are not used (currently set to zero).
        struct map_u64_u32 state_map;
        map_u64_u32_create(&state_map);

        // Soft bits and their logs for systematic/quick look-in bits.
        float soft_log_p0[n_bits+n_flush];
        float soft_log_p1[n_bits+n_flush];
        float soft_p0[n_bits+n_flush];
        float soft_p1[n_bits+n_flush];

        for (int i=0; i<n_bits; i++) {
            // qli_bit is assumed to be non-negative
            const int j = i*code->n_poly;
            const int l = (i+code->qli_bit) % n_bits;
            float p0, p1;
            if (code->poly[0] == 1) {
                // Systematic code, so simply copy the systematic bit
                p1 = encoded[j + 0];
            } else {
                // Quick look-in by XOR of first two polynomials
                p1 = (encoded[j + 0] * (1.0f - encoded[j + 1]))
                   + ((1.0f - encoded[j + 0]) * encoded[j + 1]);
            }
            p0 = 1.0f - p1;
            soft_p0[l] = p0;
            soft_p1[l] = p1;
            // Compute negative logs (with clipping)
            soft_log_p0[l] = -logf((p0 < 1e-10f)? 1e-10f: p0);
            soft_log_p1[l] = -logf((p1 < 1e-10f)? 1e-10f: p1);
        }

        // Find the starting offset which will result in the lowest uncertainty
        // (minimum entropy) over the CODE_K data/quick look-in bits following.
        float best_entropy = FLT_MAX;
        int best_start = -1;
        for (int i=0; i<n_bits; i++) {
            float entropy = 0.0f;
            for (int j=0; j<code->k; j++) {
                entropy += soft_p0[(i+j)%n_bits]*soft_log_p0[(i+j)%n_bits] +
                           soft_p1[(i+j)%n_bits]*soft_log_p1[(i+j)%n_bits];
            }
            //printf("entropy at %d: %g\n", i, entropy);
            if (entropy < best_entropy) {
                best_entropy = entropy;
                best_start = i;
            }
        }

        // Populate the heap with a set of high-probability states given the quick
        // look-in data. This is done by starting from the most likely state
        // (simple hard bit-wise decision), exploring from there by flipping one
        // bit at a time. A separate heap as well as a hash table are used for
        // keeping track of the most likely unexplored state.

        // Get best guess.
        item.depth = best_start;
        for (int i=0; i<code->k; i++) {
            const int j = (best_start+i) % n_bits;
            if (soft_p1[j] > soft_p0[j]) {
                item.history[MAX_WORDS-1] |= 1ULL << ((code->k-1)-i);
                item.score += soft_log_p1[j];
            } else {
                item.score += soft_log_p0[j];
            }
        }
        // Insert initial guess into both heaps.
        // The best-scoring state in h_init will be removed at each iteration,
        // and its nearest neighbors (Hamming distance = 1) that are not in
        // the table of explored states (state_map) will be added.  Items will
        // only be removed from h if the heap is full and there is a
        // better-scoring state available.
        mmheap_insert(h_init, &item);
        mmheap_insert(h, &item);
        map_u64_u32_add(&state_map, item.history[MAX_WORDS-1], 0);

        // TODO: how many iterations are really needed?
        // TODO: what size of h_init is needed?
        for (size_t i=0; i<stack_size*2; i++) {
            if (mmheap_delete_min(h_init, &item)) break;
            for (int j=0; j<code->k; j++) {
                const int l = (best_start + j) % n_bits;
                const uint64_t candidate =
                    item.history[MAX_WORDS-1] ^ (1ULL << ((code->k-1) - j));
                if (map_u64_u32_get_ptr(&state_map, candidate) == NULL) {
                    map_u64_u32_add(&state_map, candidate, 0);
                    stack_item new_item = item;
                    new_item.history[MAX_WORDS-1] = candidate;
                    const float new_log_p =
                        ((candidate >> ((code->k-1) - j)) & 1)?
                        soft_log_p1[l] : soft_log_p0[l];
                    const float old_log_p =
                        ((candidate >> ((code->k-1) - j)) & 1)?
                        soft_log_p0[l] : soft_log_p1[l];
                    new_item.score = item.score + new_log_p - old_log_p;

                    if (mmheap_is_full(h)) {
                        const size_t max_index = mmheap_max_index(h);
                        if (new_item.score < h->data[max_index].score) {
                            h->data[max_index] = new_item;
                            mmheap_update_max(h, max_index);
                        }
                    } else {
                        mmheap_insert(h, &new_item);
                    }
                    if (mmheap_is_full(h_init)) {
                        const size_t max_index = mmheap_max_index(h_init);
                        if (new_item.score < h_init->data[max_index].score) {
                            h_init->data[max_index] = new_item;
                            mmheap_update_max(h_init, max_index);
                        }
                    } else {
                        mmheap_insert(h_init, &new_item);
                    }
                }
            }
        }

        mmheap_destroy(h_init);
        decode_depth = n_bits*4 + code->k;
    } else {
        // Not tail-biting, so initialize empty state to zero and start
        // decoding from position zero.
        mmheap_insert(h, &item);
        decode_depth = n_bits + n_flush;
    }

    for (long int iter=0; iter < max_iter; iter++) {
        // Lowest (best) scoring hypothesis, will be replaced in-place with
        // its 0 branch.
        stack_item *branch0 = h->data + mmheap_min_index(h);
        // Make a copy of the parent node before overwriting.
        stack_item parent = *branch0;
        stack_item branch1;

        if (parent.depth == decode_depth) {
            // Found (possible) solution, copy data and return success.
            int bit0 = (tailbiting)? MAX_BITS-(n_bits+code->k) :
                                     MAX_BITS-(n_bits+n_flush);
            int bit1 = (tailbiting)? MAX_BITS-code->k : MAX_BITS-n_flush;
            int n_words = (n_bits+63)/64;
            for (int i=0; i<n_words; i++)
                decoded[i] = 0;

            //printf("decoded %"PRIx64" %"PRIx64"\n",
            //        parent.history[MAX_WORDS-2], parent.history[MAX_WORDS-1]);

            for (int bit=bit0; bit<bit1; bit++) {
                const uint64_t x =
                    (parent.history[bit/64] >> (63 - (bit%64))) & 1;
                decoded[(bit-bit0)/64] |= x << (63 - ((bit-bit0)%64));
            }
            return 0;
        }

        branch0->history[MAX_WORDS-1] = (parent.history[MAX_WORDS-1] << 1) | 0;
        branch1.history[MAX_WORDS-1]  = (parent.history[MAX_WORDS-1] << 1) | 1;
        for (int i=0; i < MAX_WORDS-1; i++)
            branch1.history[i] = branch0->history[i] =
                (parent.history[i] << 1) | ((parent.history[i+1] >> 63) & 1);

        unsigned int encoded_0 =
            encode_state(code, branch0->history[MAX_WORDS-1]);
        unsigned int encoded_1 =
            encode_state(code, branch1.history[MAX_WORDS-1]);

        branch1.depth = branch0->depth = parent.depth + 1;
        const int base = (parent.depth % (n_bits+n_flush)) << code->n_poly;
        branch0->score = parent.score + encoded_joint[base + encoded_0];
        branch1.score = parent.score + encoded_joint[base + encoded_1];

        // Delete parent and insert 0 branch into the heap. This does not
        // change the heap size.
        mmheap_update_min(h);

        if (tailbiting || (parent.depth < n_bits)) {
            // There is only a 1 branch if we have not yet reached the
            // flushing part (where all bits are assumed to be zero),
            // or if we are decoding a tail-biting code.
            if (mmheap_is_full(h)) {
                // If the heap is full, the 1 branch should re-use the space
                // of the maximum (worst) stack item.
                const size_t max_index = mmheap_max_index(h);
                h->data[max_index] = branch1;
                mmheap_update_max(h, max_index);
            } else {
                // If the heap is not full, simply insert the 1 branch.
                mmheap_insert(h, &branch1);
            }
        }
    }
    // No solution found in the time given, return failure.
    return 1;
}

/* Encoder for convolutional codes with up to MAX_BITS bits of data.
 *
 * data contains the input to be encoded, where the i:th bit can be accessed as
 * (data[i/64] >> (63 - i%64)) & 1
 *
 * encoded contains binary output symbols. There will be
 *  (n_bits+n_flush)*code->n_polys
 * bits, one per element in the array.
 *
 * n_bits is the number of original data bits, excluding the all-zero flushing
 * bits. Note that n_bits+n_flush must be at most MAX_BITS.
 *
 * n_flush is the number of zero bits at the end used for flushing, typically
 * this would be equal to the constraint length code->k. For tail-biting codes
 * this should be set to 0.
 *
 * tailbiting is 0 if tail bits (all 0) are used, or non-zero if the code is
 * tail-biting.
 */
static void encode_block(
        const conv_code *code,
        const uint64_t *data,
        int *encoded,
        size_t n_bits,
        size_t n_flush,
        int tailbiting)
{
    uint32_t state = 0;
    if (tailbiting) {
        // Code is tail-biting, so initialize state with end of data
        for (size_t i=n_bits-code->k; i<n_bits; i++) {
            state <<= 1;
            state |= (data[i/64] >> (63 - i%64)) & 1;
        }
    }
    for (size_t i=0; i<n_bits+n_flush; i++) {
        state <<= 1;
        if (i < n_bits) state |= (data[i/64] >> (63 - i%64)) & 1;

        const unsigned int encoded_bits = encode_state(code, state);
        for (int j=0; j<code->n_poly; j++)
            encoded[code->n_poly*i + j] = (encoded_bits >> j) & 1;
    }
}

