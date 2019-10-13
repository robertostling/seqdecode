# cython: language_level=3
from libc.stdint cimport uint32_t, uint64_t

cdef extern from "conv.c":
    ctypedef struct conv_code:
        int k
        int n_poly
        uint32_t poly[8]
        int qli_bit

    int max_words
    conv_code available_codes[2]

    void decode_block(
            const conv_code *code,
            const float *encoded,
            uint64_t *decoded,
            size_t n_bits,
            size_t n_flush,
            long int max_iter_per_bit,
            int log2_stack_size,
            int tailbiting)

    void encode_block(
            const conv_code *code,
            const uint64_t *data,
            int *encoded,
            size_t n_bits,
            size_t n_flush,
            int tailbiting)


def encode(data, n_bits, n_flush, tailbiting, code_id):
    # Conservatively large input and output buffers
    cdef uint64_t in_buf[128]
    cdef int out_buf[128*8*64]
    cdef conv_code code = available_codes[code_id]
    cdef int padding_bits = 0

    if 8 + n_bits + (0 if tailbiting else n_flush) >= 64*max_words:
        raise ValueError('Message too long, increase MAX_WORDS in conv.c')

    assert data < (1 << n_bits)

    if n_bits % 64:
        padding_bits = 64 - (n_bits % 64)
        data = data << padding_bits

    for i in range(0, n_bits, 64):
        in_buf[i//64] = ((data >> (n_bits+padding_bits-(64+i)))
                            & 0xffffffffffffffff)

    encode_block(&code, in_buf, out_buf, n_bits, n_flush, tailbiting)

    return [out_buf[i] for i in range((n_bits+n_flush)*code.n_poly)]


def decode(data, n_bits, n_flush, tailbiting, code_id, max_iter_per_bit,
           log2_stack_size):
    cdef float in_buf[128*8*64]
    cdef uint64_t out_buf[128]
    cdef int padding_bits = 0
    cdef conv_code code = available_codes[code_id]

    # TODO: other rules with tail-biting codes?
    if 8 + n_bits + (0 if tailbiting else n_flush) >= 64*max_words:
        raise ValueError('Message too long, increase MAX_WORDS in conv.c')

    for i, x in enumerate(data):
        in_buf[i] = x

    decode_block(
            &code, in_buf, out_buf, n_bits, n_flush, max_iter_per_bit,
            log2_stack_size, tailbiting)

    data = 0
    for i in range(0, n_bits, 64):
        data = (data << 64) | out_buf[i//64]

    if n_bits % 64:
        data = data >> (64 - (n_bits % 64))

    return data

