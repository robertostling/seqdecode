/* min-max heap implementation
 * https://cglab.ca/~morin/teaching/5408/refs/minmax.pdf
 *
 * Before including this file, define the following:
 *
 *  mmheap_item_t       -- type of heap items
 *  mmheap_key_t        -- type of keys
 *  mmheap_get_key      -- function or macro from
 *                         mmheap_item_t* to mmheap_key_t
 */

#include <stdlib.h>

typedef struct {
    size_t size;            /* current size of heap */
    size_t max_size;        /* number of allocated items */
    mmheap_item_t *data;    /* array of max_size allocated items */
} mmheap;

#define __GET_KEY(h, i) mmheap_get_key((h)->data + i)

#define __SWAP(h, i, j) { \
    mmheap_item_t tmp = h->data[i]; \
    h->data[i] = h->data[j]; \
    h->data[j] = tmp; \
    }

mmheap *mmheap_create(size_t max_size) {
    mmheap *h = malloc(sizeof(mmheap));
    if (h == NULL) return NULL;
    h->size = 0;
    h->max_size = max_size;
    h->data = malloc(max_size*sizeof(mmheap_item_t));
    if (h->data == NULL) {
        free(h);
        return NULL;
    }
    return h;
}

void mmheap_destroy(mmheap *h) {
    free(h->data);
    free(h);
}

/* Return 1 if floor(log2(i+1)) is odd, 0 if it is even */
static inline int is_max_level(size_t i) {
    return ((i+1) & 0x5555555555555555UL) <= ((i+1) & 0xaaaaaaaaaaaaaaaaUL);
}

static inline int exists(const mmheap *h, size_t i) {
    return i < h->size;
}

static inline size_t left_child(size_t i) {
    return i*2 + 1;
}

static inline size_t right_child(size_t i) {
    return i*2 + 2;
}

static inline size_t parent(size_t i) {
    return (i-1) >> 1;
}

static inline int has_grandparent(size_t i) {
    return i > 2;
}

static inline int has_parent(size_t i) {
    return i != 0;
}

#define __IS_BETTER(x,y)    ((x)<(y))
#define __IS_WORSE(x,y)     ((x)>(y))
#define __TRICKLE_FUN_NAME  mmheap_trickle_down_min
#define __BUBBLE_FUN_NAME   mmheap_bubble_up_min
#include "mmheap_trickle.c"
#undef __IS_BETTER
#undef __IS_WORSE
#undef __TRICKLE_FUN_NAME
#undef __BUBBLE_FUN_NAME

#define __IS_BETTER(x,y)    ((x)>(y))
#define __IS_WORSE(x,y)     ((x)<(y))
#define __TRICKLE_FUN_NAME  mmheap_trickle_down_max
#define __BUBBLE_FUN_NAME   mmheap_bubble_up_max
#include "mmheap_trickle.c"
#undef __IS_BETTER
#undef __IS_WORSE
#undef __TRICKLE_FUN_NAME
#undef __BUBBLE_FUN_NAME

static void mmheap_trickle_down(mmheap *h, size_t i) {
    if (is_max_level(i))
        mmheap_trickle_down_max(h, i);
    else
        mmheap_trickle_down_min(h, i);
}

static void mmheap_bubble_up(mmheap *h, size_t i) {
    if (is_max_level(i)) {
        if (has_parent(i) && (__GET_KEY(h, i) < __GET_KEY(h, parent(i)))) {
            __SWAP(h, i, parent(i));
            mmheap_bubble_up_min(h, parent(i));
        } else {
            mmheap_bubble_up_max(h, i);
        }
    } else {
        if (has_parent(i) && (__GET_KEY(h, i) > __GET_KEY(h, parent(i)))) {
            __SWAP(h, i, parent(i));
            mmheap_bubble_up_max(h, parent(i));
        } else {
            mmheap_bubble_up_min(h, i);
        }
    }
}

static int mmheap_insert(mmheap *h, const mmheap_item_t *item) {
    if (h->size == h->max_size) return 1;
    const size_t i = h->size++;
    h->data[i] = *item;
    mmheap_bubble_up(h, i);
    return 0;
}

static int mmheap_delete_min(mmheap *h, mmheap_item_t *item) {
    if (h->size == 0) return 1;
    *item = h->data[0];
    if (h->size == 1) {
        h->size = 0;
        return 0;
    }
    h->data[0] = h->data[--h->size];
    mmheap_trickle_down(h, 0);
    return 0;
}

static int mmheap_delete_max(mmheap *h, mmheap_item_t *item) {
    if (h->size == 0) return 1;
    if (h->size == 1) {
        if (item != NULL) *item = h->data[0];
        h->size = 0;
        return 0;
    }
    if (h->size == 2) {
        if (item != NULL) *item = h->data[1];
        h->size = 1;
        return 0;
    }
    if (__GET_KEY(h, 1) > __GET_KEY(h, 2)) {
        if (item != NULL) *item = h->data[1];
        h->data[1] = h->data[--h->size];
        mmheap_trickle_down(h, 1);
    } else {
        if (item != NULL) *item = h->data[2];
        h->data[2] = h->data[--h->size];
        mmheap_trickle_down(h, 2);
    }
    return 0;
}

static inline size_t mmheap_min_index(mmheap *h) {
    return 0;
}

static inline size_t mmheap_max_index(mmheap *h) {
    if (h->size < 3) return h->size - 1;
    return (__GET_KEY(h, 1) > __GET_KEY(h, 2))? 1: 2;
}

static void mmheap_update_min(mmheap *h) {
    mmheap_trickle_down(h, 0);
}

static void mmheap_update_max(mmheap *h, size_t max_index) {
    mmheap_trickle_down(h, max_index);
}

static inline int mmheap_is_full(mmheap *h) {
    return h->size == h->max_size;
}

#undef __GET_KEY
#undef __SWAP

