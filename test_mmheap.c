#include <stdio.h>
#include <math.h>

typedef struct {
    float x;
    int id;
} mmheap_item_t;

typedef float mmheap_key_t;

static inline mmheap_key_t mmheap_get_key(mmheap_item_t *item) {
    return item->x;
}

// #define mmheap_get_key(item) ((item)->x)

#include "mmheap.c"

int main(void) {
    size_t n = 10000000;
    mmheap *h = mmheap_create(2*n);
    mmheap_item_t item;

    float phi = 0.0f;
    for (int i=0; i<2*n; i++) {
        item.x = sinf(phi);
        item.id = i;
        phi += 0.1f;
        mmheap_insert(h, &item);
        //printf("inserted %2d %.3f\n", item.id, item.x);
        /*
        for (size_t j=0; j<h->size; j++) {
            printf(" %.3f", h->data[j].x);
        }
        printf("\n");
        */
    }

    printf("Insertion complete\n");

    phi = 0.0f;
    for (int i=0; i<2*n; i+=10) {
        const size_t max_index = mmheap_max_index(h);
        h->data[max_index].x = sinf(phi);
        mmheap_update_max(h, max_index);
        h->data[0].x = -sinf(phi);
        mmheap_update_min(h);
        phi += 0.01f;
    }
    
    printf("Updates complete\n");

    mmheap_delete_max(h, &item);
    float max_x = item.x;
    mmheap_delete_min(h, &item);
    float min_x = item.x;

    for (int i=1; i<n; i++) {
        mmheap_delete_max(h, &item);
        //printf("%2d %.3f\n", item.id, item.x);
        if (! (item.x <= max_x)) printf("!!! %.3f > %.3f\n", item.x, max_x);
        max_x = item.x;

        mmheap_delete_min(h, &item);
        //printf("%2d %.3f\n", item.id, item.x);
        if (! (item.x >= min_x)) printf("!!! %.3f < %.3f\n", item.x, min_x);
        min_x = item.x;

        /*
        printf("%2d %.3f\n", item.id, item.x);
        for (size_t j=0; j<h->size; j++) {
            printf(" %.3f", h->data[j].x);
        }
        printf("\n");
        //mmheap_delete_min(h, &item);
        mmheap_delete_max(h, &item);
        printf("%2d %.3f\n", item.id, item.x);
        for (size_t j=0; j<h->size; j++) {
            printf(" %.3f", h->data[j].x);
        }
        printf("\n");
        */
    }

    mmheap_destroy(h);
    return 0;
}

