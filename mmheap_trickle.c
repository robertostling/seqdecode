static inline void __TRICKLE_FUN_NAME(mmheap *h, size_t i) {
    /* if no child nodes, do nothing */
    if (exists(h, left_child(i))) {
        /* find m, the index of the best (smallest/largest) child or
         * grandchild of i. The left child is the best so far */
        size_t m = left_child(i);
        mmheap_key_t best = __GET_KEY(h, m);
        /* then explore the children of the left child, if they exist */
        if (exists(h, left_child(left_child(i))) &&
                __IS_BETTER(__GET_KEY(h, left_child(left_child(i))), best))
        {
            m = left_child(left_child(i));
            best = __GET_KEY(h, m);
        }
        if (exists(h, right_child(left_child(i))) &&
                __IS_BETTER(__GET_KEY(h, right_child(left_child(i))), best))
        {
            m = right_child(left_child(i));
            best = __GET_KEY(h, m);
        }

        /* if there is also a right child, explore that and its children */
        if (exists(h, right_child(i))) {
            if (__IS_BETTER(__GET_KEY(h, right_child(i)), best)) {
                m = right_child(i);
                best = __GET_KEY(h, m);
            }
            if (exists(h, left_child(right_child(i))) &&
                    __IS_BETTER(__GET_KEY(h, left_child(right_child(i))),
                                best))
            {
                m = left_child(right_child(i));
                best = __GET_KEY(h, m);
            }
            if (exists(h, right_child(right_child(i))) &&
                    __IS_BETTER(__GET_KEY(h, right_child(right_child(i))),
                                best))
            {
                m = right_child(right_child(i));
                best = __GET_KEY(h, m);
            }
        }
        if (parent(m) == i) {
            /* m is a child of i */
            if (__IS_BETTER(best, __GET_KEY(h, i)))
                __SWAP(h, i, m);
        } else {
            /* m is a grandchild of i */
            const mmheap_key_t worse = __GET_KEY(h, i);
            if (__IS_BETTER(best, worse)) {
                __SWAP(h, i, m);
                if (__IS_WORSE(worse, __GET_KEY(h, parent(m)))) {
                    __SWAP(h, m, parent(m));
                }
                __TRICKLE_FUN_NAME(h, m);
            }
        }
    }
}


static inline void __BUBBLE_FUN_NAME(mmheap *h, size_t i) {
    if (has_grandparent(i)) {
        const size_t m = parent(parent(i));
        if (__IS_BETTER(__GET_KEY(h, i), __GET_KEY(h, m))) {
            __SWAP(h, i, m);
            __BUBBLE_FUN_NAME(h, m);
        }
    }
}

