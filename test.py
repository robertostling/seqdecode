"""Example of sequential decoding"""

import random, math
import numpy as np
import scipy.stats

import pyconv


# NOTE: the default message length limit is 256 bits, including flush bits
#       i.e. we must ensure that data_bits+flush_bits < 256 or a ValueError
#       exception will be raised by encode() and decode()
#       This limit can be changed in conv.py by changing the default value of
#       MAX_WORDS.
#

#message = b'a message to be encoded and back'
message = b'a message to get'
data_bits = 8*len(message)
data = sum(c << (8*i) for i, c in enumerate(message))

# 0 = use flushing bits
# 1 = make code tail-biting
is_tail_biting = 1

# if is_tail_biting == 1, this should be 0
# Otherwise it should generally be equal to K (which currently is always 32)
flush_bits = 0 if is_tail_biting else 32

# ID of convolutional code, see available_codes in conv.c
# 0 is a systematic code suitable for tail-biting (K = 32, rate = 1/2)
# 1 is a non-systematic code (K = 32, rate = 1/2)
code_id = 0 if is_tail_biting else 1

# These parameters offer a tradeoff between speed and accuracy
max_iter_per_bit = 1000000
log2_stack_size = 14

# Variance of normally distributed noise to add to binary signal
noise_sigma = 0.30

# Encode into list of bits (represented as python ints, 0 or 1)
encoded = pyconv.encode(data, data_bits, flush_bits, is_tail_biting, code_id)

n_correct = 0
n_total = 100

print('Encoded %x (%d bits) into %d bits' % (data, data_bits, len(encoded)))

def distort_bit(bit, sigma):
    """Return a soft bit given the original bit and a noise level"""
    x = float(bit) + np.random.randn()*sigma
    prior = [0.5, 0.5]
    p_0 = prior[0]*scipy.stats.norm.pdf(x, loc=0)
    p_1 = prior[1]*scipy.stats.norm.pdf(x, loc=1)
    return p_1 / (p_0 + p_1)


for i in range(n_total):
    encoded_noise = [distort_bit(bit, noise_sigma) for bit in encoded]
    n_errors = sum(int(x > 0.5) != y for x, y in zip(encoded_noise, encoded))
    decoded = pyconv.decode(
            encoded_noise, data_bits, flush_bits, is_tail_biting, code_id,
            max_iter_per_bit, log2_stack_size)
    if decoded == data: n_correct += 1
    print('  Transmission %d (%d bit errors): %s' % (
        i, n_errors, 'ERROR' if decoded != data else 'OK'))

print('%d of %d messages decoded correctly' % (n_correct, n_total))

