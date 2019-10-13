# Sequential decoding of convolutional codes

This repository contains a C library with Python wrappers for convolutional
error-correcting codes. It currently implements the following:

* Sequential decoding
* Tail-biting or tail bits
* Two K=32 rate=1/2 convolutional codes ([http://portal.research.lu.se/portal/files/5311959/1058578.pdf](Johannesson & Ståhl 1999); [https://ntrs.nasa.gov/search.jsp?R=19720026674](Layland & Lushbaugh 1971) used in WSPR)

Note that decoding with default parameters is fairly slow, in general several
seconds per 256-bit message with tail-biting codes. This library is meant for
exploring weak-signal communications with low transmission rates.

For a documented example, please see `test.py`.

