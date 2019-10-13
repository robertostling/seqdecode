CFLAGS=-lm

all: pyconv.pyx conv.c hash.c mmheap.c mmheap_trickle.c natmap.c
	cythonize -a -i pyconv.pyx 

test_mmheap: test_mmheap.c

test: all
	python3 test.py

clean:
	rm -f pyconv.cpython*.so pyconv.c test_mmheap

