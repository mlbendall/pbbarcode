.PHONY: clean all
SHELL = /bin/bash -e

all: build/sw.so

build/sw.so: sw.c
	mkdir -p ./build;\
	gcc -O4 -DGETPROB -shared -fPIC sw.c -o build/sw.so
clean:
	rm -rf build

