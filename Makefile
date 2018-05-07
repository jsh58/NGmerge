PREFIX=/usr/local
DESTDIR=
CC?=gcc

NGmerge: NGmerge.c NGmerge.h
	$(CC) -g -Wall -std=gnu99 -fopenmp -O2 -o NGmerge NGmerge.c -lz

install: NGmerge
	@mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp NGmerge $(DESTDIR)$(PREFIX)/bin

clean:
	-@rm NGmerge 2>/dev/null || true
