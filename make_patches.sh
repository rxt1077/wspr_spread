#!/bin/bash

# Doppler spread calculation patch
diff -ru wsprd.c.orig wsprd.c.spread |
	sed "s:^--- wsprd.c.orig:--- wsjtx/lib/wsprd/wsprd.c:" |
	sed "s:^+++ wsprd.c.nodrift:+++ wsjtx/lib/wsprd/wsprd.c:" > spread.patch

# No drift patch
diff -ru wsprd.c.orig wsprd.c.nodrift |
	sed "s:^--- wsprd.c.orig:--- wsjtx/lib/wsprd/wsprd.c:" |
	sed "s:^+++ wsprd.c.nodrift:+++ wsjtx/lib/wsprd/wsprd.c:" > nodrift.patch

# Both patches combined (they're dissimilar enough to work)
cp wsprd.c.orig wsprd.c.tmp
patch -p3 wsprd.c.tmp spread.patch
patch -p3 wsprd.c.tmp nodrift.patch
diff -ru wsprd.c.orig wsprd.c.tmp |
	sed "s:^--- wsprd.c.orig:--- wsjtx/lib/wsprd/wsprd.c:" |
	sed "s:^+++ wsprd.c.tmp:+++ wsjtx/lib/wsprd/wsprd.c:" > all.patch
rm wsprd.c.tmp
rm wsprd.c.tmp.orig
