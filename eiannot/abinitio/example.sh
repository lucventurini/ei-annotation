#!/usr/bin/env bash

# Convert repeats into a single file

repeats_to_hints.py $REPEATS $REPHINTS

# Then convert proteins to hints

proteins_to_hints.py $PROTEINS $PROTHINTS

# Convert Mikado (HOW?!?!?!!)

mikado_to_hints.py $MIKADO $HINTS

# Convert the FLN


