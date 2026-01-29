#!/bin/bash
# Clone reference repositories for sagemath-ts development
# These are used for reading source code, not for building

set -e

REFERENCE_DIR="reference"
mkdir -p "$REFERENCE_DIR"
cd "$REFERENCE_DIR"

echo "Cloning reference repositories..."
echo "Note: These are shallow clones to save disk space."
echo ""

# SageMath main repository
if [ ! -d "sage" ]; then
    echo "Cloning sagemath/sage..."
    git clone --depth 1 https://github.com/sagemath/sage.git
else
    echo "sage/ already exists, skipping..."
fi

# cypari2 - Python interface to PARI/GP
if [ ! -d "cypari2" ]; then
    echo "Cloning sagemath/cypari2..."
    git clone --depth 1 https://github.com/sagemath/cypari2.git
else
    echo "cypari2/ already exists, skipping..."
fi

# PARI/GP - Number theory library
if [ ! -d "pari" ]; then
    echo "Cloning pari/pari..."
    git clone --depth 1 https://pari.math.u-bordeaux.fr/git/pari.git
else
    echo "pari/ already exists, skipping..."
fi

# FLINT - Fast Library for Number Theory
if [ ! -d "flint" ]; then
    echo "Cloning flintlib/flint..."
    git clone --depth 1 https://github.com/flintlib/flint.git
else
    echo "flint/ already exists, skipping..."
fi

# NTL - Number Theory Library
if [ ! -d "ntl" ]; then
    echo "Cloning libntl/ntl..."
    git clone --depth 1 https://github.com/libntl/ntl.git
else
    echo "ntl/ already exists, skipping..."
fi

echo ""
echo "Done! Reference repositories are in: $REFERENCE_DIR/"
echo ""
echo "Key paths for number theory:"
echo "  - SageMath rings:      reference/sage/src/sage/rings/"
echo "  - SageMath arith:      reference/sage/src/sage/arith/"
echo "  - SageMath EC:         reference/sage/src/sage/schemes/elliptic_curves/"
echo "  - cypari2 interface:   reference/cypari2/cypari2/"
echo "  - PARI/GP source:      reference/pari/src/"
echo "  - FLINT source:        reference/flint/src/"
echo "  - NTL source:          reference/ntl/src/"
