#!/bin/bash

echo "Creating ANGEL documentation using doxygen..."
doxygen doxy.conf
echo "Creating PDF..."
cd latex
pdflatex refman.tex
pdflatex refman.tex
cd ..
ln -s latex/refman.pdf manual.pdf
echo "Done. Look at manual.pdf or at ./html/index.html."
