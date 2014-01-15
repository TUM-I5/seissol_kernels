#!/bin/sh

for i in `find . -name "*.?pp"`; do astyle --options=astylerc $i ; done
