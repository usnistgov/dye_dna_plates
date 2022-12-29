#!/bin/bash
cd ..
python3 -m doctest src/wells.py REPRODUCEME.md -v > test/results.txt