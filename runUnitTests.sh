#!/usr/bin/sh

# remove results from prior tests
rm -f test?.*

# run test script
python ./unitTests.py

# run comparisons
echo `diff -q test1.pdb unitTests.expected/test1.pdb`
echo `diff -q test2.pdb unitTests.expected/test2.pdb`
echo `diff -1 test3.txt unitTests.expected/test3.txt`