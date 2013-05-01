#!/usr/bin/sh

# remove results from prior tests
rm -f test?.*

# run test script
python ./unitTests.py

# run comparisons
echo "If this is the last non-empty line of output, then your PeptideBuilder installation works as expected."
echo `diff -q test1.pdb unitTests.expected/test1.pdb`
echo `diff -q test2.pdb unitTests.expected/test2.pdb`
echo `diff -q test3.txt unitTests.expected/test3.txt`
echo `diff -q test4.pdb unitTests.expected/test4.pdb`
echo `diff -q test5.pdb unitTests.expected/test5.pdb`
echo `diff -q test6.pdb unitTests.expected/test6.pdb`
echo `diff -q test7.pdb unitTests.expected/test7.pdb`
