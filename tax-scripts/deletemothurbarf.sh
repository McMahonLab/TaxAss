#!/bin/bash

# delete the mothur barf lines from terminal output to check for errors.

# syntax:
# ./deletemothurbarf.sh filename.txt > newfilename.txt


filename="$1"

# comments are examples of what's removed. the \ allows multiline syntax

# 112000
egrep -v '^\d+$' $filename | \

# 3000	3000
egrep -v '^\d+\t\d+$' | \

# Processing sequence: 10390000
egrep -v '^Processing sequence: \d+$' 
