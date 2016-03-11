#!/bin/bash
# RRR 2/11/16

# sed is suuuper fucking annoying because it doesn't recognize tabs.  wtf?
# bash script to get around that and replace the first semicolon on each line with a tab.
# because I can just type a fucking tab here.


sed 's/;/	/' <custom.taxonomy.semicolons >custom.taxonomy

# and need a semicolon on the end of the line too

#fuck

exit 0