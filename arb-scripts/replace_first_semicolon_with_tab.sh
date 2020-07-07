#!/bin/bash
# RRR 2/11/16

# sed is suuuper fucking annoying because it doesn't recognize tabs.  wtf?
# bash script to get around that and replace the first semicolon on each line with a tab.
# because I can just type a fucking tab here.

# syntax: $1 is the first argument supplied (the file with all semicolons)
# and $2 is the second argument supplied (to mothur-formatted created file)
# ./replace_first_semicolon_with_tab.sh FT_semicol_noprefix.tax FT_semicol_noprefix_mothur-fmt.tax


sed 's/;/	/' <$1 >$2


exit 0