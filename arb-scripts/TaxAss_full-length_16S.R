# RRR 1/17/18
# To satisfy reviewer #2, and also I think Pat Schloss bioRxiv review:
# Concern was: in blue/red improvement plot, "truth" is TaxAss to support using TaxAss over GG alone.
# The issue is that I'm testing TaxAss, so it shouldn't be used as the reference even though
# this example is just meant to demonstrate obvious errors in the GG-only method
# Approach for a better "truth":
# Ryan aligned some new full-length 16S sequences to the FreshTrain
# Take their assignment from using full length and arb as the "truth"
# I'll chop them into a variable region and then assign them using TaxAss
# This should validate that TaxAss works, and address concerns about overclassifying missing references

