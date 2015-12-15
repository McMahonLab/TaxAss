Databases for 16S Taxonomy Assignment Workflow
===

This folder contains the taxonomy databases used for OTU assignment. The `.taxonomy` files contain tables with sequence IDs followed by taxonomic names. The `.fasta` files contain the nucleotide sequence for each sequence ID.

Custom Database
---

The custom database is the McMahon lab's curated freshwater ARB database. The most recent `.fasta` and `.taxonomy` files from this database are below, though they have not yet been incorporated in the [Freshwater Microbial Field Guide repo](https://github.com/mcmahon-uw/FWMFG).

* FWonly_10Dec2015_ready.fasta
*	FWonly_10Dec2015-taxonomy.nds

These files were sourced form an email from Trina to Robin, dated 2015-12-10 with the subject "shiniest new version of taxonomy training set".

General Database
---

The general database is the May 2013 release of GreenGenes, downloaded from the [GreeenGenes ftp](ftp://greengenes.microbio.me/greengenes_release/gg_13_5/) on 2015-08-03.
10-14-15

* gg_13_5_taxonomy
* gg_13_5.fasta

The mothur developers recommend a different version of this database, corresponding to the August 2013 release. This version is available from the [mothur website](http://mothur.org/wiki/Greengenes-formatted_databases).

****think that's a mistake on mothurs website- check on it!

Because these files are large, they are not stored in this repo.
