16S Taxonomic Assignment Workflow
===
Copyright (c) 2015, Katherine D. McMahon and Lackeys  
URL: https://mcmahonlab.wisc.edu/  
URL: https://github.com/McMahonLab/  
All rights reserved.

Overview
---

This workflow assigns taxonomy to a fasta file of OTU sequences using both a
small, custom taxonomy database and a large general database. The workflow was developed specifically for the classification of freshwater OTUs against the McMahonLab's [curated freshwater database](https://github.com/mcmahon-uw/FWMFG). Your mileage may vary.

Background
---

Our lab's 16S taxonomy pipeline uses the [RDP classifier](https://github.com/rdpstaff/classifier) to classify OTUs. as implemented in [mothur](http://www.mothur.org/) (or [qiime](http://qiime.org/)).

The original process was straightforward. Sequences were classified using our FW database, and all "lineage"-level assignments with greater than 70% bootstrap support were kept. The remaining OTUs were classified using the GreenGenes database, with 60% bootstrap support.

However, this workflow was flawed due to a misunderstanding of the bootstrap support values. The RDP classifier is a Bayesian classifier, and was developed using the (large) [RDP database](https://rdp.cme.msu.edu/). The classifier uses the reference database to calculate the probability that a novel sequence belongs to a particular OTU defined in the reference database. The OTU is classified as a member of the genus giving the highest probability score, __regardless of the value of that probability__. The process is repeated 100 times to give a bootstrap support score, which is reported by the classifier. In other words, the classifier will classify any sequence you give it: if the "true" classification is absent from the classifier, the sequence will still be classified, albeit possibly with low probability (unreported) and low bootstrap support (reported).

This procedure is acceptable for the GreenGenes or RDP databases. The databases are large: if you attempt to classify something that the database doesn't have,	it will likely be classified differently each time and the bootstrap support score will be low (e.g., the sequence will be "unclassified"). However, our database is small. If you attempt to classify a novel sequence, the classifier will report the most probable sequence, even if the probability is low. Because there are fewer sequences, the bootstrap support value is also likely to be higher. Because the classifier does not report the probabilities, you have no __a priori__ way of knowing if a sequence has been mis-classified. As a simple example, if you give our FW database an archaeal sequence, 100% of the time the sequence will be classified as bacterial, because the database only has bacteria in it.

As a consequence, we were "forcing" sequences to be classified as belonging to defined FW taxa, even if they did not match well. We may have "missed" other important (non-FW) taxa, because they were not represented in the database. Finally, we may have muddled the relationships of our key taxa by adding in unrelated ones.

Possible Solutions
---
Robin Rohwer identified and evaluated a number of possible solutions, as described below.

1. Classify in GreenGenes first, then reclassify just those phyla with known freshwater representatives. This approach was adopted by Jason Flowers. However, phyla with known freshwater representatives also contain non-FW taxa whih are not in our database.

2. Classify in GreenGenes first, then reclassify just those orders unique to freshwater. However, GreenGenes has insufficient sequences from those orders to give a confident classification, so many FW taxa would get "unclassified." As a result, we would miss many "true" FW sequences.

3. Combine the two taxonomy databases and classify all samples at once. However, if some OTUs are present in both databases (classified differently, one FW and one the Linnean equivalent), novel sequences will be split during classification (50% FW, 50% Linnean), and remain unclassified.

4. Classify our sequences with GreenGenes, remove those sequences, and the combine the databases to classify all at once. However, the taxonomy databases are highly curated and adding something new is not that simple. Further, the structure of the two databases may differ at finer taxonomic resolution. Thus, the GreenGenes classification might not be a good match if the FW sequence doesn't exist in GreenGenes, so we would be removing an unrelated sequence that we shouldn't. There is no good way to differentiate between scenarios where the FW sequence doesn't exist in GreenGenes, or exists and is just classified to the class level (FW taxonomy differs from Linnean at order and narrower).

5. Classify the sequences in the FW database with a different method, such as BLAST, since it's a small database anyway. However, the taxonomy assignment algorithm takes into consideration phylogeny, which is a better classification than that obtained by just using sequence similarity.

Proposed Solution
---

Our propose solution is as follows, and as documented in further detail in this workflow.

1. First classify novel sequences using the FW database, but define a cutoff independent of RDP classifer bootstrap support values.
2. Use a BLAST cutoff to identify highly similar sequences to those in our FW sequence database.
3. Classify those sequences using the FW database.
4. Classify the remainining sequences with GreenGenes.

Prerequisites
---
Prior to using this workflow, please make sure the following software are installed and working:

* [mothur](http://mothur.org/wiki/Main_Page)
* [BLAST](http://www.ncbi.nlm.nih.gov/books/NBK52640/)
* [Python](https://www.python.org/)
* [R](https://www.r-project.org/)
* A text editor, such as [TextWrangler](http://www.barebones.com/products/textwrangler/), [Atom](https://atom.io/), or [emacs](https://www.gnu.org/software/emacs/)

Workflow Summary
---
0. Formatting
  * Format all nucleotide sequence files such that they follow the example given in `scripts\FWonly_7Sept2015.fasta`
  * Format the taxonomic databases such that follow the example given in `scripts/FWonly_7Sept2015.taxonomy`.

1. BLAST Database Creation

    `makeblastdb -dbtype nucl -in custom.fasta -input_type fasta -parse_seqids -out custom.db`

Detailed Workflow Instructions and Notes
---

0. Formatting of nucleotide sequences and taxonomic databases  

    Reformat `fasta` files of nucleotide sequences to have no whitespace in the seqID line. BLAST does not parse whitespace well, and will consider any text after a space to be a comment. Having BLAST IDs that don't match the full comment line of the `fasta` file break the script `fetch_fastas_with_seqIDs.py`. **This file of sequences to be classified will be called `otus.fasta`.**


    Reformat your taxonomy files to be compatible with `mothur`.  `mothur` requires each line to be of the form:

      `seqID<tab>k_kingdom;p_phylum;...;s_species;`

    where taxonomic levels beneath `kingdom` are optional. The file __cannot__ contain whitespace, and the semi-colon at the end is required.

    In this workflow, we will refer to the following taxonomy files:

    * `custom` - the curated database you wish to use before primary classification.
    * `general` - the large database you with to use for classifying the remaining sequences.

    Each database will have two files:

    * `.fasta` - fasta nucleotide file of sequences and IDs
    * `.taxonomy` - IDs with taxonomic names


1. Creation of BLAST database from curated taxonomy

    Use the command `makeblastdb` to create a BLAST database out of the FW taxonomy fasta files. We need to create a database because:  
    	1. BLAST will run faster  
    	2. Having a database is necessary for some of the output formats

    Full command (typed in terminal):  

      `makeblastdb -dbtype nucl -in custom.fasta -input_type fasta -parse_seqids -out custom.db`

    What the filenames are:

    | File  | Description  |
    |---|---|
    | `custom.fasta` | Path to the small custom taxonomy file you want to use to assign for primary classification (e.g., the freshwater taxonomy database). __Note__: if the file path contains spaces, the path needs to be double-quoted: `' "/path/double quoted" '`.  |
    | `custom.db` | Name of the BLAST DB to be created. The command will generate six files with different final extensions, but all begin with this string.  

    What each flag does:

    | Flag | Description |
    |------|-------------|
    | `-dbtype nucl`	| A required argument, specifying the database input file contains nucleic acid sequencess (nucl) |
    | `-in custom.fasta` |	Path of the input file to be used for DB creation |
    | `-input_type fasta` | Specify the input file is in `fasta` format |
    | `-parse_seqids` | Include sequence IDs in the formatted DB, so that sequences can later be retrieved. Also necessary for using blast_formatter later. |
    | `-out custom.db` | Specify path for the database files to be created |

    These six files are created:  
    * custom.db.nhr  
  	* custom.db.nog  
  	* custom.db.nsd  
  	* custom.db.nsi  
  	* custom.db.nsq  
    But the string `custom.db` is all that's required to invoke the DB later.
