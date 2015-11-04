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
* [R](https://www.r-project.org/), version 3.2 or newer (required for Rscript in steps 4 and 11)
* A text editor, such as [TextWrangler](http://www.barebones.com/products/textwrangler/), [Atom](https://atom.io/), or [emacs](https://www.gnu.org/software/emacs/)

Workflow Summary
---
0. Formatting
  * Format all nucleotide sequence files such that they follow the example given in `scripts\FWonly_7Sept2015.fasta`
  * Format the taxonomic databases such that follow the example given in `scripts/FWonly_7Sept2015.taxonomy`.

1. BLAST Database Creation

    `makeblastdb -dbtype nucl -in custom.fasta -input_type fasta -parse_seqids -out custom.db`

2. Run BLAST

    `blastn -query otus.fasta -task megablast -db custom.db -out otus.custom.blast -outfmt 11 -max_target_seqs 1`

3. Reformat BLAST Database

    ` blast_formatter -archive otus.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table`

4. Filter BLAST Results (Run one script twice)

    `Rscript find_seqIDs_with_pident.R otus.custom.blast.table outputfile hitsFile cutoff TRUE`

    `Rscript find_seqIDs_with_pident.R otus.custom.blast.table outputfile hitsFile cutoff FALSE`

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

2. Run BLAST

    Use the command `blastn` to run a megablast that returns the best hit in the taxonomy database (subject) for each of your OTU sequences (queries). Megablast is optimized for finding very similar matches with sequences longer than 30 bp.  

    __Note__: This workflow was made for ~100 bp sequences and may need to be re-optimized for sequences of different lengths.

    Full Command (type in terminal):

    `blastn -query otus.fasta -task megablast -db custom.db -out otus.custom.blast -outfmt 11 -max_target_seqs 1`


    What the filenames are:

    | File  | Description  |
    |---|---|
    | otus.fasta | The query file. The `fasta` file of sequences you want classified. |
    | custom.db | The subject file. The BLAST database made from custom taxonomy sequences in Step 1. |
    |  otus.custom.blast	| The BLAST output file. Its format is unreadable by you, but will be reformatted using `blast_formatter`. |

    What each flag does:						

    | Flag | Description |
    |------|-------------|
    | -query otus.fasta |	Location of the query file. |
    | -task megablast | Optimized version of BLAST to detect high similarity hits. It is also the blastn default task. More info [here](http://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/BLAST/nucleotide_blast.html). |
    | -db custom.db | Name of the blast database created previously (without the additional file extensions).|
    | -out otus.custom.blast | Location of the BLAST output file to be created. |
    | -outfmt 11 | Format to be used with `blast_formatter`. ASN.1 format. |
    |	-max_target_seqs 1 | Only keep the best hit for each query sequence. ]

    __Note__.
      * We would like to find the longest matching sequence with a percent ID which still matches the cutoff. However, it's possible that the best hit is a shorter sequence with a higher percent ID.
      * This is something that could be examined to decide if `blastn` custom scoring should be used instead of `megablast`.
      * The R script in step 4 reports stats on alignment length vs. query length. If many alignments are shorter than the full sequence, consider evaluating beyond the first hit.


3. Reformat BLAST Results

    The `blast_formatter` function in reformats the ASN.1-formatted file and reformats it to other formats. We reformat to a custom table format to feed into `find_seqIDs_with_pident.R`. However, using `blast_formatter` you can look at your blast results from the previous step any way you'd like.

    Full Command (type in terminal):

    `blast_formatter -archive otus.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table`

      What the filenames are:

      | File  | Description  |
      |---|---|
      | otus.custom.blast	| BLAST result file generated in Step 2, in the ASN.1 blast file format. |
      | otus.custom.blast.table | BLAST results reformatted into a six-column table for `find_seqIDs_with_pident.R`.  The tab-delimited columns are: qseqid pident length qlen qstart qend.

      What each flag does:		

      | Flag | Description |
      |------|-------------|
      | -archive otus.custom.blast | Specify path to the blast result file you are reformatting. |
      | -outfmt "6 qseqid pident length qlen qstart qend" | `6` is tabular format without headers or other information between rows of data, with columns further describe below. |
      | -out otus.custom.blast.table | Specify path to the output file of formatted blast results. |

      Columns in re-formatted BLAST results.

      | Number | Name | Description |
      |--------|------|-------------|
      | 1 | qseqid | query (OTU) sequence ID |
      | 2 | pident | percent identity (# of matches / # "columns" in the sequence) |
      | 3 | length | length of alignment |
      | 4 | qlen | full length of query sequence |
      | 5 | qstart | index of beginning of alignment on query sequence |
      | 6 | qend | index of end of alignment on query sequence |

4. Filter BLAST results

      `find_seqIDs_with_pident.R` takes the formatted BLAST file and
      calculates a corrected pident value that corrects the value for the entire length of the query. The corrected pident is a worst case scenario that assumes any edge gaps are mismatches.

      Calculation:
      `corrected pident" = pident * length / (length - (qend - qstart) + qlen)`

      The R script is run twice, once to generate the sequence ID's above/equal to the "true pident" cutoff that are destined for taxonomy assignment in the small custom database, and once to generate the sequence ID's below the cutoff that are destined for taxonomy assignment in the large general database.

      Full Two Commands (type in terminal):

      `Rscript find_seqIDs_with_pident.R otus.custom.blast.table outputfile hitsFile cutoff TRUE`

      `Rscript find_seqIDs_with_pident.R otus.custom.blast.table outputfile hitsFile cutoff FALSE`

      __Note__: Rscript requires R v 3.2 or higher. The path to the R executable must be added to your `PATH` variable.

      The arguments must be in the correct order. Rscript sources the .R script using the arguments supplied after it in the terminal. Separate all arguments with a space.

      What each argument is:

      | Number | Name | Description |
      |--------|------|-------------|
      | 1 | script.R | Path to `find_seqIDs_with_pident.R` |
      | 2 | otus.custom.blast.table | Path to `otus.custom.blast.table`from Step 3 |
      |	3 | outputfile | Path the to file you are creating, the list of sequence ID's matching your criteria (T or F for meeting the cutoff). |
      |	4 | hitsfile | Path to a secondary file. This will (optinally) be used to create plot to assess the pident. |
      | 5 | Cutoff | Numeric representing the "corrected pident" to use for matches. |
      | 6 | Matches | TRUE or FALSE. TRUE: return seqID's >= cutoff, FALSE: return seqID's < cutoff |

      What each file is:

      | File | Description |
      |------|-------------|
      | find_seqIDs_with_pident.R	| Script for this step. |
      | otus.custom.blast.table	| The formatted blast output from step 3
      | outputfile | Path the to file you are creating, the list of sequence ID's matching your criteria (T or F for meeting the cutoff). These files are newline \n delimited seqIDs. |

      __Note__: You may need to choose a different "corrected pident" cutoff for your sequence data. We selected a pident that gave classifications consistent to the class level between the small and large databases.
