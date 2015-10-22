Scripts for 16S Taxonomy Assignment Workflow
===

This folder contains the scripts used for OTU assignment.


find_seqIDs_with_pident.R
---

Script that analyses the blast output and generates a list of seqIDs for  use in each classification database.  This was pulled from Robin's scripts file in her TAGs RProject. The name was 8-24-15_Filter_BLAST_for_tax_assignment.R


find_classification_disagreements.R
---

Script that searches for disagreements in name assignment between two taxonomy files. Use this to see if the chosen cutoffs are resulting in incorrect OTU assignment by comparing the OTU assignment with GG alone assignments to the OTU asignment of the chosen workflow. This was pulled from Robin's scripts file in her TAGs RProject. The name was 8-27-15_compare_FWblast_to_GG_classifications.R


fetch_fastas_with_seqIDs.py
---

Script that takes the list of seqIDs from the output of find_seqIDs_with_pident.R and generates a fasta file containing all those sequences to classify with mothur.


fetch_seqIDs_blast_removed.py
---

Script that takes the list of sequence IDs that BLAST returns and finds all the seqIDs BLAST didn't report because they were not good enough hits. Then it creates a .fasta file of these sequences that can be concatenated to the .fasta file of BLAST seqIDs under the pident cutoff.
