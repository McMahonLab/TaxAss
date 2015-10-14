10-14-15 RRR

These are all the scripts I wrote for my taxonomy assignment workflow.


8-24-15_Filter_BLAST_for_tax_assignment.R

	Robin's R script that analyses the blast output and generates a list of seqIDs for 
	use in each classification database.  This was pulled from my scripts file in my TAGs
	RProject.  In the workflow I change the name to find_seqIDs_with_pident.R
	

8-27-15_compare_FWblast_to_GG_classifications.R

	Robin's R script that searches for disagreements in name assignment between two taxonomy files.
	Use this to see if the chosen cutoffs are resulting in incorrect OTU assignment by comparing the
	OTU assignment with GG alone assignments to the OTU asignment of the chosen workflow.
	For the workflow I change its file name to find_classification_disagreements.R


fetch_fastas_with_seqIDs.py

	Robin's python script that takes the list of seqIDs from the 8-24-15_Filter_BLAST_for_tax_assignment.R 
	R script and generates a fasta file containing all those sequences to classify with mothur.


fetch_seqIDs_blast_removed.py

	Robin's python script that takes the list of sequence IDs that BLAST returns and finds
	all the seqIDs BLAST didn't report because they were not good enough hits.
	Then it creates a .fasta file of these sequences that can be concatenated to the .fasta
	file of BLAST seqIDs under the pident cutoff.

