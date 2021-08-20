The FreshTrain
===

The zipped files in this folder contain the McMahon Lab's custom taxonomy database for lake epilimia, fondly referred to as *The FreshTrain* (for Freshwater Training Set). **There are 6 files inside each zipped folder:**  

1. The FreshTrain `.taxonomy` file is a table with sequence IDs followed by taxonomic names.   
2. The FreshTrain `.fasta` file lists the nucleotide sequence for each sequence ID.  
3. The silva or greengenes `.taxonomy` file is formatted for taxass, for your convenience.
4. The silva or greengenes `.fasta` file is formatted for taxass, for your convenience.
5. The `README` file lists the starting files, commands, and output for creating the version.  
6. The `Compare_FreshTrain_and_General_Names.xlsx` is a reference of taxon name comparisons.  

There are an additional 3 files for Silva 138:  

7. The silva file `cyano_edits.taxonomy` can alternately be paired with the silva fasta file.  
8. The `README-cyano_edits` file details the edits to "clean up" but not "correct" silva's Cyanobacteria genus names.  
9. The `Compare_Edited_Cyano_and_Silva138_Names.xlsx` is a reference of _Cyanobacteria_ name comparisons.


<br> 
  
zipped file name        | description
------------------------|----------------------------------------------------------------
FreshTrain15Jun2020greengenes13_8.zip | Greengenes version  
FreshTrain15Jun2020silva128.zip | Silva 128 version  
FreshTrain15Jun2020silva132.zip | Silva 132 version  
**FreshTrain15Jun2020silva138.zip** | **Silva 138 version<br>This is the recommended, most-recent version!**

The different formats match the FreshTrain's coarse-level nomenclature to the nomenclature in the comprehensive database of choice. The FreshTrain defines lineage-clade-tribe (~family-genus-species) level phylogenies, so the phylum, class, and order names are changed in the different FreshTrain versions to be consistent with the paired comprehensive database.  

The silva/greengenes files are downloaded from the mothur website with minor changes- mainly that missing names are uniformly changed to "unnamed.parentname" and that the fasta files are unaligned (which takes a while to run). Directions and explanations of this general database formatting are in Step 0 of the TaxAss directions.  

<br>
The citation for the original FreshTrain database and the arb version of it is:  

[Newton RJ, Jones SE, Eiler A, McMahon KD, Bertilsson S. 2011. A Guide to the Natural History of Freshwater Lake Bacteria. Microbiol Mol Biol Rev 75:14–49.](https://mmbr.asm.org/content/75/1/14.full)  The original arb files are available at [github.com/McMahonLab/FWMFG](https://github.com/McMahonLab/FWMFG), but you should contact us for the most recent version.  

<br>
The citation for these taxonomy assignment-compatible formats of the FreshTrain and the TaxAss method is:  

[Rohwer RR, Hamilton JJ, Newton RJ, McMahon KD. 2018. TaxAss: Leveraging a Custom Freshwater Database Achieves Fine-Scale Taxonomic Resolution. mSphere 3:e00327-18.](https://msphere.asm.org/content/3/5/e00327-18)  
