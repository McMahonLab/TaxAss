16S Taxonomic Assignment Workflow
===
URL: https://mcmahonlab.wisc.edu/  
URL: https://github.com/McMahonLab/  

Overview
---

The TaxAss workflow assigns taxonomy to a fasta file of OTU sequences using both a small, custom taxonomy database and a large general database. TaxAss was developed specifically for the classification of freshwater 16S rRNA gene amplicon data using the McMahonLab's custom, curated Freshwater Training Set (The FreshTrain). However, TaxAss should be generalizeable to use with other custom taxonomy databases. The manuscript to cite is in preparation- please check back!


How it works
---

Directions are here: https://htmlpreview.github.io/?https://github.com/McMahonLab/TaxAss/blob/master/tax-scripts/TaxAss_Directions.html

TaxAss uses a series of python, R, and bash scripts in addition to using mothur's classify.seqs() command.  The scripts are typed into the terminal window (mac or linux), and the list of commands used to call them are located in the `workflow.txt` file. (Sorry it's not in fancy Markdown yet.) All the commands are listed at the top of the file, and then are listed below with detailed explanations of what they're doing.  The `workflow.txt` file also includes instructions about how to format your input files.  You start with a fasta file of your unique sequences or representative OTUs along with the fasta and taxonomoy files for the two databases.  There's an optional step to explore your cutoff choice that also requires an OTU abundance table.


Where stuff is
---

**The stuff you need:**  
1. Step by step directions are in `TaxAss_Directions.html` inside the `tax-scripts` folder (viewable without download from [here](https://htmlpreview.github.io/?https://github.com/McMahonLab/TaxAss/blob/master/tax-scripts/TaxAss_Directions.html)).  
2. Scripts are in the `tax-scripts` folder.  
3. The FreshTrain taxonomy files are in the `FreshTrain-files` folder.  

**The stuff you can ignore:**  
1. Scripts to process & edit the FreshTrain arb files are in the `arb-scripts` folder  
2. Scripts for making poster/paper figures are in `figure-scripts` folder.  
3. A pdf of my ISME16 poster (that might be a good overview): `poster-ISME16.pdf`  


What if it doesn't work?
---

Please tell me!  If you find a bug (or if you get confused) I want to fix it/help you.  If you're all github savvy you can submit an issue, or you can just send me an e-mail with TaxAss in the title: rrohwer@wisc.edu



