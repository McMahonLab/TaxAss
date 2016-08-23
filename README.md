16S Taxonomic Assignment Workflow
===
URL: https://mcmahonlab.wisc.edu/  
URL: https://github.com/McMahonLab/  

Overview
---

The TaxAss workflow assigns taxonomy to a fasta file of OTU sequences using both a small, custom taxonomy database and a large general database. TaxAss was developed specifically for the classification of freshwater OTUs against the McMahonLab's custom, curated Freshwater Training Set (The FreshTrain), but it should be generalizeable to use with other custom taxonomy databases. The manuscript to cite is in preparation- please check back!


How it works
---

TaxAss consists of a series of python, R, and bash scripts in addition to using mothur's classify.seqs() command.  The scripts are typed into the terminal window (mac or linux), and the commands used to call them are located in the `Clean Workflow.txt` file. (Sorry it's not in fancy Markdown yet.) The commands are listed at the top of the file, and then are listed below with detailed explanations of what they're doing.  The `Clean Workflow.txt` file also includes instructions about how to format your input files.  You start with a fasta file of your sequences or OTUs and the fasta and taxonomoy files for the two databases.  There's an optional step to explore your cutoff choices that also requires an OTU abundance table.


Where stuff is
---

The important stuff:  
1. Step by step directions are in the `Clean Workflow.txt` file.  
2. Scripts are in the `tax-scripts` folder.  
3. The FreshTrain taxonomy files are in the `database` folder.  
4. You have to download Greengenes or SILVA yourself.

The other stuff you can ignore:  
1. Scripts I use with the FreshTrain arb files are in the `arb-scripts` folder  
2. The scripts I used for making poster/paper figures are (in progress) and in `paper_figures` folder.  
3. A pdf of my ISME16 poster: `RohwerISME16poster.pdf`  


What if it doesn't work?
---

Please tell me!  If you find a bug (or if you get confused) I want to fix it/help you.  If you're all github savvy you submit an issue, or you can just send me an e-mail with TaxAss in the title: rrohwer@wisc.edu



