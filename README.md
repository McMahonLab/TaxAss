What's TaxAss??
---
TaxAss is a **Taxonomy Assignment** workflow that lets you classify 16S datasets using two different taxonomy reference databases. TaxAss was developed so that the McMahon Lab's Freshwater Training Set (the FreshTrain) could be used for taxonomy assignment alongside comprehensive databases like Greengenes or Silva. We believe it should also work with other custom databases (so if that's you, get in touch!)

How do I TaxAss?
---

**Step-by-step directions:** [tax-scripts/TaxAss_Directions.html](https://htmlpreview.github.io/?https://github.com/McMahonLab/TaxAss/blob/master/tax-scripts/TaxAss_Directions.html)

TaxAss uses a series of python, R, and bash scripts in addition to using BLAST+ and mothur's classify.seqs() command.  The scripts are sourced from the terminal window (mac or linux). You'll need to download this repository (green "Clone or download" button, top right), and then just add the tax-scripts folder to your working diriectory.

Where's the stuff I need?
---

**The stuff you need:**  
1. Step by step directions are in `TaxAss_Directions.html` inside the `tax-scripts` folder (view without download [here](https://htmlpreview.github.io/?https://github.com/McMahonLab/TaxAss/blob/master/tax-scripts/TaxAss_Directions.html)).  
2. Scripts are in the `tax-scripts` folder.  
3. The FreshTrain taxonomy files are in the `FreshTrain-files` folder.  

**The stuff you can ignore:**  
1. Scripts to process & edit the FreshTrain arb files are in the `arb-scripts` folder  
2. Scripts for making poster/paper figures are in `figure-scripts` folder.  
3. A pdf of my ISME16 poster (that might be a good overview): `poster-ISME16.pdf`  

What if it doesn't work?
---

Please tell me!  If you find a bug (or if you get confused) I want to fix it/help you.  If you're all github savvy you can submit an issue, or you can just send me an e-mail with TaxAss in the title: robin.rohwer@gmail.com

Who made TaxAss?
---
The McMahon lab studies lake microbial ecology (and some of it also studies wastewater). We are at the Univerisity of Wisconsin-Madison.  
Lab Website: https://mcmahonlab.wisc.edu/  
Lab Twitter: @mcmahonlab 
