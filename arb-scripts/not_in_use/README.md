## README- defunct Arb Scripts


name													                            | Why it sucks  
----------------------------------------------------------|-----------------------
`find_typos.R`											                      | Replaced by<br>`tax-scripts/force_consistent_nomenclature_on_database.R`
`processSilvaFromMothur.py`								                | I think this was an old version when we added<br>`p__` etc to silva to match gg/freshtrain,<br>but now we instead remove those from FreshTrain<br>`arb-scripts/remove_gg_taxa_level_prefixes.sh`.  
`reformat_silva_failed_attempt.sh`						            | self explanatory in the name, lol. pretty sure done now by <br>`tax-scripts/force_consistent_nomenclature_on_database.R`  
`remove_non-monophyletic_refs_for_SILVA-compatibility.R`  | this is also an incomplete script-<br>started it but then just did it manually.  
`replace_blanks_with_unclassified.R`					            | this is now done by<br>`tax-scripts/force_consistent_nomenclature_on_database.R`    

