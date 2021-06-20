WT-finder is a python script for identifying the wild-type peptide counterpart of the variant peptides missed by the proteomics searches. This tool searches the precursor mass of wild-type peptide within the user-defined window corresponds to its variant peptide precursor mass. As the match found, the program further looks for its daughter ion and performs the matching the theoretical derived b and y ions against experimental derived m/z values. 
The program provides two outputs.

1. An image file with the peptides b and y ion match
2. Tab delimited output

The tab delimited file provides the following information
id,TD,ScanNr,numberOfMatchingPeaks,charge,Exp.Mz,RTinSec,Theo.Mz,Delta_Mass_Error(PPM),length,peptide and protein accession

