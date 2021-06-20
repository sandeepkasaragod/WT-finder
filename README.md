WT-finder is a python script for identifying the wild-type peptide counterpart of the variant peptides missed by the proteomics searches. This tool searches the precursor mass of wild-type peptide within the user-defined window corresponds to its variant peptide precursor mass. As the match found,  the program further looks for its daughter ion and performs the matching the theoretical derived b and y ions against experimental derived m/z values. 
The program provides two outputs.

1. An image file with the peptides b and y ion match
2. Tab delimited output

The tab delimited file provides the following information

id: accession number
TD: -1 for deocy and 1 for wild-type peptides, 
ScanNr: m/z scan number
numberOfMatchingPeaks: number of b and y mathced with raw m/z peaks 
charge: charge of identified peptides (from raw file)
Exp.Mz: exprimental m/z values
RTinSec: retention time at peptides identified
Theo.Mz: theoratical m/z values
Delta_Mass_Error(PPM):
length: lenght of peptides
peptide: wild-type peptide
protein: protein accession number

