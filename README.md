WT-finder is a python script for identifying the wild-type peptide counterpart of the variant peptides missed by the proteomics searches. This tool searches the precursor mass of wild-type peptide within the user-defined window corresponds to its variant peptide precursor mass. As the match found,  the program further looks for its daughter ion and performs the matching the theoretical derived b and y ions against experimental derived m/z values. 
The program provides two outputs.

1. An image file with the peptides b and y ion match
2. Tab delimited output

The tab delimited file provides the following information

id: accession number <br />
TD: -1 for deocy and 1 for wild-type peptides <br />
ScanNr: m/z scan number <br />
numberOfMatchingPeaks: number of b and y mathced with raw m/z peaks <br />
charge: charge of identified peptides (from raw file) <br />
Exp.Mz: exprimental m/z values <br />
RTinSec: retention time at peptides identified <br />
Theo.Mz: theoratical m/z values <br />
Delta_Mass_Error(PPM): <br />
length: lenght of peptides <br />
peptide: wild-type peptide <br />
protein: protein accession number <br />

