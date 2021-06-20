WT-finder is a python script for identifying the wild-type peptide counterpart of the variant peptides missed by the proteomics searches. This tool searches the precursor mass of wild-type peptide within the user-defined window corresponds to its variant peptide precursor mass. As the match found,  the program further looks for its daughter ion and performs the matching the theoretical derived b and y ions against experimental derived m/z values. 
The program provides two outputs.

1. An image file with the peptides b and y ion match
2. Tab delimited output

The tab delimited file provides the following information

<b>id</b>: accession number <br />
<b>TD</b>: -1 for deocy and 1 for wild-type peptides <br />
<b>ScanNr</b>: m/z scan number <br />
<b>numberOfMatchingPeaks</b>: number of b and y mathced with raw m/z peaks <br />
<b>charge</b>: charge of identified peptides (from raw file) <br />
<b>Exp.Mz</b>: exprimental m/z values <br />
<b>RTinSec</b>: retention time at peptides identified <br />
<b>Theo.Mz</b>: theoratical m/z values <br />
<b>Delta_Mass_Error(PPM)</b>: <br />
<b>length</b>: lenght of peptides <br />
<b>peptide</b>: wild-type peptide <br />
<b>protein</b>: protein accession number <br />

![alt text](https://github.com/sandeepkasaragod/WT-finder/blob/main/doc/QNQHEELQNVRK_00522_E02_P003811_B0M_A00_R1.png)
![alt text](https://github.com/sandeepkasaragod/WT-finder/blob/main/doc/Example_output.png)

