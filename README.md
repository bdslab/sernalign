# SERNAlign - Structural sEquence RNA secondary structure Alignment 

@version 1.0

SERNAlign builds Structural Sequences starting from RNA secondary structures 
with arbitrary pseudoknots and computes SERNA Distance by aligning two 
Structural Sequences. Structural sequences are an abstraction of the arc diagram 
of a secondary structure in which the sequence of nucleotides is not considered 
and the topology of the arcs is represented by a numerical sequence. For each
RNA secondary structure there is only ony corresponding structural sequence. 

If you use SERNAlign please cite:

- Tesei, L., Levi, F., Quadrini, M., Merelli, E. "Alignment of RNA Secondary 
Structures with Arbitrary Pseudoknots using Structural Sequences". Submitted
for publication. 2024


## Accepted Input file formats for RNA secondary structures 
* Extended Dot-Bracket Notation. See 
<https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/rna_structure_notations.html>
* Extended Dot-Bracket Notation without sequence (only structure)
* Arc Annotated Sequence, similar to the Extended Dot-Bracket Notation format, 
in which the weak bonds are expressed as a list `(i_1,j_1);(i_2,j_2); ... ;(i_m,j_m)` 
where each index i_k, j_k belongs to the interval [1,n] (n is the length
of the primary sequence) and i_k < j_k + 1 for all k.
* Arc Annotated Sequence without sequence (only structure)
* BPSEQ format, with or without header information, see <https://www.ibi.vu.nl/programs/k2nwww/static/data_formats.html>
* CT format, with or without header information, see <https://www.ibi.vu.nl/programs/k2nwww/static/data_formats.html>

All lines starting with \# are considered comments. File format is automatically
detected from the text file, any extension is accepted. 

# Installation

Download the zip file of last version of SERNAlign from folder `download` at <https://github.com/bdslab/sernalign/>

Direct link: <https://github.com/bdslab/sernalign/tree/master/download>

and put it in any position of your drive. 

Unzip the file with the facilities of your operating system. The folder 
`SERNAlign-<VersionNumber>` is created containing the following files:

- SERNAlign.jar --- executable jar of the basic SERNAlign comparison
- SERNAlignWorkbench.jar --- executable jar for the SERNAlign workbench comparator
- examples --- folder containing sample input and output files
- INSTALL.txt --- information on SERNAlign installation
- README.md --- SERNAlign description and usage information
- COPYING.txt --- copyright information
- LICENSE --- full GNU GPL Version 3 License
- CHANGELOG.txt --- information about the evolution of SERNAlign versions

The executable jar files runs on every Linux, Windows and Mac OS platform
in which a Java SE Runtime Environment 8 is installed. 

For information and installing the Java Runtime Environment see
<http://www.oracle.com/technetwork/java/javase/downloads/index.html>

## Using SERNAlign

Open a terminal window of your operating system and use the change directory 
(cd) command to move to a folder in which the executable jar(s) and the
configuration file(s) were placed. To launch the basic SERNAlign comparator 
digit:

`> java -jar SERNAlign.jar <options>`

The following <options> can be used:

	-a,--align <input-file1 input-file2>   Align two given structures
	                                       producing an alignment and
	                                       distance
	-d,--outdist                           Output only distance, no alignment
                                           (works only with option -a)
	-h,--help                              Show usage information
	-i,--info                              Show license and other info
	-o,--out <output-file>                 Output result on the given file
                                           instead of standard output
	-s,--struct <input-file>               Produce the structural sequence
                                           corresponding to the given
                                           structure
    -n, --no-constraints                   Do not use constraints in the alignment
                                           (works only with option -a)
## Using SERNAlignWorkbench

Open a terminal window of your operating system and use the change directory 
(cd) command to move to a folder in which the executable jar(s) were placed. 
To launch the basic SERNAlignWorkbench comparator digit:

`> java -jar SERNAlignWorkbench.jar <options>`

The following <options> can be used:

	-f,--input <input-folder>     Process the files in the given folder
	-h,--help                     Show usage information
	-i,--info                     Show license and other info
	-o,--output <file-1 file-2>   Output structure descriptions on file-1 and
                                  comparison results on file-2 instead of
                                  generating the default ouput files
    -n, --no-constraints          Do not use constraints in the alignment
                                  (works only with option -a)
## SERNAlign.jar usage examples

Subfolder `examples`, distributed with SERNAlign, contains sample input files 
from the paper: Tesei, L., Levi, F., Quadrini, M., Merelli, E. "Alignment of RNA Secondary 
Structures with Arbitrary Pseudoknots using Structural Sequences". Submitted
for publication. 2024

* `> java -jar SERNAlign.jar -s examples/structS1.aas.txt`

Print on the standard output the structural sequence corresponding
to the RNA secondary structure given in the Arc Annotated Sequence file
`examples/structS1.aas.txt`.

In this particular case, the output is 

	[ 1, 1, 5, 3 ]


* `> java -jar SERNAlign.jar -a examples/structS1.aas.txt examples/structS2.aas.txt -o examples/alignmentS1S2.txt`

Write on the file `examples/alignmentS1S2.txt` the alignment and the SERNA 
distance of the structural sequences corresponding to the RNA secondary 
structures given in the Arc Annotated Sequence files `examples/structS1.aas.txt` and 
`examples/structS2.aas.txt`.

In this particular case the file contains: 

	(1, 1)(1, 1)(5, 5)(-, 5)(3, 2)

	Distance = 2

The sequence 

	(1, 1)(1, 1)(5, 5)(-, 5)(3, 2)

corresponds to the optimal alignmet 

	1 1 5 - 3
	1 1 5 5 2
	
between the structural sequences `1 1 5 3` and `1 1 5 5 2` where the edit operations
are `(1, 1)` (match), `(1, 1)` (match), `(5, 5)` (match), `(-, 5)` (insertion) 
and `(3, 2)` (mismatch/substitution). 

The SERNA distance, i.e., the total cost of these edit operations 
is 0 (match) + 0 (match) + 0 (match) + 1 (insertion) + 1 (substitution) = 2

## SERNAlignWorkbench.jar usage examples

* `> java -jar SERNAlignWorkbench.jar -f TestWorkBench1`

Processes all the files in folder "`TestWorkBench1`". Each file is read as
an RNA secondary structure with arbitrary pseudoknots. Comma-separated values 
files "`SERNAlignProcessedStructures.csv`" and "`SERNAlignComparisonResults.csv`"
are created in the folder "`TestWorkBench1`". The former contains the
description of all the structures that were found and correctly processed.
The latter contains, for each pair of processed structures, the SERNA
Distance between the two structures and execution time information.

* `>java -jar SERNAlignWorkbench.jar -f TestWorkBench1 -o stucts.csv
cmpr.csv`

Processes all the files in folder `TestWorkBench1` as above but produce
the description of processed structures in file `structs.csv` and
comparison results in file `cmpr.csv`. 

See folder `examples` for some sample input folders containing structures 
coming from public databases.

# Copyright and License

SERNAling Copyright (C) 2020  Luca Tesei, Francesca Levi, Michela Quadrini, 
Emanuela Merelli - BioShape and Data Science Lab at the University of 
Camerino, Italy - <http://www.emanuelamerelli.eu/bigdata/>

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later
version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.

# Contact Information

Please report any issue to luca.tesei@unicam.it or to Luca Tesei, Polo
Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.


