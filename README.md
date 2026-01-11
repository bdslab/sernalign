# SERNAlign - Structural sEquence RNA secondary structure Alignment

\@version 1.0

SERNAlign builds Structural Sequences starting from RNA secondary
structures with arbitrary pseudoknots and computes the *SERNA distance*
by aligning two *Structural Sequences*. Structural sequences are an
abstraction of the arc diagram of a secondary structure in which the
sequence of nucleotides is not considered and the topology of the arcs
is represented by a numerical sequence. For each RNA secondary structure
there is only one corresponding structural sequence.

If you use SERNAlign please cite:

- Tesei, L., Levi, F., Quadrini, M., Merelli, E.: Alignment of RNA
  Secondary Structures with Arbitrary Pseudoknots using Structural
  Sequences. Research Square preprint. (2024).
  [https://doi.org/10.21203/rs.3.rs-4831215/v1](https://doi.org/10.21203/rs.3.rs-4831215/v1)

## Accepted Input file formats for RNA secondary structures

- Extended Dot-Bracket Notation. See
  [https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/io/rna_structures.html](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/io/rna_structures.html)

- Extended Dot-Bracket Notation without sequence (only structure)

- Arc Annotated Sequence, similar to the Extended Dot-Bracket Notation
  format, in which the weak bonds are expressed as a list
  (i_1,j_1);(i_2,j_2); \... ;(i_m,j_m) where each index i_k, j_k belongs
  to the interval \[1,n\] (n is the length of the primary sequence) and
  i_k \< j_k + 1 for all k.

- Arc Annotated Sequence without sequence (only structure)

- BPSEQ format, with or without header information, see
  <https://www.ibi.vu.nl/programs/k2nwww/static/data_formats.html>

- CT format, with or without header information, see
  <https://www.ibi.vu.nl/programs/k2nwww/static/data_formats.html>

All lines starting with \# are considered comments. File format is
automatically detected from the text file, any file extension is
accepted.

# Installation

Download the latest release of SERNAlign from
[https://github.com/bdslab/sernalign/releases](https://github.com/bdslab/sernalign/releases)

The release contains the following files for version 1.0, 1.0 in the examples below:

- SERNAlign-v1.0.jar --- executable jar of the basic SERNAlign
  comparison

- SERNAlignWorkbench-v1.0.jar --- executable jar for the SERNAlign
  workbench comparator

- README.v1.0.md --- this file

- Source code in .zip and .tar.gz

Example structures and test structures are part of the source code
resources.

The executable jar files run on every Linux, Windows and Mac OS platform
in which a Java SE Runtime Environment at least version 8 is installed.

For information and installing the Java Runtime Environment see
<http://www.oracle.com/technetwork/java/javase/downloads/index.html>

## Using SERNAlign

Open a terminal window of your operating system and use the change
directory (cd) command to move to a folder in which the executable
jar(s) and the configuration file(s) were placed. To launch the basic
SERNAlign comparator digit:

`> java -jar SERNAlign-v1.0.jar <options>`

The following can be used:

```
  -a,--align <input-file1 input-file2>   Align two given structures
                                        producing an alignment and
                                        distance
 -d,--outdist                           Output only distance, no alignment
                                        (works only with option -a)
 -h,--help                              Show usage information
 -i,--info                              Show license and other info
 -n,--no-constraints                    Do not use constraints for the
                                        alignment (works only with option
                                        -a)
 -o,--out <output-file>                 Output result on the given file
                                        instead of standard output
 -s,--struct <input-file>               Produce the structural sequence
                                        corresponding to the given
                                        structure
```
## Using SERNAlignWorkbench

Open a terminal window of your operating system and use the change
directory (cd) command to move to a folder in which the executable
jar(s) were placed. To launch the basic SERNAlignWorkbench comparator
digit:

`> java -jar SERNAlignWorkbench-v1.0.jar <options>`

The following can be used:

```
 -f,--input <input-folder>     Process the files in the given folder
 -h,--help                     Show usage information
 -i,--info                     Show license and other info
 -j,--json                     Also generate output in JSON format
 -n,--no-constraints           Do not use constraints on the alignment
 -o,--output <file-1 file-2>   Output structure descriptions on file-1 and
                               comparison results on file-2 instead of
                               generating the default output files
```

## SERNAlign usage examples

Subfolder examples, distributed with SERNAlign source code, contains
sample input files from the paper: Tesei, L., Levi, F., Quadrini, M.,
Merelli, E. "Alignment of RNA Secondary Structures with Arbitrary
Pseudoknots using Structural Sequences".

`> java -jar SERNAlig-v1.0.jar -s examples/structS1.aas.txt`

Print on the standard output the structural sequence corresponding to
the RNA secondary structure given in the Arc Annotated Sequence file
examples/structS1.aas.txt.

In this particular case, the output is

```
[ 1, 1, 5, 3 ]
```

`> java -jar SERNAlign-v1.0.jar -a examples/structS1.aas.txt
  examples/structS2.aas.txt -o examples/alignmentS1S2.txt`

Write on the file examples/alignmentS1S2.txt the alignment and the SERNA
distance of the structural sequences corresponding to the RNA secondary
structures given in the Arc Annotated Sequence files
examples/structS1.aas.txt and examples/structS2.aas.txt.

In this particular case the file contains:

```
(1, 1)(1, 1)(5, 5)(3, 5)(-, 2)

Distance = 2
```

The sequence

```
(1, 1)(1, 1)(5, 5)(3, 5)(-, 2)
```

corresponds to the optimal alignment

```
1 1 5 3 -
1 1 5 5 2
```

between the structural sequences 1 1 5 3 and 1 1 5 5 2 where the edit
operations are (1, 1) (match), (1, 1) (match), (5, 5)
(mismatch/substitution), (3, 5), and (-, 2) (insertion).

The SERNA distance, i.e., the total cost of these edit operations is 0
(match) + 0 (match) + 0 (match) + 1 (insertion) + 1 (substitution) = 2

## SERNAlignWorkbench.jar usage examples

`> java -jar SERNAlignWorkbench-v1.0.jar -f examples/Eukaryota23S`

Processes all the files in folder "Eukaryota23S". Each file is read as
an RNA secondary structure with arbitrary pseudoknots. Comma-separated
values files "SERNAlignProcessedStructures.csv" and
"SERNAlignComparisonResults.csv" are created in the folder
"Eukaryota23S". The former contains the description of all the
structures that were found and correctly processed. The latter contains,
for each pair of processed structures, the SERNA Distance between the
two structures and execution time information.

`>java -jar SERNAlignWorkbench-v1.0.jar -f Eukaryota23S -o
  stucts.csv cmpr.csv`

Processes all the files in folder Eukaryota23S as above but produces the
description of processed structures in file structs.csv and comparison
results in file cmpr.csv.

# Copyright and License

SERNAling Copyright (C) 2024 Luca Tesei, Francesca Levi, Michela
Quadrini, Emanuela Merelli - BioShape and Data Science Lab at the
University of Camerino, Italy - <http://www.emanuelamerelli.eu/bigdata/>

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later
version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

# Contact Information

Please report any issue to luca.tesei@unicam.it or to Luca Tesei, Polo
Informatico, via Madonna delle Carceri 9, 62032 Camerino (MC) Italy.
