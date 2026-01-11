/**
 * SERNAlign - Structural sEquence RNA secondary structure Alignment
 * 
 * Copyright (C) 2024 Luca Tesei, Francesca Levi, Michela Quadrini,
 * Emanuela Merelli - BioShape and Data Science Lab at the University of 
 * Camerino, Italy - http://www.emanuelamerelli.eu/bigdata/
 *  
 * This file is part of SERNAlign.
 * 
 * SERNAlign is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * 
 * SERNAlign is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SERNAlign. If not, see <http://www.gnu.org/licenses/>.
 */
package it.unicam.cs.bdslab.sernalign;

/**
 * Collector of strings for messages to send on the console.
 * 
 * @author Luca Tesei
 *
 */
public interface CommandLineMessages {
    public static String VERSION = "1.0";

    public static String LAUNCH_COMMAND = "java -jar SERNAlign.jar";
    public static String HEADER = "\n\nSERNAling Comparator version "
	    + VERSION
	    + " - Determine the Structural Sequence associated to an RNA "
	    + "Secondary Structures with arbitrary pseudoknots and "
	    + "compute the SERNA Distance by alignign two Structural "
	    + "Sequences.\n\nDefault input file format is Extended Dot-Bracket"
	    + " Notation (EDBN). Alternatively, accepted formats are BPSEQ, CT and "
	    + "Arc Annotated Sequence (AAS), similar to the "
	    + "EDBN format in which the"
	    + " weak bonds are expressed as a list (i_1,j_1);(i_2,j_2); ... ;"
	    + "(i_m,j_m) where each index i_k, j_k belongs to the interval "
	    + "[1,n] (where n is the length of the primary sequence) and "
	    + "i_k < j_k + 1 for all k. In intput file formats EDBN and AAS the "
	    + "sequence of nucleotides is optional.\n\n"
	    + "Option -o for specifying output file is optional, if "
	    + "not specified the result is printed on the standard output.\n\n";
    public static String USAGE_EXAMPLES = "Usage examples:\n\n>"
	    + LAUNCH_COMMAND
	    + " -s structS1.aas.txt\n\nPrint on the standard output "
	    + "the structural sequence "
	    + "corresponding to the RNA secondary structure given in the "
	    + "Arc Annotated Sequence file structS1.aas.txt\n\n>" + LAUNCH_COMMAND
	    + " -a structS1.aas.txt structS2.aas.txt\n\nPrint on the standard "
	    + "output one optimal correct alignment of the structural "
	    + "sequences corresponding to the two given RNA secondary structures "
	    + "structS1.aas.txt and structS2.aas.txt. The alignment is a sequence of"
	    + "edit operations of the form (i, i) - match, (i, j) - mismatch,"
	    + " (-, i) - insertion, (i, -) - deletion\n\n";

    public static String LAUNCH_COMMAND_WB = "java -jar SERNAlignWorkbench.jar";
    // TODO modificare stringhe _WB
    public static String HEADER_WB = "\n\nSERNAling Workbench Comparator version "
	    + VERSION
	    + " - Compare all the RNA secondary structures files in a given "
	    + "input folder by computing the SERNA Distance between all possible "
	    + "pairs. All the files are expected to be in one of the recognised "
	    + "formats (see SERNAlign help). The output is given as two "
	    + "comma-separated values files describing the processed "
	    + "structures and containing the SERNA Distance calculated "
	    + "for each pair of processed structures. By default the "
	    + "output files are put in the input folder. Use option -o"
	    + " file-1 file-2 to specify different output files.\n\n";
    public static String USAGE_EXAMPLES_WB = "Usage examples:\n\n>"
	    + LAUNCH_COMMAND_WB
	    + " -f Eukaryota23S\n\nProcess all the files in folder "
	    + "\"Eukaryota23S\". Each file is read as an RNA secondary "
	    + "structure with arbitrary pseudoknots in one supported format. "
	    + "Comma-separated values files \"SERNAlignProcessedStructures.csv\" "
	    + "and \"SERNAlignComparisonResults.csv\" are created in the folder "
	    + "\"Eukaryota23S\". The former contains the description of all "
	    + "the structures that were found and correctly processed. The "
	    + "latter contains, for each pair of processed structures, the "
	    + "SERNA Distance between the two structures and execution time "
	    + "information.\n\n" + LAUNCH_COMMAND_WB
	    + " -f Eukaryota23S -o "
	    + "stucts.csv cmpr.csv \n\nProcesses all the files "
	    + "in folder \"Eukaryota23S\" as above but produce the description"
	    + " of processed structures in file \"structs.csv\" and comparison "
	    + "results in file \"cmpr.csv\".\n\n";

    public static String COPYRIGHT = "*********************************************\nSERNAling Copyright (C) 2024 Luca Tesei, Francesca Levi, Michela Quadrini, "
	    + "Emanuela Merelli - BioShape and Data Science Lab at the University of Camerino,"
	    + " Italy - http://www.emanuelamerelli.eu/bigdata/\n\n";
    public static String SHORT_NOTICE = "This program comes with ABSOLUTELY NO WARRANTY; for details use "
	    + "option '-i'. This is free software, and you are welcome to redistribute it "
	    + "under certain conditions; use option '-i' for more details.\n\n";
    public static String LONG_NOTICE = "This program is free software: you can redistribute it and/or modify "
	    + "it under the terms of the GNU General Public License as published by "
	    + "the Free Software Foundation, either version 3 of the License, or "
	    + "any later version.\n\n"
	    + "This program is distributed in the hope that it will be useful, "
	    + "but WITHOUT ANY WARRANTY; without even the implied warranty of "
	    + "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the "
	    + "GNU General Public License for more details.\n\n"
	    + "You should have received a copy of the GNU General Public License "
	    + "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n";
    public static String REPORT_TO = "Please report any issue to luca.tesei@unicam.it or to Luca Tesei, "
	    + "Polo Informatico, via Madonna delle Carceri 7, 62032 Camerino (MC) Italy.";

}
