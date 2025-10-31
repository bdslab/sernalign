/**
 * SERNAlign - Structural sEquence RNA secondary structure Alignment
 * <p>
 * Copyright (C) 2023 Luca Tesei, Francesca Levi, Michela Quadrini,
 * Emanuela Merelli - BioShape and Data Science Lab at the University of
 * Camerino, Italy - http://www.emanuelamerelli.eu/bigdata/
 * <p>
 * This file is part of SERNAlign.
 * <p>
 * SERNAlign is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 * <p>
 * SERNAlign is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with SERNAlign. If not, see <http://www.gnu.org/licenses/>.
 */
package it.unicam.cs.bdslab.sernalign;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;

//import javax.swing.JFileChooser;
//import javax.swing.JOptionPane;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * This class contains a main that runs the SERNAlign comparison algorithm
 * among all the RNA secondary structures (with arbitrary pseudoknots) in a
 * given folder.
 *
 * Two comma-separated-values files are produced with the description of the
 * processed files and the distance among all the pairs of molecules.
 * Additional information about the size of the molecules and the execution
 * times is output as well.
 *
 * @author Luca Tesei
 *
 */
public class WorkbenchComparator {

    public static void main(String[] args) {
        // Use Apache Commons CLI 1.4
        // create Options object for Command Line Definition
        Options options = new Options();

        // define command line options
        Option o1 = Option.builder("f")
                .desc("Process the files in the given folder")
                .longOpt("input").hasArg().argName("input-folder").build();
        options.addOption(o1);
        Option o2 = Option.builder("o").desc(
                        "Output structure descriptions on file-1 and comparison results on file-2 instead of generating the default ouput files")
                .longOpt("output").hasArgs().numberOfArgs(2)
                .argName("file-1 file-2").build();
        options.addOption(o2);
        Option o3 = Option.builder("i").desc("Show license and other info")
                .longOpt("info").build();
        options.addOption(o3);
        Option o4 = Option.builder("h").desc("Show usage information")
                .longOpt("help").build();
        options.addOption(o4);
        Option o5 = Option.builder("j")
                .desc("Also generate output in JSON format")
                .longOpt("json")
                .build();
        options.addOption(o5);

        options.addOption(
                "n",
                "no-constraints",
                false,
                "Do not use constraints on the alignment");

	/*Option o6 = Option.builder("c").desc(
		"Check the presence of only standard Watson-Crick and wobble base pairing (disabled by default)")
		.longOpt("chkpair").build();
	options.addOption(o6);
	Option o7 = Option.builder("e")
		.desc("Show current values of edit scores used for alignment")
		.longOpt("showscores").build();
	options.addOption(o7);
	Option o8 = Option.builder("n").desc(
		"Use the specified configuration file instead of the default one")
		.longOpt("useconffile").hasArg().argName("conf-file").build();
	options.addOption(o8);*/

        // Parse command line
        HelpFormatter formatter = new HelpFormatter();
        CommandLineParser commandLineParser = new DefaultParser();
        CommandLine cmd = null;
        try {
            cmd = commandLineParser.parse(options, args);
        } catch (ParseException e) {
            // oops, something went wrong
            System.err.println("ERROR: Command Line parsing failed.  Reason: "
                    + e.getMessage() + "\n");
            formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND_WB,
                    CommandLineMessages.HEADER_WB, options,
                    CommandLineMessages.USAGE_EXAMPLES_WB
                            + CommandLineMessages.COPYRIGHT
                            + CommandLineMessages.SHORT_NOTICE
                            + CommandLineMessages.REPORT_TO,
                    true);
            System.exit(1);
        }
        boolean noConstraints = cmd.hasOption("n");
        // Manage Option h
        if (cmd.hasOption("h")) {
            formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND_WB,
                    CommandLineMessages.HEADER_WB, options,
                    CommandLineMessages.USAGE_EXAMPLES_WB
                            + CommandLineMessages.COPYRIGHT
                            + CommandLineMessages.SHORT_NOTICE
                            + CommandLineMessages.REPORT_TO,
                    true);
            return;
        }

        // Manage Option i
        if (cmd.hasOption("i")) {
            Options optionsEmpty = new Options();
            formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND_WB, "",
                    optionsEmpty,
                    CommandLineMessages.COPYRIGHT
                            + CommandLineMessages.LONG_NOTICE
                            + CommandLineMessages.REPORT_TO
                            + "\n\nUse option -h for full usage information",
                    true);
            return;
        }

        // Manage option f
        if (cmd.hasOption("f")) {
            // Process a folder
            // Get folder file from command line
            File inputDirectory = new File(cmd.getOptionValue("f"));
            // Variables for counting execution time
            long startTimeNano = 0;
            long elapsedTimeNano = 0;
            // Maps for holding all the structures to be processed and their
            // associated
            // processing time
            Map<File, StructuralSequence> structures = new HashMap<File, StructuralSequence>();
            Map<File, Long> structuresProcessingTime = new HashMap<File, Long>();
            // List for holding all the structures files
            List<File> structuresList = new ArrayList<File>();
            // Set of skipped files
            Set<File> skippedFiles = new HashSet<File>();
            int numStructures = 1;

            // Process input files
            if (!inputDirectory.isDirectory()) {
                System.err.println("ERROR: Input file "
                        + cmd.getOptionValue("f") + " is not a folder");
                System.exit(1);
            }
            File[] fs = inputDirectory.listFiles();
            // Filter only files and put them in the list
            for (int i = 0; i < fs.length; i++)
                if (!fs[i].isDirectory())
                    if (!fs[i].isHidden())
                        structuresList.add(fs[i]);
                    else
                        System.err.println("WARNING: Skipping hidden file "
                                + fs[i].getName() + " ...");
                else
                    System.err.println("WARNING: Skipping subfolder "
                            + fs[i].getName() + " ...");

            // Order files to be processed
            Collections.sort(structuresList);

            // Output files creation
            PrintStream outputStream = null;
            PrintStream jsonOutputStream = null;
            PrintStream structuresStream = null;
            PrintStream jsonStructuresStream = null;
            String outputStreamName = inputDirectory.getAbsolutePath() + "/"
                    + "SERNAlignComparisonResults.csv";
            String structuresStreamName = inputDirectory.getAbsolutePath()
                    + "/" + "SERNAlignProcessedStructures.csv";
            String jsonOutputStreamName = inputDirectory.getAbsolutePath() + "/"
                    + "SERNAlignComparisonResults.json";
            String jsonStructuresStreamName = inputDirectory.getAbsolutePath()
                    + "/" + "SERNAlignProcessedStructures.json";


            List<String> jsonOutputEntries = new ArrayList<>();
            List<String> jsonStructuresEntries = new ArrayList<>();

            // Manage option "o"
            if (cmd.hasOption("o")) {
                String[] names = cmd.getOptionValues("o");
                structuresStreamName = names[0];
                outputStreamName = names[1];
                jsonStructuresStreamName = names[0] + ".json";
                jsonOutputStreamName = names[1] + ".json";
            }

            try {
                outputStream = new PrintStream(new File(outputStreamName));
                structuresStream = new PrintStream(
                        new File(structuresStreamName));
            } catch (FileNotFoundException e) {
                System.err.println("ERROR: creation of output file "
                        + (outputStream == null ? outputStreamName
                        : structuresStreamName)
                        + " failed");
                System.exit(3);
            }

            // JSON output

            if (cmd.hasOption("j")) {
                try {
                    jsonOutputStream = new PrintStream(new File(jsonOutputStreamName));
                    jsonStructuresStream =
                            new PrintStream(new File(jsonStructuresStreamName));
                } catch (FileNotFoundException e) {
                    System.err.println("ERROR: failed to create JSON output " +
                            "file... " + e.getMessage());
                    System.exit(4);
                }
            }

            // Write column names on the csv output files
            structuresStream.println(
                    "Num,FileName,NumberOfNucleotides,NumberOfWeakBonds,"
                            + "IsPseudoknotted,TimeToGenerateStructuralSequence[ns]");
            outputStream.println(
                    "FileName1,NumberOfNucleotides1,NumberOfWeakBonds1,IsPseudoknotted1,TimeToGenerateStructuralSequence1[ns],"
                            + "FileName2,NumberOfNucleotides2,NumberOfWeakBonds2,IsPseudoknotted2,TimeToGenerateStructuralSequence2[ns],"
                            + "MaxNumberOfNucleotides1-2,SERNADistance,TimeToCalculateSERNADistance[ns]");

            // Main Loop
            ListIterator<File> extIt = structuresList.listIterator();
            while (extIt.hasNext()) {
                // Compare the next element with all the subsequent elements
                int currentExtIndex = extIt.nextIndex();
                // Process File 1
                File f1 = extIt.next();
                // Check if skipped
                if (skippedFiles.contains(f1))
                    // skip this file
                    continue;

                // Retrieve the Structural RNA Tree for the structure 1
                StructuralSequence art1 = null;
                // Check if this structure has already been processed
                if (!structures.containsKey(f1)) {
                    // Parse the input file f1 for the secondary structure
                    RNASecondaryStructure secondaryStructure1 = null;
                    try {
                        secondaryStructure1 = RNASecondaryStructureFileReader
                                .readStructure(f1.getPath(), false);
                    } catch (IOException e) {
                        System.err.println("WARNING: Skipping file "
                                + f1.getName() + " ... " + e.getMessage());
                        // skip this structure
                        skippedFiles.add(f1);
                        continue;
                    } catch (RNAInputFileParserException e) {
                        System.err.println("WARNING: Skipping file "
                                + f1.getName() + " ... " + e.getMessage());
                        // skip this structure
                        skippedFiles.add(f1);
                        continue;
                    }
                    // Create the Structural RNA Tree and put the object into
                    // the map
                    // Build Structural Sequence and measure building time
                    startTimeNano = System.nanoTime();
                    art1 = new StructuralSequence(secondaryStructure1);
                    elapsedTimeNano = System.nanoTime() - startTimeNano;
                    // Insert Object in maps
                    structures.put(f1, art1);
                    structuresProcessingTime.put(f1, elapsedTimeNano);
                    // Output values in the structures output file
                    structuresStream.println(numStructures + "," + "\""
                            + f1.getName() + "\","
                            + art1.getSecondaryStructure().getSize() + ","
                            + art1.getSecondaryStructure().getBonds().size()
                            + ","
                            + (art1.getSecondaryStructure().isPseudoknotted()
                            ? "Yes"
                            : "No")
                            + "," + elapsedTimeNano);
                    if (jsonStructuresStream != null) {
                        String jsonStructEntry = String.format(
                                "{ \"Num\": %d, \"FileName\": \"%s\", \"NumberOfNucleotides\": %d, \"NumberOfWeakBonds\": %d, " +
                                        "\"IsPseudoknotted\": \"%s\", \"TimeToGenerateStructuralSequence_ns\": %d }",
                                numStructures,
                                f1.getName(),
                                art1.getSecondaryStructure().getSize(),
                                art1.getSecondaryStructure().getBonds().size(),
                                art1.getSecondaryStructure().isPseudoknotted() ? "Yes" : "No",
                                elapsedTimeNano
                        );
                        jsonStructuresEntries.add(jsonStructEntry);
                    }
                    numStructures++;
                } else {
                    art1 = structures.get(f1);
                }

                // Internal Loop - Compare structure 1 with all the subsequent
                // ones
                ListIterator<File> intIt = structuresList
                        .listIterator(currentExtIndex + 1);
                while (intIt.hasNext()) {
                    // Process File 2
                    File f2 = intIt.next();
                    // Check if skipped
                    if (skippedFiles.contains(f2))
                        // skip this file
                        continue;

                    // Retrieve the Structural Sequence for the structure 2
                    StructuralSequence art2 = null;
                    // Check if this structure has already been processed
                    if (!structures.containsKey(f2)) {
                        // Parse the input file f2 for the secondary structure
                        RNASecondaryStructure secondaryStructure2 = null;
                        try {
                            secondaryStructure2 = RNASecondaryStructureFileReader
                                    .readStructure(f2.getPath(), false);
                        } catch (IOException e) {
                            System.err.println(
                                    "WARNING: Skipping file " + f2.getName()
                                            + " ... " + e.getMessage());
                            // skip this structure
                            skippedFiles.add(f2);
                            continue;
                        } catch (RNAInputFileParserException e) {
                            System.err.println(
                                    "WARNING: Skipping file " + f2.getName()
                                            + " ... " + e.getMessage());
                            // skip this structure
                            skippedFiles.add(f2);
                            continue;
                        }

                        // Build Structural RNA Tree and measure building time
                        startTimeNano = System.nanoTime();
                        // Create the Structural RNA Tree and put the object
                        // into the map
                        art2 = new StructuralSequence(secondaryStructure2);
                        elapsedTimeNano = System.nanoTime() - startTimeNano;
                        // Insert Object in maps
                        structures.put(f2, art2);
                        structuresProcessingTime.put(f2, elapsedTimeNano);
                        // Output values in the structures output file
                        structuresStream.println(numStructures + "," + "\""
                                + f2.getName() + "\","
                                + art2.getSecondaryStructure().getSize() + ","
                                + art2.getSecondaryStructure().getBonds()
                                .size()
                                + ","
                                + (art2.getSecondaryStructure()
                                .isPseudoknotted() ? "Yes" : "No")
                                + "," + elapsedTimeNano);
                        if (jsonStructuresStream != null) {
                            String jsonStructEntry = String.format(
                                    "{ \"Num\": %d, \"FileName\": \"%s\", \"NumberOfNucleotides\": %d, \"NumberOfWeakBonds\": %d, " +
                                            "\"IsPseudoknotted\": \"%s\", \"TimeToGenerateStructuralSequence_ns\": %d }",
                                    numStructures,
                                    f1.getName(),
                                    art1.getSecondaryStructure().getSize(),
                                    art1.getSecondaryStructure().getBonds().size(),
                                    art1.getSecondaryStructure().isPseudoknotted() ? "Yes" : "No",
                                    elapsedTimeNano
                            );
                            jsonStructuresEntries.add(jsonStructEntry);
                        }
                        numStructures++;
                    } else {
                        art2 = structures.get(f2);
                    }

                    // Compare the two structural RNA Trees t1 and t2 to
                    // determine the distance
                    System.out.println("Processing files: " + f1.getName()
                            + " and " + f2.getName());
                    startTimeNano = System.nanoTime();
                    StructuralSequenceAligner a = new StructuralSequenceAligner(art1, art2, !noConstraints);
                    elapsedTimeNano = System.nanoTime() - startTimeNano;

                    // Write the output file
                    outputStream.println("\"" + f1.getName() + "\","
                            + art1.getSecondaryStructure().getSize() + ","
                            + art1.getSecondaryStructure().getBonds().size()
                            + ","
                            + (art1.getSecondaryStructure().isPseudoknotted()
                            ? "Yes"
                            : "No")
                            + ","
                            + structuresProcessingTime.get(f1).longValue()
                            + "," + "\"" + f2.getName() + "\","
                            + art2.getSecondaryStructure().getSize() + ","
                            + art2.getSecondaryStructure().getBonds().size()
                            + ","
                            + (art2.getSecondaryStructure().isPseudoknotted()
                            ? "Yes"
                            : "No")
                            + ","
                            + structuresProcessingTime.get(f2).longValue()
                            + ","
                            + (art1.getSecondaryStructure().getSize() > art2
                            .getSecondaryStructure().getSize()
                            ? art1.getSecondaryStructure()
                            .getSize()
                            : art2.getSecondaryStructure()
                            .getSize())
                            + "," + a.getDistance() + "," + elapsedTimeNano);

                    if (jsonOutputStream != null) {
                        String jsonEntry = String.format(
                                "{ \"File1\": \"%s\", \"Size1\": %d, \"Bonds1\": %d, \"Pseudoknotted1\": \"%s\", \"Time1_ns\": %d, " +
                                        "\"File2\": \"%s\", \"Size2\": %d, \"Bonds2\": %d, \"Pseudoknotted2\": \"%s\", \"Time2_ns\": %d, " +
                                        "\"MaxSize\": %d, \"Distance\": %d, " +
                                        "\"AlignmentTime_ns\": %d }",
                                f1.getName(), art1.getSecondaryStructure().getSize(),
                                art1.getSecondaryStructure().getBonds().size(),
                                art1.getSecondaryStructure().isPseudoknotted() ? "Yes" : "No",
                                structuresProcessingTime.get(f1),
                                f2.getName(), art2.getSecondaryStructure().getSize(),
                                art2.getSecondaryStructure().getBonds().size(),
                                art2.getSecondaryStructure().isPseudoknotted() ? "Yes" : "No",
                                structuresProcessingTime.get(f2),
                                Math.max(art1.getSecondaryStructure().getSize(), art2.getSecondaryStructure().getSize()),
                                a.getDistance(), elapsedTimeNano
                        );
                        jsonOutputEntries.add(jsonEntry);
                    }
                    // End of Internal Loop
                }
                // End of External Loop
            }

            // Close streams
            structuresStream.close();
            outputStream.close();
            if (jsonStructuresStream != null) {
                jsonStructuresStream.println("[");
                for (int i = 0; i < jsonStructuresEntries.size(); i++) {
                    if (i < jsonStructuresEntries.size() - 1) {
                        jsonStructuresStream.println(jsonStructuresEntries.get(i) + ",");
                    } else {
                        jsonStructuresStream.println(jsonStructuresEntries.get(i));
                    }
                }
                jsonStructuresStream.println("]");
                jsonStructuresStream.close();
            }
            if (jsonOutputStream != null) {
                jsonOutputStream.println("[");
                for (int i = 0; i < jsonOutputEntries.size(); i++) {
                    if (i < jsonOutputEntries.size() - 1) {
                        jsonOutputStream.println(jsonOutputEntries.get(i) + ",");
                    } else {
                        jsonOutputStream.println(jsonOutputEntries.get(i));
                    }
                }
                jsonOutputStream.println("]");
                jsonOutputStream.close();
            }
            return;
        } // End Option f

        // If no option is given, output usage
        formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND_WB,
                CommandLineMessages.HEADER_WB, options,
                CommandLineMessages.USAGE_EXAMPLES_WB
                        + CommandLineMessages.COPYRIGHT
                        + CommandLineMessages.SHORT_NOTICE
                        + CommandLineMessages.REPORT_TO,
                true);
    }

}
