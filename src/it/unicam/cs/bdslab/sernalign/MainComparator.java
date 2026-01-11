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

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * MainComparator class interacting with the user through command line
 * options.
 * 
 * @author Luca Tesei
 *
 */
public class MainComparator {

    public static void main(String[] args) {
	// Use Apache Commons CLI 1.4
	// create Options object for Command Line Definition
	Options options = new Options();

	defineCommandLineOptions(options);

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
	    formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND,
		    CommandLineMessages.HEADER, options,
		    CommandLineMessages.USAGE_EXAMPLES
			    + CommandLineMessages.COPYRIGHT
			    + CommandLineMessages.SHORT_NOTICE
			    + CommandLineMessages.REPORT_TO,
		    true);
	    System.exit(1);
	}

	// Manage Option h
	if (cmd.hasOption("h")) {
	    formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND,
		    CommandLineMessages.HEADER, options,
		    CommandLineMessages.USAGE_EXAMPLES
			    + CommandLineMessages.COPYRIGHT
			    + CommandLineMessages.SHORT_NOTICE
			    + CommandLineMessages.REPORT_TO,
		    true);
	    return;
	}

	// Manage Option i
	if (cmd.hasOption("i")) {
	    Options optionsEmpty = new Options();
	    formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND, "",
		    optionsEmpty,
		    CommandLineMessages.COPYRIGHT
			    + CommandLineMessages.LONG_NOTICE
			    + CommandLineMessages.REPORT_TO
			    + "\n\nUse option -h for full usage information",
		    true);
	    return;
	}

	// Manage option s
	if (cmd.hasOption("s")) {
	    // Determine input file
	    String inputFile = cmd.getOptionValue("s");
	    // Parse the input file for the secondary structure
	    RNASecondaryStructure secondaryStructure = null;
	    try {
		secondaryStructure = RNASecondaryStructureFileReader
			.readStructure(inputFile, false);
	    } catch (IOException e) {
		System.err.println("Input File " + inputFile + " ERROR:"
			+ e.getMessage());
		System.exit(3);
	    } catch (RNAInputFileParserException e) {
		System.err.println("Input File " + inputFile + " ERROR: "
			+ e.getMessage());
		System.exit(2);
	    }

	    // Construct the SRNSSequence
	    StructuralSequence a = new StructuralSequence(secondaryStructure);
	    // SRNSATree a = new SRNSATree(secondaryStructure);

	    String output = a.printStructuralSequence();

	    // Write Output on proper file or on standard output
	    if (cmd.hasOption("o")) {
		String outputFile = cmd.getOptionValue("o");
		try (BufferedWriter writer = new BufferedWriter(
			new FileWriter(outputFile, false))) {
		    writer.write(output + "\n");
		    writer.close();
		} catch (FileNotFoundException e) {
		    System.err.println("ERROR: Output file " + outputFile
			    + " cannot be created.");
		    System.exit(3);
		} catch (IOException e) {
		    System.err.println("ERROR: Output file" + outputFile
			    + " cannot be written.");
		    System.exit(3);
		}
	    } else
		System.out.println(output);
	    return;
	}

	// Manage Option a
	if (cmd.hasOption("a")) {
	    // Get input file names
	    String[] inputFiles = cmd.getOptionValues("a");

	    // Parse the first input file for the secondary structure
	    RNASecondaryStructure secondaryStructure1 = null;
	    try {
		secondaryStructure1 = RNASecondaryStructureFileReader
			.readStructure(inputFiles[0], false);
	    } catch (IOException e) {
		System.err.println("Input File " + inputFiles[0] + " ERROR:"
			+ e.getMessage());
		System.exit(3);
	    } catch (RNAInputFileParserException e) {
		System.err.println("Input File " + inputFiles[0] + " ERROR: "
			+ e.getMessage());
		System.exit(2);
	    }
	    // Construct structural RNA tree 1

	    StructuralSequence s1 = new StructuralSequence(
		    secondaryStructure1);

	    // Parse the second input file for the secondary structure
	    RNASecondaryStructure secondaryStructure2 = null;
	    try {
		secondaryStructure2 = RNASecondaryStructureFileReader
			.readStructure(inputFiles[1], false);
	    } catch (IOException e) {
		System.err.println("Input File " + inputFiles[1] + " ERROR:"
			+ e.getMessage());
		System.exit(3);
	    } catch (RNAInputFileParserException e) {
		System.err.println("Input File " + inputFiles[1] + " ERROR: "
			+ e.getMessage());
		System.exit(2);
	    }
	    // Construct structural sequence 2
	    StructuralSequence s2 = new StructuralSequence(
		    secondaryStructure2);

	    // Align t1 and t2, which contain two structural RNA trees
	    StructuralSequenceAligner a = new StructuralSequenceAligner(s1,
		    s2,
			!cmd.hasOption("n")
		);

	    // Produce Output
	    String output = a.printOptimalAlignment();
	    ;
	    String distance = "Distance = " + a.getDistance();

	    if (!cmd.hasOption("d")) {
		// Write Output on proper file or on standard output
		if (cmd.hasOption("o")) {
		    String outputFile = cmd.getOptionValue("o");
		    try (BufferedWriter writer = new BufferedWriter(
			    new FileWriter(outputFile, false))) {
			writer.write(output + "\n\n" + distance + "\n");
			writer.close();
		    } catch (FileNotFoundException e) {
			System.err.println("ERROR: Output file " + outputFile
				+ " cannot be created.");
			System.exit(3);
		    } catch (IOException e) {
			System.err.println("ERROR: Output file" + outputFile
				+ " cannot be written.");
			System.exit(3);
		    }
		} else
		    System.out.println(output + "\n\n" + distance);
	    } else { // Output only distance
		// Write Output on proper file or on standard output
		if (cmd.hasOption("o")) {
		    String outputFile = cmd.getOptionValue("o");
		    try (BufferedWriter writer = new BufferedWriter(
			    new FileWriter(outputFile, false))) {
			writer.write(distance + "\n");
			writer.close();
		    } catch (FileNotFoundException e) {
			System.err.println("ERROR: Output file " + outputFile
				+ " cannot be created.");
			System.exit(3);
		    } catch (IOException e) {
			System.err.println("ERROR: Output file" + outputFile
				+ " cannot be written.");
			System.exit(3);
		    }
		} else
		    System.out.println(distance);
	    }
	    return;
	}

	// If no option is given, output usage
	System.err.println("No operation specified...");
	formatter.printHelp(CommandLineMessages.LAUNCH_COMMAND,
		CommandLineMessages.HEADER, options,
		CommandLineMessages.USAGE_EXAMPLES
			+ CommandLineMessages.COPYRIGHT
			+ CommandLineMessages.SHORT_NOTICE
			+ CommandLineMessages.REPORT_TO,
		true);
    }

    private static void defineCommandLineOptions(Options options) {
	// define command line options
	Option o2 = Option.builder("s").desc(
		"Produce the structural sequence corresponding to the given "
			+ "structure")
		.longOpt("struct").hasArg().argName("input-file").build();
	options.addOption(o2);
	Option o3 = Option.builder("a").desc(
		"Align two given structures producing an alignment and distance")
		.longOpt("align").hasArgs().numberOfArgs(2)
		.argName("input-file1 input-file2").build();
	options.addOption(o3);
	Option o4 = Option.builder("o").desc(
		"Output result on the given file instead of standard output")
		.longOpt("out").hasArg().argName("output-file").build();
	options.addOption(o4);
	Option o6 = Option.builder("i").desc("Show license and other info")
		.longOpt("info").build();
	options.addOption(o6);
	Option o7 = Option.builder("h").desc("Show usage information")
		.longOpt("help").build();
	options.addOption(o7);
	Option o10 = Option.builder("d")
		.desc("Output only distance, no alignment (works only with "
			+ "option -a)")
		.longOpt("outdist").build();
	options.addOption(o10);

	options.addOption(
			"n",
			"no-constraints",
			false,
			"Do not use constraints for the alignment (works only with option -a)");
    }
}
