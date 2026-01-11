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

import static org.junit.jupiter.api.Assertions.*;

import java.io.IOException;
import java.util.Arrays;

import org.junit.jupiter.api.Test;

class StructuralSequenceTest {

    @Test
    void testGetNumericalSequence1() throws IOException {
	RNASecondaryStructure s1 = RNASecondaryStructureFileReader
		.readStructure("test/simple-example-1-aas.txt", false);
	int[] seq = { 1, 3, 1, 1, 5 };
	StructuralSequence seq1 = new StructuralSequence(s1);
	assertTrue(Arrays.equals(seq, seq1.getStructuralSequence()));
    }

    @Test
    void testGetNumericalSequence2() throws IOException {
	RNASecondaryStructure s1 = RNASecondaryStructureFileReader
		.readStructure("test/simple-example-2-aas.txt", false);
	int[] seq = { 1, 1, 1, 4 };
	StructuralSequence seq1 = new StructuralSequence(s1);
	assertTrue(Arrays.equals(seq, seq1.getStructuralSequence()));
    }

    @Test
    void testGetNumericalSequence3() throws IOException {
	RNASecondaryStructure s1 = RNASecondaryStructureFileReader
		.readStructure("test/larger-example-1-aas.txt", false);
	int[] seq = { 1, 1, 1, 5, 9, 9, 13, 14, 13, 9, 20, 21, 25, 26 };
	StructuralSequence seq1 = new StructuralSequence(s1);
	for (int i = 0; i < seq1.getStructuralSequence().length; i++)
	    System.out.println(seq1.getStructuralSequence()[i] + " ");
	assertTrue(Arrays.equals(seq, seq1.getStructuralSequence()));
    }

    @Test
    void testGetNumericalSequence4() throws IOException {
	RNASecondaryStructure s1 = RNASecondaryStructureFileReader
		.readStructure("test/larger-example-2-aas.txt", false);
	int[] seq = { 1, 2, 1, 5, 9, 9, 12, 15, 16, 15, 15, 9, 24, 26};
	StructuralSequence seq1 = new StructuralSequence(s1);
	assertTrue(Arrays.equals(seq, seq1.getStructuralSequence()));
    }

    @Test
    void testGetNumericalSequence5() throws IOException {
	RNASecondaryStructure s1 = RNASecondaryStructureFileReader
		.readStructure("test/larger-example-3-aas.txt", false);
	int[] seq = { 1, 1, 1, 5, 9, 9, 13, 14, 13, 9, 20, 21, 25, 26 };
	StructuralSequence seq1 = new StructuralSequence(s1);
	assertTrue(Arrays.equals(seq, seq1.getStructuralSequence()));
    }
}
