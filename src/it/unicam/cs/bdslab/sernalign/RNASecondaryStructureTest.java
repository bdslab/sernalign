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

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;

import org.junit.jupiter.api.Test;

/**
 * Tests for the class RNASecondaryStructure
 * 
 * @author Luca Tesei
 *
 */
class RNASecondaryStructureTest {

    @Test
    void testIsPseudoknotted() throws IOException {
	RNASecondaryStructure s1 = RNASecondaryStructureFileReader
		.readStructure("test/CRW_16S_A_C_1.db", false);
	RNASecondaryStructure s2 = RNASecondaryStructureFileReader
		.readStructure("test/CRW_16S_A_C_19.db", false);
	RNASecondaryStructure s3 = RNASecondaryStructureFileReader
		.readStructure("test/CRW_5S_A_C_20.db", false);
	RNASecondaryStructure s4 = RNASecondaryStructureFileReader
		.readStructure("test/CRW_5S_A_C_22.db", false);
	assertTrue(s1.isPseudoknotted());
	assertTrue(s2.isPseudoknotted());
	assertFalse(s3.isPseudoknotted());
	assertFalse(s4.isPseudoknotted());

    }

}
