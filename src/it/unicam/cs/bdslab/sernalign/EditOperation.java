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
 * A simple class to represent an edit operation. It is a pair (i,j), if i and
 * j are both non-null it is a match/mismatch operation, if i is null,
 * representing the gap, it is an insertion operation. If j is null, it is a
 * delete operation. i and j cannot be both null.
 */
public class EditOperation {

    private final Integer i;
    private final Integer j;

    /**
     * Create an edit operation.
     * 
     * @param i first element of the edit operation, if null corresponds to
     *          the gap
     * @param j second element of the edit operation, if null corresponds to
     *          the gap
     * 
     * @throws IllegalArgumentException if both values are null
     */
    public EditOperation(Integer i, Integer j) {
	if ((i == null) && (j == null))
	    throw new IllegalArgumentException(
		    "Attempt of creating an edit operation with two gaps");
	this.i = i;
	this.j = j;
    }

    /**
     * @return the first component of the edit operation
     */
    public Integer getI() {
	return i;
    }

    /**
     * @return the second component of the edit operation 
     */
    public Integer getJ() {
	return j;
    }

}
