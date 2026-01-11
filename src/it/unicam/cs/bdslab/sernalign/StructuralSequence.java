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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This class produces a Structural Sequence from a given RNA Secondary
 * Structure with arbitrary pseudoknots. The Structural Sequence is a
 * numerical sequence representing the topological structure of the arc
 * diagram associated to the RNA secondary structure.
 * 
 * Each RNA secondary structure has associated one and only one structural
 * sequence.
 * 
 * @author Luca Tesei
 *
 */
public class StructuralSequence {
    // original RNA secondary structure
    private final RNASecondaryStructure secondaryStructure;
    // structural sequence
    private int[] structuralSequence;

    private Context head;
    private Context tail;

    /**
     * Construct the structural sequence associated to a given RNA secondary
     * structure.
     * 
     * @param secondaryStructure an RNA secondary structure
     */
    public StructuralSequence(RNASecondaryStructure secondaryStructure) {
	this.secondaryStructure = secondaryStructure;
	this.head = null;
	this.tail = null;
	buildSequence();
    }

    /**
     * @return the original RNA secondary structure
     */
    public RNASecondaryStructure getSecondaryStructure() {
	return this.secondaryStructure;
    }

    /**
     * @return the Structural Sequence
     */
    public int[] getStructuralSequence() {
	return structuralSequence;
    }

    /**
     * @return the Structural Sequence
     */
    public List<Integer> getStructuralSequenceAsList() {
	List<Integer> result = new ArrayList<Integer>();
	for (int i = 0; i < this.structuralSequence.length; i++)
	    result.add(this.structuralSequence[i]);
	return result;
    }

    /**
     * Construct a string representation of the structural sequence.
     * 
     * @return a string representation of this structural sequence
     */
    public String printStructuralSequence() {
	StringBuffer s = new StringBuffer();
	s.append("[ ");
	for (int i = 0; i < this.structuralSequence.length - 1; i++)
	    s.append(this.structuralSequence[i] + ", ");
	s.append(this.structuralSequence[this.structuralSequence.length - 1]
		+ " ]");
	return s.toString();
    }

    /**
     * 
     * @return the length of the structural sequence, also equal to the number
     *         of weak bonds in the original structure
     */
    public int size() {
	return this.structuralSequence.length;
    }

    /*
     * Build the structural sequence associated to the RNA secondary
     * structure. To do so scan the bonds in their natural order and determine
     * for each bond the context on which its left nucleotide is inserted. The
     * context is determined on the sub-structure composed of all the bonds
     * already scanned.
     * 
     * If n is the number of bonds in the structure, this method produces the
     * structural sequence in O(n^2). In a secondary structure of length m the
     * number of weak bonds is O(m).
     */
    private void buildSequence() {
	// Recover the bonds from the structure
	List<WeakBond> bonds = this.secondaryStructure.getBonds();
	// Order the bonds, if they are not yet ordered
	Collections.sort(bonds);
	// Create the array for the structural sequence
	this.structuralSequence = new int[bonds.size()];
	int i = 0;
	for (WeakBond b : bonds) {
	    if (i == 0) {
		// just insert the first bond in the double linked list and
		// assign the first position of the sequence
		this.head = new Context(0, b.getLeft(), null, null);
		this.tail = new Context(b.getRight(), Integer.MAX_VALUE, null,
			null);
		this.head.next = new Context(b.getLeft(), b.getRight(),
			this.head, this.tail);
		this.tail.previous = this.head.next;
		this.structuralSequence[i] = 1;
	    } else {
		// find the context and (index) in which the current left
		// nucleotide is
		Context scan = this.head;
		int index = 1;
		while (true) {
		    if (contained(b.getLeft(), scan))
			break;
		    index++;
		    scan = scan.next;
		}
		// assign the current position of the structural sequence
		this.structuralSequence[i] = index;
		// split the context in two
		split(scan, b.getLeft());
		// split the tail context in two
		split(this.tail, b.getRight());
	    }
	    i++;
	}

    }

    /*
     * Split a given context by a given position creating a new element of the
     * double linked list and updating all the references and the positions of
     * the contexts. Manage the splitting of the head and the tail correctly,
     * so it can be safely called on any context.
     */
    private void split(Context scan, int pos) {
	if (scan == this.head) {
	    // splitting the head
	    Context next = scan.next;
	    scan.next = new Context(pos, scan.j, scan, next);
	    scan.j = pos;
	    next.previous = scan.next;
	} else {
	    // splitting the tail or in the middle
	    Context previous = scan.previous;
	    previous.next = new Context(scan.i, pos, previous, scan);
	    scan.i = pos;
	    scan.previous = previous.next;
	}
    }

    /*
     * Return true if the position belongs to the context.
     */
    private boolean contained(int pos, Context scan) {
	return pos >= scan.i && pos <= scan.j;
    }

    /*
     * This class is used to realize a double linked list of intervals
     * representing contexts. Each interval is composed by the starting
     * position and the ending position of a context. The first context always
     * starts at 0 and the last context always ends at Integer.MAX_VALUE.
     */
    private static class Context {
	private int i;
	private int j;
	private Context previous;
	private Context next;

	public Context(int i, int j, Context previous, Context next) {
	    this.i = i;
	    this.j = j;
	    this.previous = previous;
	    this.next = next;
	}

	@SuppressWarnings("unused")
	public Context(int i, int j) {
	    this.i = i;
	    this.j = j;
	}

    }

}
