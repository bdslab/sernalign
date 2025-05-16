/**
 * SERNAlign - Structural sEquence RNA secondary structure Alignment
 * 
 * Copyright (C) 2023 Luca Tesei, Francesca Levi, Michela Quadrini, 
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

import java.util.LinkedList;
import java.util.List;

/**
 * Structural Sequence Alignment STREQAlign
 * 
 * @author Luca Tesei
 *
 */
public class StructuralSequenceAligner {

    private final StructuralSequence x;
    private final StructuralSequence y;

    /* Constants for the alignment matrix */
    private static int STOP = -1;
    private static int DIAGONAL_MATCH_MISMATCH = 0;
    private static int UP_DELETION = 1;
    private static int LEFT_INSERTION = 2;

    /* matrix for alignment */
    private int[][] m;
    /* matrix for traceback */
    private int[][] traceback;
    /* flag to indicate if the alignment should respect structural sequences constraints */
	private final boolean constraints;

    /* best alignment as a list of edit operations */
    LinkedList<EditOperation> alignment;

    /**
     * Construct a minimum alignment to transform a structural sequence into
     * another respecting structural sequences constraints.
     * 
     * @param x the first structural sequence to align
     * @param y the second structural sequence to align
	 * @param constraints flag to indicate if the alignment should respect structural
	 *                    sequences  constraints
	 * @throws NullPointerException if one of the two sequences is null
	 *
     */
    public StructuralSequenceAligner(StructuralSequence x,
	    StructuralSequence y,
	 	boolean constraints
	) {
	if (x == null || y == null)
	    throw new NullPointerException(
		    "Tentativo di definire un aligner tra due sequenze di cui almeno una è null");
	this.x = x;
	this.y = y;
	this.m = new int[x.size() + 1][y.size() + 1];
	this.traceback = new int[x.size() + 1][y.size() + 1];
	this.constraints = constraints;
	solve();
	this.alignment = new LinkedList<EditOperation>();
	traceBack(this.x.size(), this.y.size());
    }

	/**
	 * Construct a minimum alignment to transform a structural sequence into
	 * another respecting structural sequences constraints.
	 *
	 * @param x the first structural sequence to align
	 * @param y the second structural sequence to align
	 */
	public StructuralSequenceAligner(StructuralSequence x,
									 StructuralSequence y) {
		this(x, y, true);
	}

    private void solve() {
	// initialize first raw and first column
	for (int i = 0; i < this.m.length; i++) {
	    this.m[i][0] = i;
	    this.traceback[i][0] = STOP;
	}
	for (int j = 0; j < this.m[0].length; j++) {
	    this.m[0][j] = j;
	    this.traceback[0][j] = STOP;
	}
	for (int i = 1; i < this.m.length; i++)
	    for (int j = 1; j < this.m[0].length; j++) {
		int p;
		if (this.x.getStructuralSequence()[i - 1] == this.y
			.getStructuralSequence()[j - 1])
		    p = 0; // possible match
		else
		    p = 1; // match not possible
		int valInsertion;
		int valDeletion;
		int valMatchMismatch;
		// Insertion
		if (isCorrectInPosition(this.y.getStructuralSequence()[j - 1],
			i))
		    // insertion is possible
		    valInsertion = this.m[i][j - 1] + 1;
		else
		    // insertion is not possible in this case
		    valInsertion = Integer.MAX_VALUE;
		// Deletion
		if (isCorrectInPosition(this.x.getStructuralSequence()[i - 1],
			j))
		    // deletion is possible
		    valDeletion = this.m[i - 1][j] + 1;
		else
		    // deletion is not possible in this case
		    valDeletion = Integer.MAX_VALUE;
		// Match / Mismatch
		if (isCorrectInPosition(this.y.getStructuralSequence()[j - 1],
			i)
			&& isCorrectInPosition(
				this.x.getStructuralSequence()[i - 1], j))
		    // Match / Mismatch is possible
		    valMatchMismatch = this.m[i - 1][j - 1] + p;
		else
		    // Match / Mismatch is not possible in this case
		    valMatchMismatch = Integer.MAX_VALUE;
		/*
		 * Determine the minimum of the three values
		 */
		int min = valInsertion;
		int min_direction = LEFT_INSERTION;
		if (valDeletion < min) {
		    min = valDeletion;
		    min_direction = UP_DELETION;
		}
		if (valMatchMismatch < min) {
		    min = valMatchMismatch;
		    min_direction = DIAGONAL_MATCH_MISMATCH;
		}
		// write the value in the matrix
		this.m[i][j] = min;
		this.traceback[i][j] = min_direction;
		//System.out.println(printMatrix());
	    }
    }

    private boolean isCorrectInPosition(int h, int pos) {
		// Tells whether the value h can be in position pos in a numerical
		// sequence
		return !this.constraints || h >= 1 && h <= 2 * pos - 1;
    }

	private static boolean isCorrectInPositionConstraint(int h, int pos) {
		return h >= 1 && h <= 2 * pos - 1;
	}

    /*
     * For testing purposes
     */
    protected int[][] getMatrix() {
	return this.m;
    }

    public int getDistance() {
	return this.m[this.x.size()][this.y.size()];
    }

    private void traceBack(int i, int j) {
	if (i == 0 && j == 0) {
	    // finished
	    return;
	}
	if (i == 0 && j > 0) {
	    /*
	     * insert elements into the empty sequence towards y
	     */
	    for (int k = j; k > 0; k--)
		this.alignment.addFirst(new EditOperation(null,
			this.y.getStructuralSequence()[k - 1]));
	    /* Return */
	    return;
	} else if (i > 0 && j == 0) {
	    /*
	     * delete elements of x towards the empty sequence
	     */
	    for (int k = i; k > 0; k--)
		this.alignment.addFirst(new EditOperation(
			this.x.getStructuralSequence()[k - 1], null));
	    /* Return */
	    return;
	} else { /*
		  * match/mismatch, deletion or insertion when we are
		  * considering two non empty prefixes
		  */
	    if (this.traceback[i][j] == DIAGONAL_MATCH_MISMATCH) {
		// add match / mismatch edit operation
		this.alignment.addFirst(new EditOperation(
			this.x.getStructuralSequence()[i - 1],
			this.y.getStructuralSequence()[j - 1]));
		// recursive call
		traceBack(i - 1, j - 1);
	    } else if (this.traceback[i][j] == UP_DELETION) {
		// add deletion edit operation
		this.alignment.addFirst(new EditOperation(
			this.x.getStructuralSequence()[i - 1], null));
		// recursive call
		traceBack(i - 1, j);
	    } else {
		// add insertion edit operation
		this.alignment.addFirst(new EditOperation(null,
			this.y.getStructuralSequence()[j - 1]));
		// recursive call
		traceBack(i, j - 1);
	    }
	}
    }

    private String printSeq(List<Integer> seq) {
	StringBuffer s = new StringBuffer();
	for (int i = 0; i < seq.size() - 1; i++)
	    s.append(seq.get(i) + ", ");
	s.append(seq.get(seq.size() - 1) + "");
	return s.toString();
    }

    public String printOptimalAlignmentExecution() {
	List<Integer> seq = this.x.getStructuralSequenceAsList();
	int i = this.x.size();
	int j = this.y.size();
	int ell = this.alignment.size() - 1;
	StringBuffer s = new StringBuffer();
	s.append(printSeq(seq) + "\n");
	while (i > 0 || j > 0) {
	    if (this.alignment.get(ell).getJ() == null) {
		// remove operation
		s.append("( " + this.x.getStructuralSequence()[i - 1]
			+ ", - )");
		seq.remove(seq.size() - 1);
		i--;
		ell--;
	    } else if (this.alignment.get(ell).getI() == null) {
		// insertion operation
		s.append("( " + this.x.getStructuralSequence()[i - 1] + ", "
			+ this.y.getStructuralSequence()[j - 1] + " )");
		seq.add(this.y.getStructuralSequence()[j - 1]);
		j--;
		ell--;
	    } else {
		// match/mismatch operation
		s.append("( " + this.x.getStructuralSequence()[i - 1] + ", "
			+ this.y.getStructuralSequence()[j - 1] + " )");
		seq.set(seq.size() - 1,
			this.y.getStructuralSequence()[j - 1]);
		i--;
		j--;
		ell--;
	    }
	    s.append("\n" + printSeq(seq) + "\n");
	}
	return s.toString();
    }

    public String printOptimalAlignmentConstraints() {
	return printAlignmentConstraints(this.x.size(), this.y.size(),
		this.alignment.size() - 1);
    }

    private String printAlignmentConstraints(int i, int j, int ell) {
	if (ell == -1 || i == 0 || j == 0)
	    return "";
	if (this.alignment.get(ell).getJ() == null) {
	    // deletion operation
	    return "x_" + i + " = " + this.x.getStructuralSequence()[i - 1]
		    + " <=  " + (j * 2 - 1) + " = C_" + j + "\n"
		    + printAlignmentConstraints(i - 1, j, ell - 1);
	}

	if (this.alignment.get(ell).getI() == null) {
	    // insertion operation
	    return "y_" + j + " = " + this.y.getStructuralSequence()[j - 1]
		    + " <=  " + (i * 2 - 1) + " = C_" + i + "\n"
		    + printAlignmentConstraints(i, j - 1, ell - 1);
	}

	// match/mismatch operation
	return "x_" + i + " = " + this.x.getStructuralSequence()[i - 1]
		+ " <=  " + (j * 2 - 1) + " = C_" + j + " and " + "y_" + j
		+ " = " + this.y.getStructuralSequence()[j - 1] + " <=  "
		+ (i * 2 - 1) + " = C_" + i + "\n"
		+ printAlignmentConstraints(i - 1, j - 1, ell - 1);
    }

    public boolean checkOptimalAlignment() {
	return checkAlignment(this.x.size(), this.y.size(),
		this.alignment.size() - 1);
    }

    private boolean checkAlignment(int i, int j, int ell) {
	if (ell == -1 || i == 0 || j == 0)
	    return true;
	if (this.alignment.get(ell).getJ() == null) {
	    // deletion operation
	    return isCorrectInPositionConstraint(this.x.getStructuralSequence()[i - 1],
		    j) && checkAlignment(i - 1, j, ell - 1);
	}

	if (this.alignment.get(ell).getI() == null) {
	    // insertion operation
	    return isCorrectInPositionConstraint(this.y.getStructuralSequence()[j - 1],
		    i) && checkAlignment(i, j - 1, ell - 1);
	}

	// match/mismatch operation
	return isCorrectInPositionConstraint(this.x.getStructuralSequence()[i - 1], j)
		&& isCorrectInPositionConstraint(this.y.getStructuralSequence()[j - 1],
			i)
		&& checkAlignment(i - 1, j - 1, ell - 1);
    }

    public String printOptimalAlignment() {
	StringBuffer s = new StringBuffer();
	String sI = null;
	String sJ = null;
	for (EditOperation o : this.alignment) {
	    if (o.getI() == null)
		sI = "-";
	    else
		sI = "" + o.getI();
	    if (o.getJ() == null)
		sJ = "-";
	    else
		sJ = "" + o.getJ();
	    s.append("(" + sI + ", " + sJ + ")");
	}
	return s.toString();
    }

	public String printMatrix() {
		StringBuilder matrix = new StringBuilder();
		for (int i = 0; i < this.m.length; i++) {
			for (int j = 0; j < this.m[0].length; j++) {
				matrix.append(this.m[i][j]);
				if (j < this.m[0].length - 1) {
					matrix.append(", ");
				}
			}
			matrix.append("\n");
		}
		return matrix.toString();
	}


    public List<EditOperation> getOptimalAlignment() {
	return this.alignment;
    }

}
