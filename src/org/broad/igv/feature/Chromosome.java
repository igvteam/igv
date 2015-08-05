/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


/*
Chromosome.java
 *
Created on June 29, 2007, 8:53 AM
 *
To change this template, choose Tools | Template Manager
and open the template in the editor.
 */
package org.broad.igv.feature;

import java.util.Arrays;
import java.util.List;

/**
 * Simple representation of a chromosome.  Basically a name, length, and optionally a list of cytobands.
 *
 * @author jrobinso
 */
public class Chromosome {
    private String name;
    private int index;  // Order in the chromosome (for convenience)
    private int length = 0;
    private List<Cytoband> cytobands;

    public Chromosome(int index, String name, int length) {
        this.index = index;
        this.name = name;
        this.length = length;

        // Create a single "cytoband" to represent the entire chromosome.  This can be overriden explicitly
        // if a cytoband file is loaded
        final Cytoband cytoband = new Cytoband(name);
        cytoband.setStart(0);
        cytoband.setEnd(length);
        cytobands = Arrays.asList(cytoband);

    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int ii) {
        this.index = ii;
    }

    /**
     * @return List of cytobands for this chromosome, if any.  Can be null.
     */
    public List<Cytoband> getCytobands() {
        return cytobands;
    }


    public void setCytobands(List<Cytoband> cytobands) {
        this.cytobands = cytobands;
    }


    /**
     * /**
     * Return the length of the chromosome, which is the end of the last cytoband
     */
    public int getLength() {
        return length;
    }

    public String getName() {
        return name;
    }

    public String toString() {
        return name;
    }

    public boolean equals(Object obj) {
        return ((Chromosome)obj).getIndex() == index;
    }
}
