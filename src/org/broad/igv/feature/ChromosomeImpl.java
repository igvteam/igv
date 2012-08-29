/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
public class ChromosomeImpl implements Chromosome {
    private String name;
    private int index;  // Order in the chromosome (for convenience)
    private int length = 0;
    private List<Cytoband> cytobands;

    public ChromosomeImpl(int index, String name, int length) {
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


    @Override
    public int getIndex() {
        return index;
    }

    @Override
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
