/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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

import java.util.ArrayList;
import java.util.List;

/**
 * Simple representation of a chromosome.  Basically a name, length, and optionally a list of cytobands.
 *
 * @author jrobinso
 */
public class ChromosomeImpl implements Chromosome {
    private String name;
    private int length = 0;
    private List<Cytoband> cytobands;

    public ChromosomeImpl(String name) {
        this.name = name;
    }

    public ChromosomeImpl(String name, int length) {
        this.name = name;
        this.length = length;
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

}
