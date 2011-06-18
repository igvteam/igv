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
 * @author jrobinso
 */
public class Chromosome {
    private String name;
    private int centrimere = 0;
    private List<Cytoband> cytobands;
    private int length = 0;

    public Chromosome(String name) {

        this.name = name;
        cytobands = new ArrayList<Cytoband>();
    }

    public int getCentrimere() {
        return centrimere;
    }

    public void setCentrimere(int centrimere) {
        this.centrimere = centrimere;
    }

    public List<Cytoband> getCytobands() {
        return cytobands;
    }

    public void setCytobands(List<Cytoband> cytobands) {
        this.cytobands = cytobands;
    }

    /**
     * Add a cytoband.  If the band is a centrimere and the centrimere location
     * has not been initialized set it.
     */
    public void addCytoband(Cytoband band) {

        cytobands.add(band);
        if ((band.getType() == 'c') && (centrimere == 0)) {
            centrimere = band.getEnd();
        }
        length = Math.max(length, band.getEnd());
    }

    /**
     * Return the length of the chromosome, which is the end of the last cytoband
     */
    public int getLength() {
        return length;
    }

    public String getName() {
        return name;
    }

}
