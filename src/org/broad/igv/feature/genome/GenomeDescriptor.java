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
 * GenomeType.java
 *
 * Created on November 8, 2007, 4:20 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

/**
 * @author eflakes
 */
public abstract class GenomeDescriptor {

    private String name;
    //private int version;
    private boolean chrNamesAltered;
    private String id;
    protected String cytoBandFileName;
    protected String geneFileName;
    protected String chrAliasFileName;
    private String geneTrackName;
    private String url;
    private String sequenceLocation;
    private boolean chromosomesAreOrdered = false;
    private boolean fasta = false;
    private boolean fastaDirectory = false;
    private String [] fastaFileNames;

    public GenomeDescriptor(String name,
                            boolean chrNamesAltered,
                            String id,
                            String cytoBandFileName,
                            String geneFileName,
                            String chrAliasFileName,
                            String geneTrackName,
                            String sequenceLocation,
                            boolean chromosomesAreOrdered,
                            boolean fasta,
                            boolean fastaDirectory,
                            String fastaFileNameString) {
        this.chrNamesAltered = chrNamesAltered;
        this.name = name;
        this.id = id;
        this.cytoBandFileName = cytoBandFileName;
        this.geneFileName = geneFileName;
        this.chrAliasFileName = chrAliasFileName;
        this.geneTrackName = geneTrackName;
        this.sequenceLocation = sequenceLocation;
        this.chromosomesAreOrdered = chromosomesAreOrdered;
        this.fasta = fasta;
        this.fastaDirectory = fastaDirectory;

        if(fastaFileNameString != null) {
            fastaFileNames = fastaFileNameString.split(",");
        }

        // Fix for legacy .genome files
        if (sequenceLocation != null && sequenceLocation.startsWith("/")) {
            if (!(new File(sequenceLocation)).exists()) {
                String tryThis = sequenceLocation.replaceFirst("/", "");
                if ((new File(tryThis)).exists()) {
                    this.sequenceLocation = tryThis;
                }
            }
        }
    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    // Used to determine feature file type, really only extension is needed

    public String getGeneFileName() {
        return geneFileName;
    }

    public String getGeneTrackName() {
        return geneTrackName;
    }

    public String[] getFastaFileNames() {
        return fastaFileNames;
    }

    public abstract InputStream getCytoBandStream() throws IOException;

    public abstract InputStream getGeneStream() throws IOException;

    public abstract InputStream getChrAliasStream() throws IOException;

    /**
     * Setter provided vor unit tests.
     *
     * @param sequenceLocation
     */
    public void setSequenceLocation(String sequenceLocation) {
        this.sequenceLocation = sequenceLocation;
    }

    public String getSequenceLocation() {
        return sequenceLocation;
    }

    @Override
    public String toString() {
        return name;
    }

    private boolean isFileGZipFormat(String fileName) {

        if (fileName == null) {
            return false;
        }

        if (fileName.toLowerCase().endsWith(".gz")) {
            return true;
        } else {
            return false;
        }
    }


    public boolean isChromosomesAreOrdered() {
        return chromosomesAreOrdered;
    }

    public boolean isChrNamesAltered() {
        return chrNamesAltered;
    }

    public String getUrl() {
        return url;
    }

    public void setUrl(String url) {
        this.url = url;
    }

    public boolean isFasta() {
        return fasta;
    }

    public boolean hasCytobands() {
        return cytoBandFileName != null && cytoBandFileName.length() > 0;
    }

    public abstract void close();
}
