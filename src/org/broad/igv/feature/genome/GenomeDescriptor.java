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
    private boolean hasCustomSequenceLocation;
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
                            boolean hasCustomSequenceLocation,
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
        this.hasCustomSequenceLocation = hasCustomSequenceLocation;
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

    public boolean hasCustomSequenceLocation() {
        return hasCustomSequenceLocation;
    }
}
