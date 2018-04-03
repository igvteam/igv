/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
 * Author: Jim Robinson
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

package org.broad.igv.tools.sort;

import org.apache.log4j.Logger;
import org.broad.igv.feature.tribble.MUTCodec;
import org.broad.igv.gwas.GWASParser;
import org.broad.igv.track.GFFFeatureSource;
import org.broad.igv.util.ResourceLocator;

import java.io.File;

/**
 * Created by jrobinson on 5/10/16.
 */
public class SorterFactory {

    private static Logger log = Logger.getLogger(SorterFactory.class);

    public static Sorter getSorter(File inputFile, File outputFile) {

        String shortFN = inputFile.getName().toLowerCase();
        if (shortFN.endsWith(".txt")) {
            shortFN = shortFN.substring(0, shortFN.length() - 4);
        }
        if (shortFN.endsWith(".cn") || shortFN.endsWith(".xcn") || shortFN.endsWith(".snp") || shortFN.endsWith(".igv")) {
            return new CNSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".sam")) {
            return new SAMSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".aligned") || shortFN.endsWith(".bed") || shortFN.endsWith(".bedgraph") || shortFN.endsWith(".bdg") ) {
            return new BedSorter(inputFile, outputFile);
        } else if (GFFFeatureSource.isGFF(shortFN)) {
            return new GFFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".vcf")) {
            return new VCFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".psl") || shortFN.endsWith(".pslx")) {
            return new BedSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".eqtl")) {
            return new EQTLSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".snp")) {
            return new GenericSorter(inputFile, outputFile, 1, 2);
        } else if (GWASParser.isGWASFile(shortFN)) {
            return new GWASSorter(inputFile, outputFile);
        } else if (MUTCodec.isMutationAnnotationFile(new ResourceLocator(inputFile.getAbsolutePath()))) {
            return new MUTSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".interaction")) {
            return new InteractionSorter(inputFile, outputFile);
        } else if(shortFN.endsWith(".bam")) {
            return new BAMSorter(inputFile, outputFile);
        } else {
            log.error("Unknown file type or sorting not supported for: " + inputFile.getName());
            return null;
        }
    }
}
