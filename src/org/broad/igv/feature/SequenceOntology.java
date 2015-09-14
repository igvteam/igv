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

package org.broad.igv.feature;

import java.util.HashSet;
import java.util.Set;

/**
 * A static class for defining identifiers from the Sequence Ontology project and related efforts around GFF, EMBL, and
 * NCBI formats.
 * <p/>
 * Reference:    http://www.sequenceontology.org/
 */
public class SequenceOntology {


    public static Set<String> fivePrimeUTRTypes = new HashSet<String>();
    public static Set<String> threePrimeUTRTypes = new HashSet<String>();
    public static Set<String> utrTypes = new HashSet<String>();
    public static Set<String> cdsTypes = new HashSet<String>();
    public static Set<String> exonTypes = new HashSet<String>();
    public static Set<String> transcriptParts = new HashSet();
    public static Set<String> mrnaParts = new HashSet();
    public static Set<String> geneParts = new HashSet();


    static {
        fivePrimeUTRTypes.add("five_prime_UTR");
        fivePrimeUTRTypes.add("5'-UTR");
        fivePrimeUTRTypes.add("5'-utr");
        fivePrimeUTRTypes.add("5UTR");

        threePrimeUTRTypes.add("three_prime_UTR");
        threePrimeUTRTypes.add("3'-utr");
        threePrimeUTRTypes.add("3'-UTR");
        threePrimeUTRTypes.add("3UTR");

        utrTypes.addAll(SequenceOntology.fivePrimeUTRTypes);
        utrTypes.addAll(SequenceOntology.threePrimeUTRTypes);
        cdsTypes.add("CDS");
        cdsTypes.add("cds");

        exonTypes.add("exon");
        exonTypes.add("coding_exon");

        mrnaParts.addAll(utrTypes);
        mrnaParts.addAll(cdsTypes);
        mrnaParts.addAll(exonTypes);

        transcriptParts.addAll(mrnaParts);
        transcriptParts.add("intron");

        geneParts.addAll(transcriptParts);
        geneParts.add("transcript");
        geneParts.add("processed_transcript");
        geneParts.add("mrna");
        geneParts.add("mRNA");

    }

    public static boolean isCoding(String type) {
        return cdsTypes.contains(type) || type.equals("coding_exon");
    }
}
