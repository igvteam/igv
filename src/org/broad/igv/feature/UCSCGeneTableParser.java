/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.feature;

//~--- JDK imports ------------------------------------------------------------

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;

import java.util.List;

/**
 * @author jrobinso
 */
public class UCSCGeneTableParser extends UCSCParser {

    private int cdEndColumn = 7;
    private int cdStartColumn = 6;
    private int chrColumn = 2;
    private int endColumn = 5;
    private int endsBufferColumn = 10;
    private int exonCountColumn = 8;
    private int idColumn = 1;
    private int nameColumn = 0;
    private int startColumn = 4;
    private int startsBufferColumn = 9;
    private int strandColumn = 3;

    public enum Type {

        REFFLAT, GENEPRED, UCSCGENE
    }

    ;
    private Type type;

    /**
     * Constructs ...
     *
     * @param type
     */
    public UCSCGeneTableParser(Genome genome, Type type) {
        super(genome);
        
        this.type = type;
        switch (type) {
            case REFFLAT:
                break;
            case UCSCGENE:
                idColumn = 0;
                chrColumn = 1;
                strandColumn = 2;
                startColumn = 3;
                endColumn = 4;
                cdStartColumn = 5;
                cdEndColumn = 6;
                exonCountColumn = 7;
                startsBufferColumn = 8;
                endsBufferColumn = 9;
                nameColumn = 10;
                break;
            case GENEPRED:
                nameColumn = 12;
                break;

            //#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID
            //  0      1       2        3      4        5          6         7          8          9           10
        }
    }

    @Override
    protected IGVFeature parseLine(String[] tokens, int tokenCount) {


        if (tokenCount <= strandColumn) {
            return null;
        }

        String identifier = new String(tokens[idColumn].trim());
        String name = null;
        if (tokenCount > nameColumn && tokens[nameColumn] != null) {
            name = new String(tokens[nameColumn]);
        }

        if (name == null || name.length() == nameColumn) {
            name = identifier;
        }

        String chrToken = tokens[chrColumn].trim();
        String chr = genome == null ? chrToken : genome.getChromosomeAlias(chrToken);
        int start = Integer.parseInt(tokens[startColumn]) - startBase;
        int end = Integer.parseInt(tokens[endColumn]);
        String strandString = tokens[strandColumn];
        Strand strand = Strand.NONE;
        if (strandString != null) {
            if (strandString.trim().equals("+")) {
                strand = Strand.POSITIVE;
            } else if (strandString.trim().equals("-")) {
                strand = Strand.NEGATIVE;
            }
        }

        BasicFeature gene = new BasicFeature(chr, start, end, strand);

        gene.setName(name);
        gene.setIdentifier(identifier);


        // Coding information is optional
        if (tokenCount > 8) {
            createExons(tokens, tokenCount, gene, chr, strand);
        }
        return gene;
    }

    private void createExons(String[] tokens, int tokenCount, BasicFeature gene, String chr,
                             Strand strand)
            throws NumberFormatException {

        int cdStart = Integer.parseInt(tokens[cdStartColumn])  - startBase;
        int cdEnd = Integer.parseInt(tokens[cdEndColumn]);

        int exonCount = Integer.parseInt(tokens[exonCountColumn]);
        String[] startsBuffer = new String[exonCount];
        String[] endsBuffer = new String[exonCount];
        ParsingUtils.split(tokens[startsBufferColumn], startsBuffer, ',');
        ParsingUtils.split(tokens[endsBufferColumn], endsBuffer, ',');


        if (startsBuffer.length == endsBuffer.length) {
            int exonNumber = (strand == Strand.NEGATIVE ? exonCount : 1);
            for (int i = 0; i < startsBuffer.length; i++) {
                int exonStart = Integer.parseInt(startsBuffer[i])  - startBase;
                int exonEnd = Integer.parseInt(endsBuffer[i]);
                Exon exon = new Exon(chr, exonStart, exonEnd, strand);
                exon.setCodingStart(cdStart);
                exon.setCodingEnd(cdEnd);
                exon.setNumber(exonNumber);
                gene.addExon(exon);
                if (strand == Strand.NEGATIVE) {
                    exonNumber--;
                } else {
                    exonNumber++;
                }
            }
        }

        List<Exon> exons = gene.getExons();

        // Walk through exons setting mRNA start position
        int start = strand == Strand.POSITIVE ? 0 : exonCount - 1;
        int end = strand == Strand.POSITIVE ? exonCount : -1;
        int increment = strand == Strand.POSITIVE ? 1 : -1;

        int mRNABase = 0;
        for (int i = start; i != end; i += increment) {
            Exon exon = exons.get(i);
            exon.setMrnaBase(mRNABase);
            mRNABase += exon.getCodingLength();
        }


        if (type == Type.GENEPRED && tokenCount > 15) {
            try {

                String[] frameBuffer = new String[exonCount];
                ParsingUtils.split(tokens[15], frameBuffer, ',');
                for (int i = 0; i < frameBuffer.length; i++) {
                    int exonFrame = Integer.parseInt(frameBuffer[i].trim());
                    if (exonFrame == -1) {
                        exons.get(i).setUTR(true);
                    } else {
                        int phase = (exonFrame == 0) ? 0 : strandColumn - exonFrame;
                        exons.get(i).setPhase(phase);
                    }
                }
            } catch (Exception e) {

                // Ignore -- not getting the reading frame is not the end of the world.
                //log.info("Could not parse frame buffer: ", e);
            }

        }
    }
}
