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

package org.broad.igv.feature.tribble;

import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.StringUtils;

import java.util.List;

/**
 * @author Jim Robinson
 * @date 11/19/11
 */
public class UCSCGeneTableCodec extends UCSCCodec<BasicFeature> {

    private int nameColumn = 0;
    private int idColumn = 1;
    private int chrColumn = 2;
    private int strandColumn = 3;
    private int startColumn = 4;
    private int endColumn = 5;
    private int cdStartColumn = 6;
    private int cdEndColumn = 7;
    private int exonCountColumn = 8;
    private int startsBufferColumn = 9;
    private int endsBufferColumn = 10;
    private int frameBufferColumn = 15;

    public enum Type {

        REFFLAT, GENEPRED, UCSCGENE
    }

    private Genome genome;
    private Type type;


    public UCSCGeneTableCodec(Type type, Genome genome) {
        super(BasicFeature.class);

        this.genome = genome;
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
        }

    }

    /**
     * Decode a line as a Feature.
     *
     * @param line the input line to decode
     * @return Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    public BasicFeature decode(String line) {
        if (line.startsWith("#")) {
            //Header line
            readHeaderLine(line);
            return null;
        }

        line = line.replaceAll("\"", "");
        String[] tokens = Globals.singleTabMultiSpacePattern.split(line);
        int tokenCount = tokens.length;

        if (tokenCount <= strandColumn) {
            return null;
        }

        String identifier = tokens[idColumn].trim();
        String name = null;
        if (tokenCount > nameColumn && tokens[nameColumn] != null) {
            name = tokens[nameColumn];
        }

        if (name == null || name.length() == nameColumn) {
            name = identifier;
        }

        String chrToken = tokens[chrColumn].trim();
        String chr = genome == null ? StringUtils.intern(chrToken) : genome.getChromosomeAlias(chrToken);

        int start = Integer.parseInt(tokens[startColumn]);
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

        if (tokenCount > 7) {
            gene.setThickStart(Integer.parseInt(tokens[6]));
            gene.setThickEnd(Integer.parseInt(tokens[7]));
        }

        // Coding information is optional
        if (tokenCount > 8) {
            createExons(tokens, tokenCount, gene, chr, strand);
        }
        return gene;
    }

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(String path) {
        return true;
    }


    private void createExons(String[] tokens, int tokenCount, BasicFeature gene, String chr,
                             Strand strand)
            throws NumberFormatException {

        int cdStart = Integer.parseInt(tokens[cdStartColumn]);
        int cdEnd = Integer.parseInt(tokens[cdEndColumn]);

        int exonCount = Integer.parseInt(tokens[exonCountColumn]);
        String[] startsBuffer = Globals.commaPattern.split(tokens[startsBufferColumn]);
        String[] endsBuffer = Globals.commaPattern.split(tokens[endsBufferColumn]);

        if (startsBuffer.length == endsBuffer.length) {
            int exonNumber = (strand == Strand.NEGATIVE ? exonCount : 1);
            for (int i = 0; i < startsBuffer.length; i++) {
                int exonStart = Integer.parseInt(startsBuffer[i]);
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

                String[] frameBuffer = Globals.commaPattern.split(tokens[frameBufferColumn]);
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
