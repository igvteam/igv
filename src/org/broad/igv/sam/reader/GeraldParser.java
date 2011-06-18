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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.GeraldAlignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 */
public class GeraldParser implements AlignmentParser {

    private static Logger log = Logger.getLogger(GeraldParser.class);

    private static final int READ_NAME_COLUMN = 7;
    private static final int PASSING_FILTER_COLUMN = 21;
    private static final int QUALITIES_COLUMN = 9;
    private static final int SINGLE_READ_ALIGNMENT_SCORE_COLUMN = 15;
    private static final int PAIRED_READ_ALIGNMENT_SCORE_COLUMN = 16;
    public static final int ALIGNMENT_START_COLUMN = 12;
    public static final int MATE_CHROMOSOME_COLUMN = 17;
    public static final int MATE_OFFSET_COLUMN = 19;
    public static final int READ_COLUMN = 8;
    public static final int READ_NUMBER = 7;
    public static final int CHROMOSOME_COLUMN = 10;
    private static final int MATE_STRAND_COLUMN = 20;
    private static final int READ_STRAND_COLUMN = 13;
    private static final int REQUIRED_EXPORT_COLUMNS = PASSING_FILTER_COLUMN + 1;
    private String[] fields = new String[REQUIRED_EXPORT_COLUMNS];
    SolexaQualityConverter solexaToPhredQualityConverter = SolexaQualityConverter.getSingleton();
    Genome genome;
    /**
     * Byte typed variables for all normal bases.
     */
    public static final byte a = 'a';
    public static final byte c = 'c';
    public static final byte g = 'g';
    public static final byte t = 't';
    public static final byte A = 'A';
    public static final byte C = 'C';
    public static final byte G = 'G';
    public static final byte T = 'T';

    public GeraldParser()

    {
        genome = IGV.getInstance().getGenomeManager().getCurrentGenome();
    }

    public GeraldAlignment readNextRecord(AsciiLineReader reader) {
        String nextLine = null;
        try {
            nextLine = reader.readLine();
        } catch (IOException e) {
            log.error("Error reading line", e);
            return null;
        }
        if (nextLine == null) {
            return null;
        }
        return createGeraldAlignment(nextLine);
    }

    private GeraldAlignment createGeraldAlignment(String nextLine) {


        int nTokens = ParsingUtils.split(nextLine, fields, '\t');
        // TODO -- what to do if nTokens < 22?
        StringBuffer readName = new StringBuffer(20);
        for (int c = 0; c < READ_NAME_COLUMN; c++) {
            if (c > 0 && fields[c] != null && fields[c].length() > 0) {
                readName.append(':');
            }
            readName.append(fields[c]);
        }
        GeraldAlignment alignment = new GeraldAlignment(readName.toString());

        String readNumberString = fields[READ_NUMBER];
        boolean isPaired = readNumberString != null && readNumberString.length() > 0 && nTokens > 18;

        String mappingQString = null;
        if (isPaired) {
            mappingQString = fields[PAIRED_READ_ALIGNMENT_SCORE_COLUMN];
        }
        if (mappingQString == null || mappingQString.length() == 0) {
            mappingQString = fields[SINGLE_READ_ALIGNMENT_SCORE_COLUMN];
        }
        if (mappingQString != null && mappingQString.length() > 0) {
            alignment.setMappingQuality(Integer.parseInt(mappingQString));
        } else {
            alignment.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY);
        }


        String readChr = mapChromosome(fields[CHROMOSOME_COLUMN]);
        int alignmentStart = Integer.parseInt(fields[ALIGNMENT_START_COLUMN]) - 1;

        String strandString = fields[READ_STRAND_COLUMN];
        boolean isNegativeStrand = strandString.equals("R");
        alignment.setNegativeStrand(isNegativeStrand);

        String readString = fields[READ_COLUMN];
        byte[] phredQualities = fields[QUALITIES_COLUMN].getBytes();
        solexaToPhredQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(phredQualities);
        if (isNegativeStrand) {
            readString = reverseComplement(readString);
            reverseArray(phredQualities);
        }
        byte[] read = readString.getBytes();
        alignment.setReads(readChr, alignmentStart, read, phredQualities);

        if (isPaired) {
            String mateChr = mapChromosome(fields[MATE_CHROMOSOME_COLUMN]);
            int mateStart = -1;
            if (mateChr == null || mateChr.length() == 0) {
                mateChr = readChr;
                String mateOffsetString = fields[MATE_OFFSET_COLUMN].trim();
                if (mateOffsetString != null && mateOffsetString.length() > 0) {
                    int inferredInsertSize = Integer.parseInt(mateOffsetString);
                    alignment.setInferredInsertSize(inferredInsertSize);
                    mateStart = alignment.getStart() + inferredInsertSize;
                }
            }
            String mateStrandString = fields[MATE_STRAND_COLUMN];
            boolean mateStrandIsNegative = mateStrandString.equals("R");
            alignment.setMate(new ReadMate(mateChr, mateStart, mateStrandIsNegative, true));
        }

        String passedFilterString = fields[PASSING_FILTER_COLUMN];
        boolean passingFilter = passedFilterString == null || passedFilterString.equals("Y");
        alignment.setPassedFilter(passingFilter);

        return alignment;
    }


    // TODO -- mapping

    public String mapChromosome(String seq) {

        if (seqChrMap == null) {
            loadChrMap();
        }
        String chr = seqChrMap.containsKey(seq) ? seqChrMap.get(seq) : seq;
        return genome == null ? chr : genome.getChromosomeAlias(chr);

    }

    static Map<String, String> seqChrMap;

    // TODO -- move this ot a utility class,  "loadMap(map, file)"

    static synchronized void loadChrMap() {
        seqChrMap = new HashMap();
        File samDir = new File(Globals.getIgvDirectory(), "sam");
        if (samDir.exists()) {
            File mapFile = new File(samDir, "sequence.map");
            if (!mapFile.exists()) {
                mapFile = new File(samDir, "sequence.map.txt");
            }
            if (mapFile.exists()) {
                BufferedReader br = null;
                try {
                    br = new BufferedReader(new FileReader(mapFile));
                    String nextLine;
                    while ((nextLine = br.readLine()) != null) {
                        String[] tokens = nextLine.split("\t");
                        if (tokens.length > 1) {
                            seqChrMap.put(tokens[0], tokens[1]);
                        }
                    }
                } catch (IOException e) {
                } finally {
                    if (br != null) {
                        try {
                            br.close();
                        } catch (IOException ex) {
                        }
                    }
                }
            }
        }

    }


    /**
     * clone the above method as necessary for non-object types
     */
    public static void reverseArray(byte[] array) {
        for (int left = 0, right = array.length - 1; left < right; left++, right--) {
            // exchange the first and last
            byte temp = array[left];
            array[left] = array[right];
            array[right] = temp;
        }
    }

    /**
     * Calculate the reverse complement of the specified sequence
     * (Stolen from Reseq)
     *
     * @param sequenceData
     * @return reverse complement
     */
    public static String reverseComplement(final String sequenceData) {
        final byte[] bases = net.sf.samtools.util.StringUtil.stringToBytes(sequenceData);
        reverseComplement(bases);
        return net.sf.samtools.util.StringUtil.bytesToString(bases);
    }

    /**
     * Returns the complement of a single byte.
     */
    public static final byte complement(final byte b) {
        switch (b) {
            case a:
                return t;
            case c:
                return g;
            case g:
                return c;
            case t:
                return a;
            case A:
                return T;
            case C:
                return G;
            case G:
                return C;
            case T:
                return A;
            default:
                return b;
        }
    }

    /**
     * Reverses and complements the bases in place.
     */
    public static void reverseComplement(final byte[] bases) {
        final int lastIndex = bases.length - 1;

        int i, j;
        for (i = 0, j = lastIndex; i < j; ++i, --j) {
            final byte tmp = complement(bases[i]);
            bases[i] = complement(bases[j]);
            bases[j] = tmp;
        }
        if (bases.length % 2 == 1) {
            bases[i] = complement(bases[i]);
        }
    }
}
