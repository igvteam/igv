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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.AlignmentUtils;
import org.broad.igv.sam.GeraldAlignment;
import org.broad.igv.sam.ReadMate;
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
    SolexaQualityConverter solexaToPhredQualityConverter = SolexaQualityConverter.getSingleton();
    //Genome genome;
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

    Genome genome;

    public GeraldParser() {
        //this.genome = genome;
        //   genome = GenomeManager.getInstance().getCurrentGenome();
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

        String[] fields = Globals.tabPattern.split(nextLine, -1);
        int nTokens = fields.length;
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
            readString = AlignmentUtils.reverseComplement(readString);
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

        String passedFilterString = null;
        if (fields.length > PASSING_FILTER_COLUMN) {
            passedFilterString = fields[PASSING_FILTER_COLUMN];
        }

        boolean passingFilter = passedFilterString == null || "Y".equals(passedFilterString);
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
        File samDir = DirectoryManager.getSamDirectory();
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

}
