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

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.collections.MultiMap;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * Codec for .mut and .maf mutation files
 *
 * @author jrobinso
 *         Date: 10/24/12
 *         Time: 12:05 PM
 */
public class MUTCodec extends AsciiFeatureCodec<Mutation> {

    private static Logger log = Logger.getLogger(MUTCodec.class);

    private String path;  // for error messages
    private boolean isMAF;
    private String[] headers;
    private String[] samples;
    private Genome genome;
    private int chrColumn;
    private int startColumn;
    private int endColumn;
    private int sampleColumn;
    private int typeColumn;
    private int refAlleleColumn;
    private int tumorAllele1Column;
    private int tumorAllele2Column;


    public MUTCodec(String path, Genome genome) {
        super(Mutation.class);
        this.path = path;
        this.genome = genome;
        try {
            LineIterator reader = new LineIteratorImpl(new AsciiLineReader(ParsingUtils.openInputStream(path)));
            readActualHeader(reader);
        } catch (IOException e) {
            log.error(e.getMessage(), e);
        }
    }

    public Void readActualHeader(LineIterator reader) {
        String nextLine = null;

        while (reader.hasNext()) {
            nextLine = reader.peek();
            if (nextLine.startsWith("#")) {
                reader.next();
                if (nextLine.startsWith("#samples")) {
                    String[] tokens = Globals.whitespacePattern.split(nextLine, 2);
                    if (tokens.length < 2) {
                        log.error("Error parsing sample header in mutation file: " + path);
                    } else {
                        samples = Globals.commaPattern.split(tokens[1]);
                        for (int i = 0; i < samples.length; i++) {
                            samples[i] = samples[i].trim();
                        }
                    }
                }
                continue;
            }

            String[] tokens = nextLine.split("\t");
            if (tokens.length > 4) {

                headers = Globals.tabPattern.split(nextLine);
                isMAF = headers.length > 15 && headers[0].equalsIgnoreCase("Hugo_Symbol");
                setColumns(isMAF);
                return null;
            }
        }
        throw new RuntimeException("Unexpected end-of-file (no header line): " + path);
    }

    public String[] getSamples() {
        return samples;
    }

    public int getChrColumn() {
        return chrColumn;
    }

    public int getStartColumn() {
        return startColumn;
    }

    @Override
    public Mutation decode(String line) {

        if(line.startsWith("#") || line.startsWith("Hugo_Symbol")) return null;

        String[] tokens = Globals.tabPattern.split(line);

        String chr = genome == null ? tokens[chrColumn].trim() : genome.getChromosomeAlias(tokens[chrColumn].trim());

        int start;
        try {
            start = Integer.parseInt(tokens[startColumn].trim());
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (startColumn + 1) + " must be a numeric value.", path);
        }

        int end;
        try {
            end = Integer.parseInt(tokens[endColumn].trim());
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (endColumn + 1) + " must be a numeric value.", path);
        }


        // MAF files use the 1-based inclusive convention for coordinates.  The convention is not
        // specified for MUT files, and it appears both conventions have been used.  We can detect
        // the convention used for single base mutations by testing start == end.
        if (isMAF || (start == end)) {
            start--;
        }

        String sampleId = tokens[sampleColumn].trim();
        String type = tokens[typeColumn].trim();

        MultiMap<String, String> attributes = new MultiMap();
        int n = Math.min(headers.length, tokens.length);
        for (int i = 0; i < n; i++) {
            String key = headers[i];
            String value = tokens[i];
            if (value.length() > 0) {
                attributes.put(key, value);
            }
        }


        Mutation mut = new Mutation(sampleId, chr, start, end, type);
        mut.setAttributes(attributes);

        if (refAlleleColumn > 0) {
            mut.setRefAllele(tokens[refAlleleColumn].trim());
        }
        if (tumorAllele1Column > 0) {
            mut.setAltAllele1(tokens[tumorAllele1Column].trim());
        }
        if (tumorAllele2Column > 0) {
            mut.setAltAllele2(tokens[tumorAllele2Column].trim());
        }

        return mut;
    }


    private void setColumns(boolean isMAF) {

        this.isMAF = isMAF;
        if (isMAF) {
            chrColumn = 4;
            startColumn = 5;
            endColumn = 6;
            sampleColumn = 15;
            typeColumn = 8;
            refAlleleColumn = 10;
            tumorAllele1Column = 11;
            tumorAllele2Column = 12;
        } else {
            chrColumn = 0;
            startColumn = 1;
            endColumn = 2;
            sampleColumn = 3;
            typeColumn = 4;
            refAlleleColumn = -1;
            tumorAllele1Column = -1;
            tumorAllele2Column = -1;
        }
    }

    /**
     * @param path
     * @return
     * @throws IOException
     */
    public static boolean isMutationAnnotationFile(String path) {

        String ext;
        {
            //Only want pathLC when checking extension, so we limit its scope
            String pathLC = path.toLowerCase();
            if(pathLC.endsWith(".maf.annotated")) {
                // TCGA extension
                return true;
            }

            ext = ParsingUtils.getIGVExtension(pathLC);
        }
        if (ext.equals("mut")) {
            return true;
        } else if(ext.equals("maf")) {

            BufferedReader reader = null;
            try {
                reader = ParsingUtils.openBufferedReader(path);
                if (reader == null) {
                    return false;
                }

                String nextLine;
                while ((nextLine = reader.readLine()) != null && nextLine.startsWith("#")) {
                    if (nextLine.startsWith("#")) {
                        continue;
                    }
                }

                String[] tokens = nextLine.split("\t");
                return tokens.length > 15 && tokens[0].equalsIgnoreCase("Hugo_Symbol");
            } catch (IOException e) {
                log.error("Error reading: " + path, e);
                return false;
            } finally {
                if (reader != null) {
                    try {
                        reader.close();
                    } catch (IOException e) {
                        log.error("Error closing: " + path, e);
                    }
                }
            }
        }
        else {
            return false;
        }


    }

}
