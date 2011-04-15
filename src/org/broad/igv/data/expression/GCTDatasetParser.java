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
/*
 * GCTDatasetParser.java
 *
 * Parser for GCT and related file formats (e.g. RES).
 *
 * Created on October 18, 2007, 2:33 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.data.expression;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.collections.IntArrayList;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.exceptions.ProbeMappingException;
import org.broad.igv.feature.*;
import org.broad.igv.tools.StatusMonitor;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * TODO -- handle case with probe file
 *
 * @author jrobinso
 */
public class GCTDatasetParser {

    private static Logger log = Logger.getLogger(GCTDatasetParser.class);
    static String LOCUS_START_DELIMITER = "|@";
    static String LOCUS_END_DELIMITER = "|";

    enum FileType {

        RES, GCT, MAPPED, TAB, MET, DCHIP, MAGE_TAB
    }

    ResourceLocator dataFileLocator;

    FileType type;
    int dataStartColumn;
    int probeColumn;
    int descriptionColumn;
    Genome genome;

    /**
     * Map chr -> longest feature
     */
    Map<String, Integer> longestProbeMap;

    /**
     * Map colum heading -> index for effecient reverse lookup
     */
    private Map<String, Integer> headingIndexMap = new HashMap();
    Map<String, List<Row>> rowMap = new HashMap();
    StatusMonitor statusMonitor;
    GeneManager geneManager = null;
    GeneToLocusHelper locusHelper;


    /**
     * Constructs ...
     *
     * @param resFile
     * @param probeFile
     * @param genome
     */
    public GCTDatasetParser(ResourceLocator resFile, String probeFile, Genome genome) throws IOException {

        this.dataFileLocator = resFile;
        this.genome = genome;
        this.geneManager = GeneManager.getGeneManager(genome.getId());
        longestProbeMap = new HashMap();
        locusHelper = new GeneToLocusHelper(probeFile, genome.getId());

    }


    /**
     * Constructs ...
     *
     * @param resFile
     * @param probeFile
     * @param genome
     */
    public GCTDatasetParser(File resFile, String probeFile, Genome genome) throws IOException {
        this(new ResourceLocator(resFile.getAbsolutePath()), probeFile, genome);
    }

    /*
     */

    public static boolean parsableMAGE_TAB(ResourceLocator file) throws IOException {
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(file);
            String nextLine = null;

            //skip first row
            reader.readLine();

            //check second row for MAGE_TAB identifiers
            if ((nextLine = reader.readLine()) != null && (nextLine.contains("Reporter REF") || nextLine.contains("Composite Element REF") || nextLine.contains("Term Source REF") || nextLine.contains("CompositeElement REF") || nextLine.contains("TermSource REF") || nextLine.contains("Coordinates REF"))) {
                int count = 0;
                // check if this mage_tab data matrix can be parsed by this class
                while ((nextLine = reader.readLine()) != null && count < 5) {
                    nextLine = nextLine.trim();
                    if (nextLine.startsWith("SNP_A") || nextLine.startsWith("CN_")) {
                        return false;
                    }

                    count++;
                }
                return true;
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        return false;
    }

    public GCTDataset createDataset() {
        GCTDataset dataset = new GCTDataset(genome);
        parse(dataset);
        return dataset;
    }

    /**
     * Parse the file and return a Dataset
     *
     * @return
     */
    public void parse(GCTDataset dataset) {

        // Create a buffer for the string split utility.  We use  a custom utility as opposed
        // to String.split() for performance.
        String[] tokens = new String[100000];

        String fn = dataFileLocator.getPath().toLowerCase();
        if (fn.endsWith(".txt") || fn.endsWith(".tab") || fn.endsWith(".xls") || fn.endsWith(".gz")) {
            fn = fn.substring(0, fn.lastIndexOf("."));
        }


        descriptionColumn = -1;    // Default - no description column
        if (fn.endsWith("res")) {
            type = FileType.RES;
            dataStartColumn = 2;
            probeColumn = 1;
            descriptionColumn = 0;
        } else if (fn.endsWith("gct")) {
            type = FileType.GCT;
            dataStartColumn = 2;
            probeColumn = 0;
            descriptionColumn = 1;
        } else if (fn.endsWith("mapped")) {
            type = FileType.MAPPED;
            dataStartColumn = 4;
            probeColumn = 0;
        } else if (fn.endsWith("met")) {
            type = FileType.MET;
            dataStartColumn = 4;
            probeColumn = 0;
        } else if (fn.endsWith("dchip")) {
            type = FileType.DCHIP;
            dataStartColumn = 1;
            probeColumn = 0;
            descriptionColumn = -1;
        } else if (dataFileLocator.getDescription() != null &&
                dataFileLocator.getDescription().equals("MAGE_TAB")) {
            type = FileType.MAGE_TAB;
            descriptionColumn = -1;
            dataStartColumn = 1;
            probeColumn = 0;
        } else {
            type = FileType.TAB;
            dataStartColumn = 1;
            probeColumn = 0;
        }


        boolean hasCalls = (type == FileType.RES);
        boolean hasDescription = (descriptionColumn >= 0);

        //
        //Assume data is gene expression for now
        dataset.setType(TrackType.GENE_EXPRESSION);

        AsciiLineReader reader = null;
        String nextLine = null;
        String[] columnHeadings = null;
        try {
            reader = ParsingUtils.openAsciiReader(dataFileLocator);

            String headerLine = null;

            // Skip header rows
            if (type == FileType.GCT) {
                nextLine = reader.readLine();
                if (nextLine.startsWith("#")) {
                    parseComment(nextLine, dataset);
                }
                nextLine = reader.readLine();
                if (nextLine.startsWith("#")) {
                    parseComment(nextLine, dataset);
                }
                headerLine = reader.readLine();
            } else if (type != FileType.MAGE_TAB) {

                // Skip meta data, if any
                while ((nextLine = reader.readLine()).startsWith("#") && (nextLine != null)) {
                    parseComment(nextLine, dataset);
                }
                headerLine = nextLine;
            } else {
                headerLine = reader.readLine();
            }

            // Parse column headings
            int skip = hasCalls ? 2 : 1;
            int nTokens = ParsingUtils.split(headerLine, tokens, '\t');

            int nColumns = (nTokens - dataStartColumn) / skip;
            ArrayList columnHeadingsObj = new ArrayList();
            for (int i = 0; i < nColumns; i++) {
                String heading = tokens[dataStartColumn + i * skip].replace('\"', ' ').trim();
                if (type == FileType.MAGE_TAB) {
                    if (!columnHeadingsObj.contains(heading)) {
                        columnHeadingsObj.add(heading);
                        headingIndexMap.put(heading, columnHeadingsObj.size() - 1);
                    }
                } else {
                    columnHeadingsObj.add(heading);
                    headingIndexMap.put(heading, i);
                }
            }

            columnHeadings = (String[]) columnHeadingsObj.toArray(new String[0]);
            dataset.setColumnHeadings(columnHeadings);

            nColumns = columnHeadings.length;

            //parse quantitation type column header

            IntArrayList valuesIndices = new IntArrayList(nColumns);
            if (type == FileType.MAGE_TAB) {
                nextLine = reader.readLine();
                nTokens = ParsingUtils.split(nextLine, tokens, '\t');
                for (int i = dataStartColumn; i < nTokens; i++) {
                    String heading = tokens[i].replace('\"', ' ').trim();

                    //Check for tcga data column headings
                    if (heading.contains("Beta_Value") || heading.contains("Beta value") || heading.contains("log2 Signal") || heading.contains("Signal") || heading.contains("unc_DWD_Batch_adjusted")) {
                        valuesIndices.add(i);
                    }
                    if (heading.contains("Gene symbol") || heading.contains("Gene_Symbol")) {
                        descriptionColumn = i;
                        hasDescription = true;
                    }
                }

                if (nColumns != valuesIndices.size()) {
                    // Throw parser error
                }
            }

            // If format is RES skip the two lines following the header
            if (type == FileType.RES) {
                reader.readLine();
                reader.readLine();
            }

            int lineCount = 0;

            while ((nextLine = reader.readLine()) != null) {
                nTokens = ParsingUtils.split(nextLine, tokens, '\t');
                String probeId = new String(tokens[probeColumn]);
                float[] values = new float[nColumns];
                char[] calls = hasCalls ? new char[nColumns] : (char[]) null;

                String description = (hasDescription && (nTokens > descriptionColumn))
                        ? new String("|@" + tokens[descriptionColumn] + "|") : null;

                if (type == FileType.MAGE_TAB && probeId.startsWith("cg")) {
                    dataset.setType(TrackType.DNA_METHYLATION);
                }

                for (int i = 0; i < nColumns; i++) {
                    try {
                        int dataIndex = -1;
                        if (type == FileType.MAGE_TAB) {
                            dataIndex = valuesIndices.get(i);
                        } else {
                            dataIndex = dataStartColumn + i * skip;
                        }

                        // If we are out of value tokens, or the cell is blank, assign NAN to the cell.
                        if ((dataIndex >= nTokens) || (tokens[dataIndex].length() == 0)) {
                            values[i] = Float.NaN;
                        } else {
                            values[i] = Float.parseFloat(tokens[dataIndex]);
                        }

                    } catch (NumberFormatException numberFormatException) {

                        // This s an expected condition.  IGV uses NaN to
                        // indicate non numbers (missing data values)
                        values[i] = Float.NaN;
                    }

                    // We ignore calls, just skip them if present
                    if (hasCalls) {
                        calls[i] = tokens[3 + i * skip].charAt(0);
                    }
                }
                addRow(probeId, description, values, calls);
                lineCount++;

                // This method is designed to be interruptable (canceled by
                // user.  Check every 1000 lines for an interrupt.
                if (lineCount == 1000) {
                    checkForInterrupt();
                    lineCount = 0;
                    if (statusMonitor != null) {
                        statusMonitor.incrementPercentComplete(1);
                    }
                }
            }    // End loop through lines

        } catch (FileNotFoundException ex) {
            throw new RuntimeException(ex);
        } catch (InterruptedException e) {
            throw new RuntimeException("Operation cancelled");
        } catch (Exception e) {
            if (nextLine != null && reader.getCurrentLineNumber() != 0) {
                throw new ParserException(e.getMessage(), e, reader.getCurrentLineNumber(), nextLine);
            } else {
                throw new RuntimeException(e);
            }
        }
        finally {
            if (reader != null) {
                reader.close();
            }
        }


        // Sort row for each chromosome by start location
        sortRows();

        // Update dataset
        for (String chr : rowMap.keySet()) {
            dataset.setStartLocations(chr, getStartLocations(chr));
            dataset.setEndLocations(chr, getEndLocations(chr));
            dataset.setFeatureNames(chr, getProbes(chr));
            for (String heading : columnHeadings) {
                dataset.setData(heading, chr, getData(heading, chr));
            }
        }

        dataset.setLongestFeatureMap(longestProbeMap);

        if ((dataset == null) || dataset.isEmpty()) {
            String genomeId = genome == null ? "" : genome.getId();
            throw new ProbeMappingException(fn, genomeId);
        }
    }


    public void addRow(String probeId, String description, float[] values, char[] calls) {

        List<Locus> loci = locusHelper.getLoci(probeId, description);
        System.out.println("description: " + description);
        System.out.println("loci" + loci);
        if (loci != null) {
            for (Locus locus : loci) {
                if ((locus != null) && locus.isValid()) {
                    addRow(probeId, locus, values);
                }
            }
        }
    }

    private void addRow(String probeId, Locus locus, float[] values) {

        //if (log.isDebugEnabled()) {
        //    log.debug("Adding row for probe " + probeId + " to chromosome " + locus.getChr());
        //}

        List<Row> rows = rowMap.get(locus.getChr());
        if (rows == null) {
            rows = new ArrayList();
            rowMap.put(locus.getChr(), rows);
        }

        String chr = locus.getChr();
        int length = locus.getEnd() - locus.getStart();
        if (longestProbeMap.containsKey(chr)) {
            longestProbeMap.put(chr, Math.max(longestProbeMap.get(chr), length));
        } else {
            longestProbeMap.put(chr, length);
        }

        rows.add(new Row(probeId, locus.getChr(), locus.getStart(), locus.getEnd(), values));

    }

    /**
     * Sort all row collections by ascending start location
     */
    private void sortRows() {
        Comparator<Row> c = new Comparator<Row>() {

            public int compare(GCTDatasetParser.Row arg0, GCTDatasetParser.Row arg1) {
                return arg0.start - arg1.start;
            }
        };
        for (List<Row> rows : rowMap.values()) {
            Collections.sort(rows, c);
        }
    }

    public String[] getProbes(String chr) {
        List<Row> rows = rowMap.get(chr);
        String[] labels = new String[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            labels[i] = rows.get(i).feature;
        }
        return labels;

    }


    public int[] getStartLocations(String chr) {

        List<Row> rows = rowMap.get(chr);
        int[] startLocations = new int[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            startLocations[i] = rows.get(i).start;
        }
        return startLocations;
    }


    public int[] getEndLocations(String chr) {

        List<Row> rows = rowMap.get(chr);
        int[] endLocations = new int[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            endLocations[i] = rows.get(i).end;
        }
        return endLocations;
    }


    public float[] getData(String heading, String chr) {

        int columnIndex = this.headingIndexMap.get(heading);
        List<Row> rows = rowMap.get(chr);

        float[] data = new float[rows.size()];
        for (int i = 0; i < rows.size(); i++) {
            data[i] = rows.get(i).values[columnIndex];
        }
        return data;

    }


    /**
     * Note:  This is an exact copy of the method in IGVDatasetParser.  Refactor to merge these
     * two parsers, or share a common base class.
     *
     * @param comment
     * @param dataset
     */
    private void parseComment(String comment, GCTDataConsumer dataset) {

        String tmp = comment.substring(1, comment.length());
        if (tmp.startsWith("track")) {
            dataset.setTrackLine(tmp);

        } else {
            String[] tokens = tmp.split("=");
            if (tokens.length != 2) {
                return;
            }
            String key = tokens[0].trim().toLowerCase();
            if (key.equals("name")) {
                dataset.setName(tokens[1].trim());
            } else if (key.equals("type")) {

                try {
                    dataset.setType(TrackType.valueOf(tokens[1].trim().toUpperCase()));
                } catch (Exception exception) {
                    // Ignore
                }
            }
        }
    }

    private void checkForInterrupt() throws InterruptedException {
        Thread.sleep(1);    // <- check for interrupted thread
    }

    class Row {

        String feature;
        String chr;
        int start;
        int end;
        float[] values;

        Row(String feature, String chr, int start, int end, float[] values) {
            this.feature = feature;
            this.chr = chr;
            this.start = start;
            this.end = end;
            this.values = values;
        }
    }
}
