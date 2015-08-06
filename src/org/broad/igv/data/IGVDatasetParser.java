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

package org.broad.igv.data;

//~--- non-JDK imports --------------------------------------------------------

import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.FloatArrayList;
import org.broad.igv.util.collections.IntArrayList;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Class description
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/11/11
 */
public class IGVDatasetParser {

    private static Logger log = Logger.getLogger(IGVDatasetParser.class);
    private ResourceLocator dataResourceLocator;
    private int chrColumn = -1;
    private int startColumn = -1;
    private int endColumn = -1;
    private int firstDataColumn = -1;
    private int lastDataColumn = -1;
    private int probeColumn = -1;
    private boolean hasEndLocations = false;
    private boolean hasCalls = false;
    private Genome genome;
    private IGV igv;


    private int startBase = 0;

    public IGVDatasetParser(ResourceLocator copyNoFile, Genome genome) {
        this.dataResourceLocator = copyNoFile;
        this.genome = genome;
        this.igv = IGV.hasInstance() ? IGV.getInstance() : null;
    }

    private void setColumnDefaults() {
        String tmp = (dataResourceLocator.getPath().endsWith(".txt")
                ? dataResourceLocator.getPath().substring(0,
                dataResourceLocator.getPath().length() - 4) : dataResourceLocator.getPath()).toLowerCase();

        if (tmp.endsWith(".igv")) {
            chrColumn = 0;
            startColumn = 1;
            endColumn = 2;
            probeColumn = 3;
            firstDataColumn = 4;
            hasEndLocations = true;
            hasCalls = false;
        } else if (tmp.endsWith(".xcn") || tmp.endsWith("cn") || tmp.endsWith(".snp") || tmp.endsWith(".loh")) {
            probeColumn = 0;
            chrColumn = 1;
            startColumn = 2;
            endColumn = -1;
            firstDataColumn = 3;
            hasEndLocations = false;
            hasCalls = tmp.endsWith(".xcn") || tmp.endsWith(".snp");
        } else if (tmp.endsWith(".expr")) {
            //gene_id	bundle_id	chr	left	right	FPKM	FPKM_conf_lo	FPKM_conf_hi
            probeColumn = 0;
            chrColumn = 2;
            startColumn = 3;
            endColumn = 4;
            startBase = 1;
            firstDataColumn = 5;
            lastDataColumn = 5;
            hasEndLocations = true;

        } else {
            // TODO -- popup dialog and ask user to define columns,  and csv vs tsv?
            throw new ParserException("Unknown file type: ", 0);
        }
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
                        return true;
                    }

                    count++;
                }
                return false;
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        return false;
    }

    /**
     * Scan the datafile for chromosome breaks.
     *
     * @param dataset
     * @return
     */
    public List<ChromosomeSummary> scan(IGVDataset dataset) {

        int estLineCount = ParsingUtils.estimateLineCount(dataResourceLocator.getPath());
        Map<String, Integer> longestFeatureMap = new HashMap();

        float dataMin = 0;
        float dataMax = 0;

        InputStream is = null;
        AsciiLineReader reader = null;
        String nextLine = null;
        ChromosomeSummary chrSummary = null;
        List<ChromosomeSummary> chrSummaries = new ArrayList();
        String[] headings = null;
        WholeGenomeData wgData = null;
        int nRows = 0;

        int headerRows = 0;
        int count = 0;

        boolean logNormalized;
        try {

            int skipColumns = hasCalls ? 2 : 1;

            // BufferedReader reader = ParsingUtils.openBufferedReader(dataResourceLocator);
            is = ParsingUtils.openInputStreamGZ(dataResourceLocator);
            reader = new AsciiLineReader(is);

            // Infer datatype from extension.  This can be overriden in the
            // comment section
            if (isCopyNumberFileExt(dataResourceLocator.getPath())) {
                dataset.setTrackType(TrackType.COPY_NUMBER);
                dataset.getTrackProperties().setWindowingFunction(WindowFunction.mean);
            } else if (isLOHFileExt(dataResourceLocator.getPath())) {
                dataset.setTrackType(TrackType.LOH);
                dataset.getTrackProperties().setWindowingFunction(WindowFunction.mean);
            } else {
                dataset.getTrackProperties().setWindowingFunction(WindowFunction.mean);
            }

            // Parse comments and directives, if any
            nextLine = reader.readLine();
            while (nextLine.startsWith("#") || (nextLine.trim().length() == 0)) {
                headerRows++;

                if (nextLine.length() > 0) {
                    parseDirective(nextLine, dataset);
                }
                nextLine = reader.readLine();
            }

            if (chrColumn < 0) {
                setColumnDefaults();
            }


            // Parse column headings
            String[] data = nextLine.trim().split("\t");

            // Set last data column
            if (lastDataColumn < 0) {
                lastDataColumn = data.length - 1;
            }

            headings = getHeadings(data, skipColumns);

            dataset.setDataHeadings(headings);

            // Infer if the data is logNormalized by looking for negative data values.
            // Assume it is not until proven otherwise
            logNormalized = false;

            wgData = new WholeGenomeData(headings);

            int chrRowCount = 0;

            // Update
            int updateCount = 5000;
            long lastPosition = 0;
            while ((nextLine = reader.readLine()) != null) {

                if (igv != null && ++count % updateCount == 0) {
                    igv.setStatusBarMessage("Loaded: " + count + " / " + estLineCount + " (est)");
                }
                // Distance since last sample

                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                int nTokens = tokens.length;
                if (nTokens > 0) {
                    String thisChr = genome.getChromosomeAlias(tokens[chrColumn]);
                    if (chrSummary == null || !thisChr.equals(chrSummary.getName())) {
                        // Update whole genome and previous chromosome summary, unless this is
                        // the first chromosome
                        if (chrSummary != null) {
                            updateWholeGenome(chrSummary.getName(), dataset, headings, wgData);
                            chrSummary.setNDataPoints(nRows);
                        }

                        // Shart the next chromosome
                        chrSummary = new ChromosomeSummary(thisChr, lastPosition);
                        chrSummaries.add(chrSummary);
                        nRows = 0;
                        wgData = new WholeGenomeData(headings);
                        chrRowCount = 0;

                    }
                    lastPosition = reader.getPosition();

                    int location = -1;
                    try {
                        location = ParsingUtils.parseInt(tokens[startColumn]) - startBase;

                    } catch (NumberFormatException numberFormatException) {
                        log.error("Column " + tokens[startColumn] + " is not a number");
                        throw new ParserException("Column " + (startColumn + 1) +
                                " must contain an integer value." + " Found: " + tokens[startColumn],
                                count + headerRows, nextLine);
                    }

                    int length = 1;
                    if (hasEndLocations) {
                        try {
                            length = ParsingUtils.parseInt(tokens[endColumn].trim()) - location + 1;

                        } catch (NumberFormatException numberFormatException) {
                            log.error("Column " + tokens[endColumn] + " is not a number");
                            throw new ParserException("Column " + (endColumn + 1) +
                                    " must contain an integer value." + " Found: " + tokens[endColumn],
                                    count + headerRows, nextLine);
                        }
                    }

                    updateLongestFeature(longestFeatureMap, thisChr, length);

                    if (wgData.locations.size() > 0 && wgData.locations.get(wgData.locations.size() - 1) > location) {
                        throw new ParserException("File is not sorted, .igv and .cn files must be sorted by start position." +
                                " Use igvtools (File > Run igvtools..) to sort the file.", count + headerRows);
                    }

                    wgData.locations.add(location);

                    for (int idx = 0; idx < headings.length; idx++) {
                        int i = firstDataColumn + idx * skipColumns;

                        float copyNo = i < tokens.length ? readFloat(tokens[i]) : Float.NaN;

                        if (!Float.isNaN(copyNo)) {
                            dataMin = Math.min(dataMin, copyNo);
                            dataMax = Math.max(dataMax, copyNo);
                        }
                        if (copyNo < 0) {
                            logNormalized = true;
                        }
                        String heading = headings[idx];
                        wgData.data.get(heading).add(copyNo);
                    }

                    nRows++;

                }
                chrRowCount++;
            }

            dataset.setLongestFeatureMap(longestFeatureMap);

        } catch (ParserException pe) {
            throw pe;
        } catch (FileNotFoundException e) {
            // DialogUtils.showError("SNP file not found: " + dataSource.getCopyNoFile());
            log.error("File not found: " + dataResourceLocator);
            throw new RuntimeException(e);
        } catch (Exception e) {
            log.error("Exception when loading: " + dataResourceLocator.getPath(), e);
            if (nextLine != null && (count + headerRows != 0)) {
                throw new ParserException(e.getMessage(), e, count + headerRows, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    log.error("Error closing IGVDataset stream", e);
                }
            }
        }

        // Update last chromosome
        if (chrSummary != null) {
            updateWholeGenome(chrSummary.getName(), dataset, headings, wgData);
            chrSummary.setNDataPoints(nRows);
        }

        dataset.setLogNormalized(logNormalized);
        dataset.setDataMin(dataMin);
        dataset.setDataMax(dataMax);

        return chrSummaries;
    }

    private void updateLongestFeature(Map<String, Integer> longestFeatureMap, String thisChr, int length) {
        if (longestFeatureMap.containsKey(thisChr)) {
            longestFeatureMap.put(thisChr, Math.max(longestFeatureMap.get(thisChr), length));
        } else {
            longestFeatureMap.put(thisChr, length);
        }
    }

    private float readFloat(String token) {
        float copyNo = Float.NaN;
        try {
            if (token != null) {
                copyNo = Float.parseFloat(token);
            }

        } catch (NumberFormatException e) {
            // This is an expected condition.
        }
        return copyNo;

    }

    /**
     * Load data for a single chromosome.
     *
     * @param chrSummary
     * @param dataHeaders
     * @return
     */
    public ChromosomeData loadChromosomeData(ChromosomeSummary chrSummary, String[] dataHeaders) {

        // InputStream is = null;
        try {
            int skipColumns = hasCalls ? 2 : 1;

            // Get an estimate of the number of snps (rows).  THIS IS ONLY AN ESTIMATE
            int nRowsEst = chrSummary.getNDataPts();

            SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(dataResourceLocator.getPath());
            is.seek(chrSummary.getStartPosition());
            AsciiLineReader reader = new AsciiLineReader(is);


            // Create containers to hold data
            IntArrayList startLocations = new IntArrayList(nRowsEst);
            IntArrayList endLocations = (hasEndLocations ? new IntArrayList(nRowsEst) : null);
            List<String> probes = new ArrayList(nRowsEst);

            Map<String, FloatArrayList> dataMap = new HashMap();
            for (String h : dataHeaders) {
                dataMap.put(h, new FloatArrayList(nRowsEst));
            }

            // Begin loop through rows
            String chromosome = chrSummary.getName();
            boolean chromosomeStarted = false;
            String nextLine = reader.readLine();

            while ((nextLine != null) && (nextLine.trim().length() > 0)) {

                if (!nextLine.startsWith("#")) {
                    try {
                        String[] tokens = Globals.tabPattern.split(nextLine, -1);

                        String thisChromosome = genome.getChromosomeAlias(tokens[chrColumn].trim());
                        if (thisChromosome.equals(chromosome)) {
                            chromosomeStarted = true;

                            // chromosomeData.setMarkerId(nRows, tokens[0]);

                            // The probe.  A new string is created to prevent holding on to the entire row through a substring reference
                            String probe = new String(tokens[probeColumn]);
                            probes.add(probe);

                            int start = ParsingUtils.parseInt(tokens[startColumn].trim()) - startBase;
                            if (hasEndLocations) {
                                endLocations.add(ParsingUtils.parseInt(tokens[endColumn].trim()));
                            }

                            startLocations.add(start);

                            if(tokens.length <= firstDataColumn + (dataHeaders.length - 1)*skipColumns){
                                String msg = "Line has too few data columns: " + nextLine;
                                log.error(msg);
                                throw new RuntimeException(msg);
                            }

                            for (int idx = 0; idx < dataHeaders.length; idx++) {
                                int i = firstDataColumn + idx * skipColumns;
                                float copyNo = i <= lastDataColumn ? readFloat(tokens[i]) : Float.NaN;
                                String heading = dataHeaders[idx];
                                dataMap.get(heading).add(copyNo);
                            }


                        } else if (chromosomeStarted) {
                            break;
                        }

                    } catch (NumberFormatException numberFormatException) {

                        // Skip line
                        log.info("Skipping line (NumberFormatException) " + nextLine);
                    }
                }

                nextLine = reader.readLine();
            }

            // Loop complete
            ChromosomeData cd = new ChromosomeData(chrSummary.getName());
            cd.setProbes(probes.toArray(new String[]{}));
            cd.setStartLocations(startLocations.toArray());
            if (hasEndLocations) {
                cd.setEndLocations(endLocations.toArray());
            }

            for (String h : dataHeaders) {
                cd.setData(h, dataMap.get(h).toArray());
            }

            return cd;

        } catch (IOException ex) {
            log.error("Error parsing cn file", ex);
            throw new RuntimeException("Error parsing cn file", ex);
        }

    }

    /**
     * Note:  This is an exact copy of the method in ExpressionFileParser.  Refactor to merge these
     * two parsers, or share a common base class.
     *
     * @param comment
     * @param dataset
     */
    private void parseDirective(String comment, IGVDataset dataset) {

        String tmp = comment.substring(1, comment.length());
        if (tmp.startsWith("track")) {
            ParsingUtils.parseTrackLine(tmp, dataset.getTrackProperties());

        } else if (tmp.startsWith("columns")) {
            parseColumnLine(tmp);
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
                    dataset.setTrackType(TrackType.valueOf(tokens[1].trim().toUpperCase()));
                } catch (Exception exception) {

                    // Ignore
                }
            } else if (key.equals("coords")) {

                startBase = Integer.parseInt(tokens[1].trim());

            }
        }
    }

    private boolean isCopyNumberFileExt(String filename) {
        String tmp = ((filename.endsWith(".txt") || filename.endsWith(".tab") || filename.endsWith(".xls")
                ? filename.substring(0, filename.length() - 4) : filename)).toLowerCase();
        return tmp.endsWith(".cn") || tmp.endsWith(".xcn") || tmp.endsWith(".snp");
    }

    private boolean isLOHFileExt(String filename) {
        String tmp = (filename.endsWith(".txt") || filename.endsWith(".tab") || filename.endsWith(".xls")
                ? filename.substring(0, filename.length() - 4) : filename);
        return tmp.endsWith(".loh");
    }

    /**
     * Return the sample headings for the copy number file.
     *
     * @param tokens
     * @param skipColumns
     * @return
     */
    public String[] getHeadings(String[] tokens, int skipColumns) {
        return getHeadings(tokens, skipColumns, false);
    }

    /**
     * Return the sample headings for the copy number file.
     *
     * @param tokens
     * @param skipColumns
     * @param removeDuplicates , whether to remove any duplicate headings
     * @return
     */
    public String[] getHeadings(String[] tokens, int skipColumns, boolean removeDuplicates) {

        ArrayList headings = new ArrayList();
        String previousHeading = null;
        for (int i = firstDataColumn; i <= lastDataColumn; i += skipColumns) {
            if (removeDuplicates) {
                if (previousHeading != null && tokens[i].equals(previousHeading) || tokens[i].equals("")) {
                    continue;
                }

                previousHeading = tokens[i];
            }

            headings.add(tokens[i].trim());
        }

        return (String[]) headings.toArray(new String[0]);
    }


    private void parseColumnLine(String tmp) {
        String[] tokens = tmp.split("\\s+");
        String neg_colnum = "Error parsing column line: " + tmp + "<br>Column numbers must be > 0";
        if (tokens.length > 1) {
            for (int i = 1; i < tokens.length; i++) {
                String[] kv = tokens[i].split("=");
                if (kv.length == 2) {
                    if (kv[0].toLowerCase().equals("chr")) {
                        int c = Integer.parseInt(kv[1]);
                        if (c < 1) {
                            MessageUtils.showMessage(neg_colnum);
                        } else {
                            chrColumn = c - 1;
                        }
                    } else if (kv[0].toLowerCase().equals("start")) {
                        int c = Integer.parseInt(kv[1]);
                        if (c < 1) {
                            MessageUtils.showMessage(neg_colnum);
                        } else {
                            startColumn = c - 1;
                        }


                    } else if (kv[0].toLowerCase().equals("end")) {
                        int c = Integer.parseInt(kv[1]);
                        if (c < 1) {
                            MessageUtils.showMessage(neg_colnum);
                        } else {
                            endColumn = c - 1;
                            hasEndLocations = true;
                        }


                    } else if (kv[0].toLowerCase().equals("probe")) {
                        int c = Integer.parseInt(kv[1]);
                        if (c < 1) {
                            MessageUtils.showMessage(neg_colnum);
                        } else {
                            probeColumn = c - 1;
                        }


                    } else if (kv[0].toLowerCase().equals("data")) {
                        // examples  4,  4-8,  4-
                        String[] se = kv[1].split("-");
                        int c = Integer.parseInt(se[0]);
                        if (c < 1) {
                            MessageUtils.showMessage(neg_colnum);
                        } else {
                            this.firstDataColumn = c - 1;
                        }
                        if (se.length > 1) {

                            c = Integer.parseInt(se[1]);
                            if (c < 1) {
                                MessageUtils.showMessage(neg_colnum);
                            } else {
                                this.lastDataColumn = c - 1;
                            }
                        }
                    }
                }
            }
        }
    }


    private void updateWholeGenome(String currentChromosome, IGVDataset dataset, String[] headings,
                                   IGVDatasetParser.WholeGenomeData wgData) {


        if (!genome.getHomeChromosome().equals(Globals.CHR_ALL)) {
            return;
        }

        // Update whole genome data
        int[] locations = wgData.locations.toArray();
        if (locations.length > 0) {
            Map<String, float[]> tmp = new HashMap(wgData.data.size());
            for (String s : wgData.headings) {
                tmp.put(s, wgData.data.get(s).toArray());
            }


            GenomeSummaryData genomeSummary = dataset.getGenomeSummary();
            if (genomeSummary == null) {
                genomeSummary = new GenomeSummaryData(genome, headings);
                dataset.setGenomeSummary(genomeSummary);
            }
            genomeSummary.addData(currentChromosome, locations, tmp);

        }
    }

    class WholeGenomeData {

        String[] headings;
        IntArrayList locations = new IntArrayList(50000);
        Map<String, FloatArrayList> data = new HashMap();

        WholeGenomeData(String[] headings) {
            this.headings = headings;
            for (String h : headings) {
                data.put(h, new FloatArrayList(50000));
            }
        }

        int size() {
            return locations.size();
        }
    }

}
