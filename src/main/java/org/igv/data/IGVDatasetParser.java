package org.igv.data;


import org.igv.logging.*;
import org.igv.Globals;
import org.igv.exceptions.ParserException;
import org.igv.feature.genome.Genome;
import org.igv.track.DataType;
import org.igv.track.WindowFunction;
import org.igv.ui.IGV;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;
import org.igv.util.collections.FloatArrayList;
import org.igv.util.collections.IntArrayList;
import htsjdk.tribble.readers.Utf8LineReader;

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

    private static Logger log = LogManager.getLogger(IGVDatasetParser.class);
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

        String format = dataResourceLocator.getFormat();
        switch (format) {
            case "igv":
                chrColumn = 0;
                startColumn = 1;
                endColumn = 2;
                probeColumn = 3;
                firstDataColumn = 4;
                hasEndLocations = true;
                hasCalls = false;
                break;
            case "cn":
            case "snp":
            case "xcn":
            case "loh":
                probeColumn = 0;
                chrColumn = 1;
                startColumn = 2;
                endColumn = -1;
                firstDataColumn = 3;
                hasEndLocations = false;
                hasCalls = format.equals("xcn") || format.equals("snp");
                break;
            case ".expr":
                //gene_id	bundle_id	chr	left	right	FPKM	FPKM_conf_lo	FPKM_conf_hi
                probeColumn = 0;
                chrColumn = 2;
                startColumn = 3;
                endColumn = 4;
                startBase = 1;
                firstDataColumn = 5;
                lastDataColumn = 5;
                hasEndLocations = true;
                break;
            default:
                // TODO -- popup dialog and ask user to define columns,  and csv vs tsv?
                throw new ParserException("Unknown file type: ", 0);

        }
    }

    /**
     * Parse the datafile, loading the data for all chromosomes.
     *
     * @param dataset
     * @return the data for each chromosome, in file order
     */
    public List<ChromosomeData> parse(IGVDataset dataset) {

        int estLineCount = ParsingUtils.estimateLineCount(dataResourceLocator.getPath());
        Map<String, Integer> longestFeatureMap = new HashMap();

        float dataMin = 0;
        float dataMax = 0;

        InputStream is = null;
        Utf8LineReader reader = null;
        String nextLine = null;
        String currentChr = null;
        List<ChromosomeData> chrDataList = new ArrayList();
        String[] headings = null;
        ChromosomeBuffer buffer = null;

        int headerRows = 0;
        int count = 0;

        boolean logNormalized;
        try {

            int skipColumns = hasCalls ? 2 : 1;

            // BufferedReader reader = ParsingUtils.openBufferedReader(dataResourceLocator);
            is = ParsingUtils.openInputStreamGZ(dataResourceLocator);
            reader = new Utf8LineReader(is);

            // Infer datatype from extension.  This can be overriden in the
            // comment section
            if (isCopyNumberFileExt(dataResourceLocator.getPath())) {
                dataset.setTrackType(DataType.COPY_NUMBER);
                dataset.getTrackProperties().setWindowingFunction(WindowFunction.mean);
            } else if (isLOHFileExt(dataResourceLocator.getPath())) {
                dataset.setTrackType(DataType.LOH);
                dataset.getTrackProperties().setWindowingFunction(WindowFunction.mean);
            } else {
                dataset.getTrackProperties().setWindowingFunction(WindowFunction.mean);
            }

            // Parse comments and directives, if any
            nextLine = reader.readLine().trim();
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

            // Update
            int updateCount = 5000;
            while ((nextLine = reader.readLine()) != null) {

                if (igv != null && ++count % updateCount == 0) {
                    igv.setStatusBarMessage("Loaded: " + count + " / " + estLineCount + " (est)");
                }
                // Distance since last sample

                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                int nTokens = tokens.length;
                if (nTokens > 0) {
                    String thisChr = genome.getCanonicalChrName(tokens[chrColumn]);
                    if (!thisChr.equals(currentChr)) {
                        // Finish the previous chromosome, unless this is the first one
                        if (currentChr != null) {
                            addChromosome(currentChr, buffer, dataset, chrDataList);
                        }
                        for (ChromosomeData cd : chrDataList) {
                            if (cd.getChr().equals(thisChr)) {
                                throw new ParserException("Rows for " + thisChr + " are not contiguous.  Data in .igv and .cn files " +
                                        "must be grouped by chromosome, and within each chromosome sorted by start position.  " +
                                        "Use igvtools (File > Run igvtools..) to sort the file.", count + headerRows);
                            }
                        }

                        // Start the next chromosome
                        currentChr = thisChr;
                        buffer = new ChromosomeBuffer(headings);
                    }

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
                            int end = ParsingUtils.parseInt(tokens[endColumn].trim());
                            length = end - location + 1;
                            buffer.ends.add(end);

                        } catch (NumberFormatException numberFormatException) {
                            log.error("Column " + tokens[endColumn] + " is not a number");
                            throw new ParserException("Column " + (endColumn + 1) +
                                    " must contain an integer value." + " Found: " + tokens[endColumn],
                                    count + headerRows, nextLine);
                        }
                    }

                    updateLongestFeature(longestFeatureMap, thisChr, length);

                    if (buffer.starts.size() > 0 && buffer.starts.get(buffer.starts.size() - 1) > location) {
                        throw new ParserException("File is not sorted, .igv and .cn files must be sorted by start position." +
                                " Use igvtools (File > Run igvtools..) to sort the file.", count + headerRows);
                    }

                    buffer.starts.add(location);
                    buffer.probes.add(probeColumn >= 0 ? tokens[probeColumn] : null);

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
                        buffer.data.get(heading).add(copyNo);
                    }
                }
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

        // Add the last chromosome
        if (currentChr != null) {
            addChromosome(currentChr, buffer, dataset, chrDataList);
        }

        dataset.setLogNormalized(logNormalized);
        dataset.setDataMin(dataMin);
        dataset.setDataMax(dataMax);

        return chrDataList;
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
                    dataset.setTrackType(DataType.valueOf(tokens[1].trim().toUpperCase()));
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


    /**
     * Add the buffered rows of a chromosome to the list, and to the whole genome summary if the genome has one.
     */
    private void addChromosome(String chr, ChromosomeBuffer buffer, IGVDataset dataset, List<ChromosomeData> chrDataList) {

        ChromosomeData cd = new ChromosomeData(chr);
        cd.setProbes(buffer.probes.toArray(new String[]{}));
        cd.setStartLocations(buffer.starts.toArray());
        if (buffer.ends != null) {
            cd.setEndLocations(buffer.ends.toArray());
        }
        Map<String, float[]> data = new HashMap(buffer.data.size());
        for (String h : buffer.headings) {
            float[] values = buffer.data.get(h).toArray();
            cd.setData(h, values);
            data.put(h, values);
        }
        chrDataList.add(cd);

        if (genome.getHomeChromosome().equals(Globals.CHR_ALL) && cd.getStartLocations().length > 0) {
            GenomeSummaryData genomeSummary = dataset.getGenomeSummary();
            if (genomeSummary == null) {
                genomeSummary = new GenomeSummaryData(genome, buffer.headings);
                dataset.setGenomeSummary(genomeSummary);
            }
            genomeSummary.addData(chr, cd.getStartLocations(), cd.getEndLocations(), data);
        }
    }

    /**
     * The rows of one chromosome while it is being read.
     */
    class ChromosomeBuffer {

        String[] headings;
        IntArrayList starts = new IntArrayList(50000);
        IntArrayList ends;
        List<String> probes = new ArrayList(50000);
        Map<String, FloatArrayList> data = new HashMap();

        ChromosomeBuffer(String[] headings) {
            this.headings = headings;
            this.ends = hasEndLocations ? new IntArrayList(50000) : null;
            for (String h : headings) {
                data.put(h, new FloatArrayList(50000));
            }
        }
    }

}
