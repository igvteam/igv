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

package org.broad.igv.tools.sort;

import jargs.gnu.CmdLineParser;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.ChromosomeNameComparator;
import org.broad.igv.feature.tribble.MUTCodec;
import org.broad.igv.gwas.GWASParser;
import org.broad.igv.track.GFFFeatureSource;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: nazaire
 * Date: Jun 2, 2009
 */
public abstract class Sorter {

    static private Logger log = Logger.getLogger(Sorter.class);

    static int MAX_RECORDS_IN_RAM = 500000;
    protected File inputFile;

    private File outputFile;
    private boolean writeStdOut = false;
    private int maxRecords = MAX_RECORDS_IN_RAM;

    /**
     * Directory used for storing temporary data files
     */
    private File tmpDir;
    static final String usageString = "igvtools sort <inputFile> [outputFile]";
    protected Comparator<SortableRecord> comparator = getDefaultComparator();

    public static Sorter getSorter(String[] argv) {

        if (argv.length < 2) {
            System.out.println(usageString);
            System.exit(-1);
        }

        CmdLineParser parser = new CmdLineParser();
        CmdLineParser.Option tmpDirOption = parser.addStringOption('t', "tmpDir");
        CmdLineParser.Option maxRecordsOption = parser.addStringOption('m', "maxRecords");

        try {
            parser.parse(argv);
        } catch (CmdLineParser.OptionException e) {
            String msg = "Error parsing command line " + e.getMessage();
            log.error(msg, e);
        }
        String[] nonOptionArgs = parser.getRemainingArgs();

        File inputFile = new File(nonOptionArgs[0]);
        if (!inputFile.exists()) {
            String msg = "Error: " + inputFile.getAbsolutePath() + " does not exist.";
            log.error(msg);
            throw new RuntimeException(msg);
        }

        File outputFile = new File(nonOptionArgs[1]);

        Sorter sorter = getSorter(inputFile, outputFile);

        String tmpDirName = (String) parser.getOptionValue(tmpDirOption);
        if (tmpDirName != null) {
            File tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                String msg = "Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.";
                log.error(msg);
                throw new RuntimeException(msg);
            }
            sorter.setTmpDir(tmpDir);
        }

        int mr = MAX_RECORDS_IN_RAM;
        String maxRecordsString = (String) parser.getOptionValue(maxRecordsOption);
        if (maxRecordsOption != null) {
            try {
                mr = Integer.parseInt(maxRecordsString);
            } catch (NumberFormatException e) {
                log.warn("Warning: max records is not an integer: (" + maxRecordsString + ").  Setting" +
                        "max records to " + MAX_RECORDS_IN_RAM);
                mr = MAX_RECORDS_IN_RAM;
            }
        }

        sorter.setMaxRecords(mr);
        return sorter;
    }

    public static Sorter getSorter(File inputFile, File outputFile) {
        String shortFN = inputFile.getName().toLowerCase();
        if (shortFN.endsWith(".txt")) {
            shortFN = shortFN.substring(0, shortFN.length() - 4);
        }
        if (shortFN.endsWith(".cn") || shortFN.endsWith(".xcn") || shortFN.endsWith(".snp") || shortFN.endsWith(".igv")) {
            return new CNSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".sam")) {
            return new SAMSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".aligned") || shortFN.endsWith(".bed") || shortFN.endsWith(".bedgraph")) {
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
        } else {
            log.error("Unknown file type or sorting not supported for: " + inputFile.getName());
            return null;
        }
    }

    /**
     * @param inputFile
     * @param outputFile If null, we write to stdout
     */
    public Sorter(File inputFile, File outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.writeStdOut = outputFile == null;
        this.tmpDir = new File(System.getProperty("java.io.tmpdir"), System.getProperty("user.name"));

        System.setProperty("snappy.disable", "true");
        if (!tmpDir.exists()) {
            tmpDir.mkdir();
        }
    }

    public void run() throws IOException {

        FileInputStream fis = null;
        PrintWriter writer = null;

        try {
            fis = new FileInputStream(inputFile);
            Writer rawWriter;
            if (writeStdOut) {
                rawWriter = new OutputStreamWriter(System.out);
            } else {
                rawWriter = new FileWriter(this.outputFile);
            }
            writer = new PrintWriter(new BufferedWriter(rawWriter));

            SortableRecordCodec codec = new SortableRecordCodec();

            SortingCollection cltn = SortingCollection.newInstance(SortableRecord.class, codec, comparator, maxRecords, tmpDir);

            Parser parser = getParser();
            AsciiLineReader reader = new AsciiLineReader(fis);

            String firstDataRow = writeHeader(reader, writer);
            if (firstDataRow != null) {
                cltn.add(parser.createRecord(firstDataRow));
            }

            SortableRecord next = null;
            while ((next = parser.readNextRecord(reader)) != null) {
                cltn.add(next);
            }


            CloseableIterator<SortableRecord> iter = cltn.iterator();
            while (iter.hasNext()) {
                SortableRecord al = iter.next();
                writer.println(al.getText());

            }
            iter.close();
        } finally {
            if (fis != null) fis.close();
            if (writer != null) writer.close();
        }
    }

    void setComparator(Comparator<SortableRecord> comparator) {
        this.comparator = comparator;
    }

    public static Comparator<SortableRecord> getDefaultComparator() {
        Comparator<SortableRecord> comp = new Comparator<SortableRecord>() {
            private Comparator<String> nameComparator = ChromosomeNameComparator.get();

            public int compare(SortableRecord o1, SortableRecord o2) {

                int nameComp = nameComparator.compare(o1.getChromosome(), o2.getChromosome());
                if (nameComp == 0) {
                    return o1.getStart() - o2.getStart();
                } else {
                    return nameComp;
                }
            }
        };
        return comp;
    }

    abstract Parser getParser() throws IOException;

    /**
     * Write the header to the output file. Since many readers can't help but read
     * one feature line, that line should be returned and will then be treated as a record
     *
     * @param reader
     * @param writer
     * @return
     * @throws IOException
     */
    abstract String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException;

    public void setTmpDir(File tmpDir) {
        this.tmpDir = tmpDir;
    }

    public void setMaxRecords(int maxRecords) {
        this.maxRecords = maxRecords;
    }

    public void setWriteStdOut(boolean writeStdOut) {
        this.writeStdOut = writeStdOut;
    }
}
