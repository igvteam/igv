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

package org.broad.igv.tools.sort;

import jargs.gnu.CmdLineParser;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;
import org.broad.igv.gwas.GWASParser;
import org.broad.igv.track.GFFFeatureSource;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: nazaire
 * Date: Jun 2, 2009
 */
public abstract class Sorter {

    static int MAX_RECORDS_IN_RAM = 500000;
    File inputFile;
    File outputFile;
    private int maxRecords = MAX_RECORDS_IN_RAM;
    private File tmpDir;
    static final String usageString = "igvtools sort <inputFile> [outputFile]";

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
            System.err.println("Error parsing command line " + e.getMessage());
        }
        String[] nonOptionArgs = parser.getRemainingArgs();

        File inputFile = new File(nonOptionArgs[0]);
        if (!inputFile.exists()) {
            System.out.println("Error: " + inputFile.getAbsolutePath() + " does not exist.  Exiting");
            System.exit(-1);
        }

        File outputFile = new File(nonOptionArgs[1]);

        Sorter sorter = getSorter(inputFile, outputFile);

        String tmpDirName = (String) parser.getOptionValue(tmpDirOption);
        if (tmpDirName != null) {
            File tmpDir = new File(tmpDirName);
            if (!tmpDir.exists()) {
                System.err.println("Error: tmp directory: " + tmpDir.getAbsolutePath() + " does not exist.");
                System.exit(-1);
            }
            sorter.setTmpDir(tmpDir);
        }

        int mr = MAX_RECORDS_IN_RAM;
        String maxRecordsString = (String) parser.getOptionValue(maxRecordsOption);
        if (maxRecordsOption != null) {
            try {
                mr = Integer.parseInt(maxRecordsString);
            } catch (NumberFormatException e) {
                System.out.println("Warning: max records is not an integer: (" + maxRecordsString + ").  Setting" +
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
        } else if (shortFN.endsWith(".aligned") || shortFN.endsWith(".bed")) {
            return new BedSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".sorted")) {
            return new SortedTxtSorter(inputFile, outputFile);
        } else if (GFFFeatureSource.isGFF(shortFN)) {
            return new GFFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".vcf")) {
            return new VCFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".psl") || shortFN.endsWith(".pslx")) {
            return new BedSorter(inputFile, outputFile);
        } else if (GWASParser.isGWASFile(shortFN)) {
            return new GWASSorter(inputFile, outputFile);
        } else {
            System.out.println("Unknown file type or sorting not supported for: " + inputFile.getName());
            return null;
        }
    }

    public Sorter(File inputFile, File outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.tmpDir = new File(System.getProperty("java.io.tmpdir"), System.getProperty("user.name"));
        // Disable "Snappy"
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
            writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            SortableRecordCodec codec = new SortableRecordCodec();

            Comparator<SortableRecord> comp = new Comparator<SortableRecord>() {

                public int compare(SortableRecord o1, SortableRecord o2) {


                    String chr1 = o1.getChromosome().replaceFirst("chr", "");
                    String chr2 = o2.getChromosome().replaceFirst("chr", "");
                    int s1 = Integer.MAX_VALUE;
                    try {
                        s1 = Integer.parseInt(chr1);
                    } catch (Exception e) {
                        // ignore
                    }
                    int s2 = Integer.MAX_VALUE;
                    try {
                        s2 = Integer.parseInt(chr2);
                    } catch (Exception e) {
                        // ignre
                    }


                    int t1 = s1 - s2;
                    if (t1 == 0) {
                        chr1 = chr1.replace("M", "Z");
                        chr2 = chr2.replace("M", "Z");
                        t1 = chr1.compareTo(chr2);
                    }
                    if (t1 == 0) {
                        return (int) (o1.getStart() - o2.getStart());
                    } else {
                        return t1;
                    }
                }
            };


            SortingCollection cltn = SortingCollection.newInstance(SortableRecord.class, codec, comp, maxRecords, tmpDir);

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
            try {
                fis.close();
                writer.close();
            } catch (IOException ex) {
            }
        }
    }

    abstract Parser getParser() throws IOException;

    abstract String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException;

    /**
     * @param tmpDir the tmpDir to set
     */
    public void setTmpDir(File tmpDir) {
        this.tmpDir = tmpDir;
    }

    public void setMaxRecords(int maxRecords) {
        this.maxRecords = maxRecords;
    }
}
