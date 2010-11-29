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

package org.broad.igv.tools.sort;

import jargs.gnu.CmdLineParser;
import org.broad.tribble.readers.AsciiLineReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

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
    private int maxRecords;
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
            }
            catch (NumberFormatException e) {
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
        } else if (shortFN.endsWith(".gff")) {
            return new GFFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".vcf")) {
            return new VCFSorter(inputFile, outputFile);
        } else if (shortFN.endsWith(".psl") || shortFN.endsWith(".pslx")) {
            return new BedSorter(inputFile, outputFile);
        } else {
            System.out.println("Unknown file type or sorting not supported for: " + inputFile.getName());
            return null;
        }
    }

    public Sorter(File inputFile, File outputFile) {
        this.inputFile = inputFile;
        this.outputFile = outputFile;
        this.tmpDir = new File(System.getProperty("java.io.tmpdir"));
    }

    public void run() {

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
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                fis.close();
                writer.close();
            } catch (IOException ex) {
            }
        }
    }

    abstract Parser getParser();

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
