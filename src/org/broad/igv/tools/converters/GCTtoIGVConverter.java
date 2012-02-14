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

package org.broad.igv.tools.converters;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;
import org.apache.log4j.Logger;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.data.expression.GeneToLocusHelper;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tools.sort.SortableRecord;
import org.broad.igv.tools.sort.SortableRecordCodec;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.Comparator;
import java.util.List;


/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GCTtoIGVConverter {


    private static Logger log = Logger.getLogger(GCTtoIGVConverter.class);



    /**
     * Parse the file and output in ".igv" format
     *
     * @return
     */
    public static void convert(ResourceLocator resourceLocator, File outputFile, String probeResource,
                        int maxRecords, File tmpDir, Genome genome) throws IOException {

        ExpressionFileParser.FileType type = ExpressionFileParser.determineType(resourceLocator);

        GeneToLocusHelper locusHelper = new GeneToLocusHelper(probeResource);


        BufferedReader reader = null;
        PrintWriter writer = null;

        SortingCollection cltn = getSortingCollection(maxRecords, tmpDir);
        try {
            reader = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(resourceLocator.getPath())));
            writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

            ExpressionFileParser.FormatDescriptor formatDescriptor = ExpressionFileParser.parseHeader (reader, type, null);
            String [] dataHeadings = formatDescriptor.getDataHeaders();


            // Need a better way to determine type!
            String dataType = resourceLocator.getPath().contains("methylation") ? TrackType.DNA_METHYLATION.toString()
                    : TrackType.GENE_EXPRESSION.toString();

            writer.println("#type=" + dataType);
            writer.print("Chr\tStart\tEnd\tProbe");
            for (String s : dataHeadings) {
                writer.print("\t" + s);
            }
            writer.println();

            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {

                // A gct row can map to multiple loci, normally this indicates a problem with the probe
                DataRow row = new DataRow(nextLine, formatDescriptor);
                String probe = row.getProbe();
                List<Locus> loci = locusHelper.getLoci(probe, row.getDescription(), genome.getId());
                if (loci == null || loci.isEmpty()) {
                    System.out.println("No locus found for: " + probe + "  " + row.getDescription());
                } else {
                    for (Locus locus : loci) {
                        String igvLine = locus.getChr() + "\t" + locus.getStart() + "\t" + locus.getEnd() + "\t" + probe +
                                "\t" + row.getData();
                        cltn.add(new SortableRecord(locus.getChr(), locus.getStart(), igvLine));
                    }
                }
            }


            // Ouputput the sorted file
            CloseableIterator<SortableRecord> iter = cltn.iterator();
            while (iter.hasNext()) {
                SortableRecord al = iter.next();
                writer.println(al.getText());

            }
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (writer != null) {
                writer.close();
            }
        }

    }


    static SortingCollection getSortingCollection(int maxRecords, File tmpDir) {


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


        return SortingCollection.newInstance(SortableRecord.class, codec, comp, maxRecords, tmpDir);
    }


    /**
     * Represents a row of data from a GCT or mage-tab file.  Using this class if more effecient than tokeninzing
     * the entire line.  Some GCT files have over a thousand columns and we're only interested in the first 2
     */
    static class DataRow {

        private String probe;
        private String description;
        private String data;

        DataRow(String string, ExpressionFileParser.FormatDescriptor formatDescriptor) {

            String [] tokens = string.split("\t");
            probe = tokens[formatDescriptor.getProbeColumn()];

            int descriptionColumn = formatDescriptor.getDescriptionColumn();
            description = descriptionColumn < 0 ? "" : tokens[descriptionColumn];

            StringBuffer dataBuffer = new StringBuffer();
            final int[] dataColumns = formatDescriptor.getDataColumns();
            dataBuffer.append(tokens[dataColumns[0]]);
            for(int i=1; i<dataColumns.length; i++) {
                dataBuffer.append('\t');
                dataBuffer.append(tokens[dataColumns[i]]);
            }
            data = dataBuffer.toString();

        }

        private String getProbe() {
            return probe;
        }

        public String getDescription() {
            return description;
        }

        public String getData() {
            return data;
        }
    }

}
