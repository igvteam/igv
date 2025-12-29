package org.igv.tools.converters;

import htsjdk.samtools.util.CloseableIterator;
import org.igv.logging.*;
import org.igv.data.expression.ExpressionFileParser;
import org.igv.data.expression.GeneToLocusHelper;
import org.igv.feature.Locus;
import org.igv.feature.genome.Genome;
import org.igv.tools.sort.SortableRecord;
import org.igv.tools.sort.SortableRecordCodec;
import org.igv.track.TrackType;
import org.igv.util.ParsingUtils;
import org.igv.util.ResourceLocator;
import org.igv.util.collections.SortingCollection;

import java.io.*;
import java.util.Comparator;
import java.util.List;


/**
 * @author jrobinso
 * @date Oct 9, 2010
 */
public class GCTtoIGVConverter {


    private static Logger log = LogManager.getLogger(GCTtoIGVConverter.class);



    /**
     * Parse the file and output in ".igv" format
     *
     * @return
     */
    public static void convert(ResourceLocator resourceLocator, File outputFile, String probeResource,
                        int maxRecords, File tmpDir, Genome genome) throws IOException {

        ExpressionFileParser.FileType type = ExpressionFileParser.determineType(resourceLocator);

        GeneToLocusHelper locusHelper = new GeneToLocusHelper(probeResource, genome);


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
                    log.warn("No locus found for: " + probe + "  " + row.getDescription());
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
