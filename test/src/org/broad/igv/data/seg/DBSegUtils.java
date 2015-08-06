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

package org.broad.igv.data.seg;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.junit.Ignore;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date 10/15/11
 */
@Ignore
public class DBSegUtils {

    public static void main(String[] args) throws IOException {
        // String seg = "http://www.broadinstitute.org/igvdata/tcga/gbmsubtypes/Broad.080528.subtypes.seg.gz";
        // String oFile = "gbm_subtypes.sql";
        // segToDB(seg, oFile);
        generateSampleInfoInserts("http://www.broadinstitute.org/igvdata/tcga/gbmsubtypes/sampleTable.txt");
    }

    /**
     * Convert a seg file to a DB table insert (Maggie's table format)
     * <p/>
     * CREATE TABLE `CNV` (
     * `Chromosome` varchar(255) DEFAULT NULL,
     * `Start` varchar(255) DEFAULT NULL,
     * `Stop` varchar(255) DEFAULT NULL,
     * `Event` varchar(255) DEFAULT NULL,
     * `Length` varchar(255) DEFAULT NULL,
     * `% of CNV Overlap` varchar(255) DEFAULT NULL,
     * `Probe Median` varchar(255) DEFAULT NULL,
     * `% Heterozygous` varchar(255) DEFAULT NULL,
     * `Sample` varchar(255) DEFAULT NULL
     * ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
     * <p/>
     * -- ----------------------------
     * -- Records of CNV
     * -- ----------------------------
     * INSERT INTO `CNV` VALUES ('chr6', '166,262,887', '170,899,992', 'Allelic Imbalance', '4637106', '52.33536872682417', '0.48076730966567993', '87.35632183908046', '44-65-05-2_T_Sty');
     *
     * @param segFile
     * @param outputFile
     */
    static void segToDB(String segFile, String outputFile) throws IOException {

        ResourceLocator loc = new ResourceLocator(segFile);
        SegmentedAsciiDataSet ds = new SegmentedAsciiDataSet(loc, null);
        SegmentFileParser parser = new SegmentFileParser(loc);
        parser.loadSegments(loc, null);

        Map<String, String> sampleMap = loadSampleMappings("http://www.broadinstitute.org/igvdata/tcga/gbmsubtypes/cn.sampleMappings.txt");
        PrintWriter pw = new PrintWriter(new FileWriter(outputFile));

        Set<String> chrs = ds.getChromosomes();
        List<String> samples = ds.getSampleNames();
        for (String chr : chrs) {
            for (String id : samples) {
                String sample = sampleMap.get(id);
                List<LocusScore> segments = ds.getSegments(id, chr);
                if (segments != null && sample != null) {
                    for (LocusScore ls : segments) {


                        Segment seg = (Segment) ls;
                        pw.print("INSERT INTO `CNV` VALUES ('");
                        pw.print(chr + "', '");
                        pw.print(seg.getStart() + "', '");
                        pw.print(seg.getEnd() + "', '");
                        pw.print("CN', '");
                        int length = seg.getEnd() - seg.getStart();
                        pw.print(length + "', '");
                        pw.print("0.0', '");
                        pw.print(seg.getScore() + "', '");
                        pw.print("0', '");
                        pw.println(sample + "');");

                    }
                }
            }
        }
        pw.println("commit;");
        pw.close();

    }

    static Map<String, String> loadSampleMappings(String file) throws IOException {

        BufferedReader reader = ParsingUtils.openBufferedReader(file);
        Map<String, String> map = new HashMap<String, String>();
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            String[] tokens = nextLine.split("\\s+");
            map.put(tokens[0], tokens[1]);

        }

        reader.close();
        return map;

    }


    static void generateSampleInfoInserts(String file) throws IOException {

        BufferedReader reader = ParsingUtils.openBufferedReader(file);
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            if (tokens.length > 6) {
                System.out.print("INSERT INTO `SAMPLE_INFO` VALUES (");
                for (int i = 0; i < tokens.length; i++) {
                    String val = tokens[i].equals("NA") ? "" : tokens[i];
                    System.out.print("'" + val + "'");
                    if (i < tokens.length - 1) {
                        System.out.print(", ");
                    }
                }
                System.out.println(");");
            }

        }
        reader.close();
    }
}
