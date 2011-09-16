package org.broad.igv.hic.tools;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.Comparator;

/**
 * @author jrobinso
 * @date Aug 17, 2010
 */
public class AlignmentsSorter {


    public static void main(String[] args) throws IOException {
        String ifile = "/Users/jrobinso/projects/hi-c/test/data/selected_formatted.txt";
        String ofile = "/Users/jrobinso/projects/hi-c/test/data/selected_formatted.sorted.txt";
        File tmpdir = new File("/Users/jrobinso/tmp");
        sort(ifile, ofile, tmpdir);

    }

    public static void sort(String inputFile, String outputFile, File tmpdir) throws IOException {

        int maxRecordsInRam = 500000;

        Comparator<AlignmentRecord> comp = new Comparator<AlignmentRecord>() {
            public int compare(AlignmentRecord rec1, AlignmentRecord rec2) {
                int c = rec1.chr1.compareTo(rec2.chr1);
                if (c != 0) {
                    return c;
                } else {
                    return rec1.pos1 - rec2.pos1;
                }

            }
        };


        SortingCollection<AlignmentRecord> cltn =
                SortingCollection.newInstance(AlignmentRecord.class, new Codec(), comp, maxRecordsInRam, tmpdir);


        // Parse input, create records, and add to sorting collection
        BufferedReader br = null;
        try {
            br = ParsingUtils.openBufferedReader(inputFile);

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                try {
                    String[] tokens = nextLine.split("\\s+");
                    AlignmentRecord rec = new AlignmentRecord();
                    rec.line = nextLine;
                    rec.chr1 = tokens[1].trim();
                    rec.pos1 = Integer.parseInt(tokens[2]);
                    cltn.add(rec);
                } catch (Exception e) {
                    System.err.println("Skipping line: " + nextLine + "  " + e.getMessage());
                }

            }
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }
        cltn.doneAdding();

        // Print results of sort to output file
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            CloseableIterator<AlignmentRecord> iter = cltn.iterator();
            while (iter.hasNext()) {
                AlignmentRecord rec = iter.next();
                pw.println(rec.line);

            }
            //for (AlignmentRecord rec : cltn) {
            //    pw.println(rec.line);
            //}
        } finally {
            if (pw != null) pw.close();
        }


    }


    static class AlignmentRecord {
        String chr1;
        int pos1;
        String line;
    }


    static class Codec implements SortingCollection.Codec<AlignmentRecord> {

        PrintWriter os;
        BufferedReader is;

        public void setOutputStream(OutputStream outputStream) {
            os = new PrintWriter(outputStream);
        }

        public void setInputStream(InputStream inputStream) {
            is = new BufferedReader(new InputStreamReader(inputStream), 100);
        }

        public void encode(AlignmentRecord rec) {

            os.println(rec.line);
            os.flush();
        }

        public AlignmentRecord decode() {
            AlignmentRecord rec = new AlignmentRecord();
            try {
                rec.line = is.readLine();
                if (rec.line == null) {
                    return null;
                }
                String[] tokens = rec.line.split("\\s+");
                rec.chr1 = tokens[1].trim();
                rec.pos1 = Integer.parseInt(tokens[2]);
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            return rec;
        }

        public SortingCollection.Codec<AlignmentRecord> clone() {
            Codec clone = new Codec();
            clone.is = is;
            clone.os = os;
            return clone;

        }
    }
}
