package org.broad.igv.util.encode;

import org.broad.igv.Globals;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;

import java.awt.*;
import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 10/31/13
 *         Time: 12:16 PM
 */
public class UCSCEncodeUtils {

    static HashSet<String> labs = new HashSet<String>();
    static HashSet<String> dataTypes = new HashSet<String>();
    static HashSet<String> cells = new HashSet<String>();
    static HashSet<String> antibodies = new HashSet<String>();
    static HashSet<String> fileTypes = new HashSet<String>();
    static HashSet<String> allHeaders = new LinkedHashSet<String>();

    private static List<String> rnaChipQualifiers = Arrays.asList("CellTotal", "Longnonpolya", "Longpolya",
            "NucleolusTotal", "ChromatinTotal", "ChromatinTotal", "NucleoplasmTotal");

    public static void main(String[] args) throws IOException {


//        List<EncodeFileRecord> records = new ArrayList();
//        parseFilesDotTxt(args[0], records);
//        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(args[1])));
//
//        pw.print("path");
//        for (String h : EncodeTableModel.columnHeadings) {
//            pw.print("\t");
//            pw.print(h);
//        }
//        pw.println();
//
//        for (EncodeFileRecord rec : records) {
//            pw.print(rec.getPath());
//            for (String h : EncodeTableModel.columnHeadings) {
//                pw.print("\t");
//                String value = rec.getAttributeValue(h);
//                pw.print(value == null ? "" : value);
//            }
//            pw.println();
//        }
//        pw.close();

        updateEncodeTableFile(args[0], args[1]);

    }

    private static List<EncodeFileRecord> parseTableFile(String url) throws IOException {

        List<EncodeFileRecord> records = new ArrayList<EncodeFileRecord>(20000);

        BufferedReader reader = null;

        reader = ParsingUtils.openBufferedReader(url);

        String[] headers = Globals.tabPattern.split(reader.readLine());

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {
            if (!nextLine.startsWith("#")) {
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                String path = tokens[0];
                Map<String, String> attributes = new HashMap<String, String>();
                for (int i = 0; i < headers.length; i++) {
                    String value = tokens[i];
                    if (value.length() > 0) {
                        attributes.put(headers[i], value);
                    }
                }
                records.add(new EncodeFileRecord(path, attributes));
            }

        }
        return records;
    }


    private static void updateEncodeTableFile(String inputFile, String outputFile) throws IOException {

        List<EncodeFileRecord> records = parseUCSCMasterFile(inputFile);

        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));

        pw.print("path");
        for (String h : EncodeTableModel.columnHeadings) {
            pw.print("\t");
            pw.print(h);
        }
        pw.println();

        for (EncodeFileRecord rec : records) {
            pw.print(rec.getPath());
            for (String h : EncodeTableModel.columnHeadings) {
                pw.print("\t");
                String value = rec.getAttributeValue(h);
                pw.print(value == null ? "" : value);
            }
            pw.println();
        }
        pw.close();
    }

    private static List<EncodeFileRecord> parseUCSCMasterFile(String url) throws IOException {

        List<EncodeFileRecord> records = new ArrayList<EncodeFileRecord>();

        BufferedReader reader = null;
        reader = ParsingUtils.openBufferedReader(url);

        String rootPath = reader.readLine();

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            if (!nextLine.startsWith("#")) {
                String dir = rootPath + nextLine;
                String filesDotTxt = dir + "/files.txt";

                try {
                    if (HttpUtils.getInstance().resourceAvailable(new URL(filesDotTxt))) {
                        parseFilesDotTxt(filesDotTxt, records);
                    }
                } catch (IOException e) {
                    // e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        }

//
//
//        System.out.println("Labs");
//        for (String dt : labs) System.out.println(dt);
//
//        System.out.println();
//        System.out.println("Data Types");
//        for (String dt : dataTypes) System.out.println(dt);
//
//        System.out.println();
//        System.out.println("Antibodies");
//        for (String dt : antibodies) System.out.println(dt);
//
//        System.out.println();
//        System.out.println("cells");
//        for (String dt : cells) System.out.println(dt);

        for (String dt : fileTypes) System.out.println(dt);

        return records;

    }

    static HashSet knownFileTypes = new HashSet(Arrays.asList("bam", "bigBed", "bed", "bb", "bw", "bigWig", "gtf", "broadPeak", "narrowPeak", "gff"));

    public static void parseFilesDotTxt(String url, List<EncodeFileRecord> fileRecords) throws IOException {


        BufferedReader reader = null;


        reader = ParsingUtils.openBufferedReader(url);
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            String[] tokens = Globals.tabPattern.split(nextLine);
            if (tokens.length < 2) continue;

            String fn = tokens[0];

            String[] attributes = Globals.semicolonPattern.split(tokens[1]);

            LinkedHashMap<String, String> kvalues = new LinkedHashMap<String, String>();
            for (String tk : attributes) {

                String[] kv = Globals.equalPattern.split(tk);
                if (kv.length > 1) {
                    kvalues.put(kv[0].trim(), kv[1].trim());
                    allHeaders.add(kv[0].trim());
                }

            }

            // Hack for RnaChip -- need this to disambiguate them
            if ("RnaChip".equals(kvalues.get("dataType"))) {
                for (String qual : rnaChipQualifiers) {
                    if (fn.contains(qual)) {
                        kvalues.put("antibody", qual);
                    }
                }
            }

            EncodeFileRecord df = new EncodeFileRecord(url.replace("files.txt", fn), kvalues);

            if (knownFileTypes.contains(df.getFileType())) {
                fileRecords.add(df);
            }

            dataTypes.add(df.getAttributeValue("dataType"));
            antibodies.add(df.getAttributeValue("antibody"));
            cells.add(df.getAttributeValue("cell"));
            labs.add(df.getAttributeValue("lab"));
            fileTypes.add(df.getFileType());

        }

        reader.close();

    }

    /*
File types
bam
bigBed
shortFrags
csqual
spikeins
bai
pdf
bed
matrix
bigWig
tab
bed9
bedCluster
peptideMapping
csfasta
gtf
fastq
broadPeak
narrowPeak
gff
bedRrbs
bedRnaElements
tgz
bedLogR
peaks
*/


}
