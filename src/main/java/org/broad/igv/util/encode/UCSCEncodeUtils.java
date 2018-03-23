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

    static String[] columnHeadings = {"cell", "dataType", "antibody", "view", "replicate", "type", "lab"};

    private static void updateEncodeTableFile(String inputFile, String outputFile) throws IOException {

        List<EncodeFileRecord> records = new ArrayList<EncodeFileRecord>();

        BufferedReader reader = null;
        reader = ParsingUtils.openBufferedReader(inputFile);

        String rootPath = reader.readLine();

        String hub = null;
        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            if (nextLine.startsWith("#")) {
                if(nextLine.startsWith("#hub=")) {
                    hub = nextLine.substring(5);
                }
            }
            else {
                String dir = nextLine.equals(".") ? rootPath : rootPath + nextLine;
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
        for (String dt : fileTypes) System.out.println(dt);


        outputRecords(outputFile, records, hub);
    }

    private static void outputRecords(String outputFile, List<EncodeFileRecord> records, String hub) throws IOException {
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
        pw.print("path");
        for (String h : columnHeadings) {
            pw.print("\t");
            pw.print(h);
        }
        if(hub != null) {
            pw.print("\thub");
        }
        pw.println();

        for (EncodeFileRecord rec : records) {
            pw.print(rec.getPath());
            for (String h : columnHeadings) {
                pw.print("\t");
                String value = rec.getAttributeValue(h);
                pw.print(value == null ? "" : value);
            }
            if(hub != null) {
                pw.print("\t" + hub);
            }
            pw.println();
        }
        pw.close();
    }

    static HashSet knownFileTypes = new HashSet(Arrays.asList(
            "bam", "bigBed", "bed", "bb", "bw", "bigWig", "gtf", "broadPeak", "narrowPeak", "gappedPeak", "gff"));

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

            String path = fn.startsWith("http") ? fn : url.replace("files.txt", fn);

            EncodeFileRecord df = new EncodeFileRecord(path, kvalues);

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


}
