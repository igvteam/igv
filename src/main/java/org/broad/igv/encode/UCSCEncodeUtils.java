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

package org.broad.igv.encode;

import org.broad.igv.Globals;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * Date: 10/31/13
 * Time: 12:16 PM
 */
public class UCSCEncodeUtils {

    private static Set<String> labs = new HashSet<>();
    private static Set<String> dataTypes = new HashSet<>();
    private static Set<String> cells = new HashSet<>();
    private static Set<String> antibodies = new HashSet<>();
    private static Set<String> fileTypes = new HashSet<>();
    private static Set<String> allHeaders = new LinkedHashSet<>();

    private static List<String> rnaChipQualifiers = Arrays.asList("CellTotal", "Longnonpolya", "Longpolya",
            "NucleolusTotal", "ChromatinTotal", "ChromatinTotal", "NucleoplasmTotal");

    public static void main(String[] args) throws IOException {

        updateEncodeTableFile(args[0], args[1]);

    }


    static String[] columnHeadings = {"cell", "dataType", "antibody", "view", "replicate", "type", "lab"};

    private static void updateEncodeTableFile(String inputFile, String outputFile) throws IOException {

        List<EncodeFileRecord> records = new ArrayList<>();

        try (BufferedReader reader = ParsingUtils.openBufferedReader(inputFile)) {
            String rootPath = reader.readLine();

            String hub = null;
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("#")) {
                    hub = nextLine.startsWith("#hub=") ? nextLine.substring(5) : hub;
                } else {
                    String dir = nextLine.equals(".") ? rootPath : rootPath + nextLine;
                    String filesDotTxt = dir + "/files.txt";
                    if (HttpUtils.getInstance().resourceAvailable(filesDotTxt)) {
                        try {
                            parseFilesDotTxt(filesDotTxt, records);
                        } catch (IOException e) {
                            // e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                        }
                    }
                }

            }
            fileTypes.forEach(System.out::println);

            outputRecords(outputFile, records, hub);
        }
    }

    private static void outputRecords(String outputFile, List<EncodeFileRecord> records, String hub) throws IOException {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)))) {
            StringBuilder sb = new StringBuilder("path");
            for (String h : columnHeadings) {
                sb.append("\t").append(h);
            }
            if (hub != null) {
                sb.append("\thub");
            }
            pw.println(sb.toString());

            for (EncodeFileRecord rec : records) {
                sb = new StringBuilder(rec.getPath());
                for (String h : columnHeadings) {
                    sb.append("\t").append(Optional.ofNullable(rec.getAttributeValue(h)).orElse(""));
                }
                if (hub != null) {
                    sb.append("\t").append(hub);
                }
                pw.println(sb.toString());
            }
        }
    }

    static HashSet knownFileTypes = new HashSet(Arrays.asList(
            "bam", "bigBed", "bed", "bb", "bw", "bigWig", "gtf", "broadpeak", "narrowpeak", "gappedpeak", "regionpeak", "gff"));

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

            String ftype = df.getFileType() == null ? null : df.getFileType().toLowerCase();
            if (knownFileTypes.contains(ftype)) {
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
