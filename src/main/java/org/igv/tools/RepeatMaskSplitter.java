package org.igv.tools;

import org.igv.Globals;
import org.igv.util.ParsingUtils;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Properties;

/**
 * Splits a repeat mask file downloaded from UCSC into multiple files,  one per repeat class.
 * Assumes downloaded columns as follows (use table browser, and "select columns" option
 * <p/>
 * genoName  genoStart  genoEnd  strand  repName repClass repFamily
 * <p/>
 * Assumes file is sorted by chromosome
 *
 * @author jrobinso
 */
public class RepeatMaskSplitter {

    public static void main(String[] args) {
        File file = new File(args[0]);
        split(file);
    }

    public static void split(File inputFile) {

        int binCol = 0;
        int millDivCol = 2;
        int millDelCol = 3;
        int millInsCol = 4;
        int chrCol = 5;
        int startCol = 6;
        int endCol = 7;
        int strandCol = 9;
        int nameCol = 10;
        int classCol = 11;
        int famCol = 12;

        Map<String, LinkedHashMap<String, String>> fileMappings = new HashMap();

        AsciiLineReader reader = null;
        HashMap<String, PrintWriter> writers = new HashMap();
        PrintWriter allWriter = null;
        try {
            reader = new AsciiLineReader(new FileInputStream(inputFile));
            // Skip header
            reader.readLine();
            String nextLine;
            File dir = inputFile.getParentFile();

            allWriter = new PrintWriter(new BufferedWriter(new FileWriter("rmsk.bed")));
            allWriter.println("track name=\"Repeat Masker\" \" gffTags=\"on\"");

            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                String chr = tokens[chrCol];

                String repClass = tokens[classCol];
                if (repClass.contains("?")) {
                    continue;
                }
                String filename = repClass + ".bed";

                // Get or create file writer for the class
                PrintWriter pw = writers.get(filename);
                if (pw == null) {
                    File outputFile = new File(dir, filename);
                    pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
                    pw.println("track name=\"" + repClass + "\" gffTags=\"on\"");
                    writers.put(filename, pw);
                }

                String nm = tokens[nameCol];
                String fam = tokens[famCol];

                String name = "Name=" + nm + ";Class=" + repClass + ";Family=" + fam;

                pw.print(chr);
                pw.print("\t");
                pw.print(Integer.parseInt(tokens[startCol]));
                pw.print("\t");
                pw.print(Integer.parseInt(tokens[endCol]));
                pw.print("\t");
                pw.print(name);
                pw.print("\t");
                pw.print(tokens[strandCol]);
                pw.println();

                allWriter.print(chr);
                allWriter.print("\t");
                allWriter.print(Integer.parseInt(tokens[startCol]));
                allWriter.print("\t");
                allWriter.print(Integer.parseInt(tokens[endCol]));
                allWriter.print("\t");
                allWriter.print(name);
                allWriter.print("\t");
                allWriter.print(tokens[strandCol]);
                allWriter.println();

            }

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            reader.close();
           // allWriter.close();
            closeWriters(writers);
        }

    }

    private static void closeWriters(HashMap<String, PrintWriter> writers) {
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
        writers.clear();
    }
}

