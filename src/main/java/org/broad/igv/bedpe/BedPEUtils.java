package org.broad.igv.bedpe;

import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.*;

/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPEUtils {

    public static void main(String[] args) throws IOException {

        juiceboxToBedPE(args[0], args[1]);

    }

    public static void interactionToBedPE(String ifile, String ofile) throws IOException {

        BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = new BufferedReader(new FileReader(ifile));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(ofile)));

            String nextLine;
            while ((nextLine = br.readLine()) != null) {

                String[] tokens = ParsingUtils.WHITESPACE_PATTERN.split(nextLine);

                if (tokens.length == 3) {

                    String[] t1 = tokens[0].split(":");
                    String chr1 = t1[0];

                    t1 = t1[1].split("-");
                    String start1 = t1[0];
                    String end1 = t1[1];

                    t1 = tokens[1].split(":");
                    String chr2 = t1[0];

                    t1 = t1[1].split("-");
                    String start2 = t1[0];
                    String end2 = t1[1];

                    String score = tokens[2];
                    String name = tokens[0] + "->" + tokens[1];

                    pw.println(chr1 + "\t" + start1 + "\t" + end1 + "\t" + chr2 + "\t" + start2 + "\t" + end2 + "\t" + name + "\t" + score);

                } else {
                    System.out.println("Skipping line: " + nextLine);
                }
            }
        } finally {
            if (pw != null) pw.close();
            if (br != null) br.close();
            ;
        }
    }

    public static void juiceboxToBedPE(String ifile, String ofile) throws IOException {

        BufferedReader br = null;
        PrintWriter pw = null;

        try {
            br = new BufferedReader(new FileReader(ifile));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(ofile)));

            String nextLine;
            br.readLine(); // Skipping first line

            pw.println("#columns color=11;thickness=12");

            while ((nextLine = br.readLine()) != null) {

                if (nextLine.startsWith("#")) {
                    pw.println(nextLine);
                } else {
                    String[] tokens = ParsingUtils.WHITESPACE_PATTERN.split(nextLine);

                    if (tokens.length >= 6) {

                        for (int i = 0; i < 6; i++) {
                            pw.print(tokens[i] + "\t");
                        }
                        for (int i = 0; i < 4; i++) {
                            pw.print(".\t");
                        }
                        pw.println(tokens[6] + "\t" + "2\t");


                    } else {
                        System.out.println("Skipping line: " + nextLine);
                    }
                }
            }
        } finally {
            if (pw != null) pw.close();
            if (br != null) br.close();
            ;
        }
    }

    public static void packFeatures(List<BedPE> features, int maxRows) {

        int gap = 2;

        int start = features.get(0).getStart();
        int end = features.get(features.size() - 1).getStart();

        int bucketStart = start;
        int nextStart = bucketStart;

        Map<Integer, Stack<BedPE>> bucketList = new HashMap<>();
        for(BedPE alignment : features) {

            int buckListIndex = Math.max(0, alignment.getStart() - bucketStart);
            if (bucketList.get(buckListIndex) ==  null) {
                bucketList.put(buckListIndex, new Stack<>());
            }
            bucketList.get(buckListIndex).add(alignment);
        }


        int row = 0;
        int allocatedCount = 0;
        int lastAllocatedCount = -1;

        while (allocatedCount < features.size() && row <= maxRows) {

            while (nextStart <= end) {

                Stack<BedPE> bucket = null;

                int index = 0;
                while (bucket == null && nextStart <= end) {
                    index = nextStart - bucketStart;
                    if (bucketList.get(index) == null) {
                        ++nextStart;                     // No buckets at this index
                    } else {
                        bucket = bucketList.get(index);
                    }
                } // while (bucket)

                if (bucket == null) {
                    break;
                }

                BedPE feature = bucket.pop();
                if (0 ==  bucket.size()) {
                    bucketList.remove(index);
                }

                feature.setRow(row);

                nextStart = feature.getEnd() + gap;
                ++allocatedCount;

            } // while (nextStart)

            row++;
            nextStart = bucketStart;

            if (allocatedCount == lastAllocatedCount) break;   // Protect from infinite loops

            lastAllocatedCount = allocatedCount;

        } // while (allocatedCount)

    }


}
