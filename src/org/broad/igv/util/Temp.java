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

package org.broad.igv.util;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jan 21, 2010
 * Time: 10:43:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class Temp {


    public static void main(String[] args) throws IOException, InterruptedException {
        fixDicty();
    }

    public static void fixDicty() throws IOException, InterruptedException {

        String dir = "/Users/jrobinso/IGV/dicty";
        Map<String, String> chrMap = new HashMap();
        chrMap.put("DDB0169550", "M");
        chrMap.put("DDB0215018", "3F");
        chrMap.put("DDB0215151", "2F");
        chrMap.put("DDB0220052", "BF");
        chrMap.put("DDB0232428", "1");
        chrMap.put("DDB0232429", "2");
        chrMap.put("DDB0232430", "3");
        chrMap.put("DDB0232431", "4");
        chrMap.put("DDB0232432", "5");
        chrMap.put("DDB0232433", "6");
        chrMap.put("DDB0237465", "R");


        /*Iterator<Map.Entry<String, String>> iter = chrMap.entrySet().iterator();
        Map.Entry<String, String> first = iter.next();
        sedExpression.append("sed 's/" + first.getKey() + "/" + first.getValue() + "/g'  dicty.fa");
        while(iter.hasNext()) {
            first = iter.next();
            sedExpression.append(" | ");
            sedExpression.append("sed 's/" + first.getKey() + "/" + first.getValue() + "/g'");
        }

       sedExpression.append(" > dicty_chr.fa");
       System.out.println(sedExpression.toString());
        */

        File gffDir = new File(dir, "gff3");

        for (File f : gffDir.listFiles()) {
            if (f.getName().endsWith("gff")) {
                StringBuffer sedExpression = new StringBuffer();
                Iterator<Map.Entry<String, String>> iter = chrMap.entrySet().iterator();
                Map.Entry<String, String> first = iter.next();
                sedExpression.append("sed 's/" + first.getKey() + "/" + first.getValue() + "/g' " + f.getAbsolutePath());
                while (iter.hasNext()) {
                    first = iter.next();
                    sedExpression.append(" | ");
                    sedExpression.append("sed 's/" + first.getKey() + "/" + first.getValue() + "/g'");
                }
                String outputFile = f.getAbsolutePath().replace("chromosome_", "chr");
                String cmdLine = sedExpression.toString() + " > " + outputFile;
                Process p = Runtime.getRuntime().exec(cmdLine);
                p.waitFor();
                System.out.println(cmdLine);

                System.out.println(p.exitValue());
            }
        }

        /*
        Process p = Runtime.getRuntime().exec("cat");

        OutputStream os = p.getOutputStream();

        PrintWriter pw = new PrintWriter(new OutputStreamWriter(os));
        pw.println("abcd") ;
        pw.flush();
        os.close();

        InputStream is = p.getInputStream();
        */

    }


    public static void makeSampleInfo() throws Exception {

        Set<String> sec = new HashSet(Arrays.asList("TCGA-02-0010", "TCGA-02-0028", "TCGA-02-0102",
                "TCGA-02-0114", "TCGA-08-0525"));

        Set<String> treat = new HashSet(Arrays.asList("TCGA-02-0001", "TCGA-02-0007",
                "TCGA-02-0010",
                "TCGA-02-0014",
                "TCGA-02-0021",
                "TCGA-02-0024",
                "TCGA-02-0028",
                "TCGA-02-0043",
                "TCGA-02-0057",
                "TCGA-02-0058",
                "TCGA-02-0080",
                "TCGA-02-0083",
                "TCGA-02-0089",
                "TCGA-02-0099",
                "TCGA-02-0102",
                "TCGA-02-0107",
                "TCGA-02-0113",
                "TCGA-02-0114",
                "TCGA-02-0116",
                "TCGA-08-0517",
                "TCGA-08-0525"));

        Set<String> hyper = new HashSet(Arrays.asList("TCGA-02-0010",
                "TCGA-02-0014",
                "TCGA-02-0028",
                "TCGA-02-0043",
                "TCGA-02-0083",
                "TCGA-02-0099",
                "TCGA-02-0114",
                "TCGA-02-0015",
                "TCGA-02-0070"));


        BufferedReader br = new BufferedReader(new FileReader("Public_IGV_TCGA_sample_info_20080827.txt"));
        PrintWriter pw = new PrintWriter(new FileWriter("Sample_info.txt"));

        String nextLine = br.readLine();
        pw.println(nextLine + "\tTreated\tPrimary/Secondary\tHypermutated");

        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            String sample = tokens[0].length() > 12 ? tokens[0].substring(0, 12) : tokens[0];
            String t = treat.contains(sample) ? "Y" : "";
            String s = sec.contains(sample) ? "Secondary" : "Primary";
            String h = hyper.contains(sample) ? "Y" : "";
            pw.println(nextLine + "\t" + t + "\t" + s + "\t" + h);
        }

        pw.close();
        br.close();


    }
}
