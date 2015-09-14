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

package org.broad.igv.util;

import org.broad.igv.Globals;

import java.io.*;
import java.util.*;

/**
 * Created by jrobinso on 8/28/14.
 * <p/>
 * 1-off utilties for manipulating 1000 genomes project files
 */
public class OneKGUtils {


    static final String prefix = "http://1000genomes.s3.amazonaws.com/";
    static final String prefix2 = "http://www.broadinstitute.org/igvdata/1KG/b37/";

    /**
     * Input files contain lines of the form
     * <p/>
     * data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
     */

    static void createResourceXMLs(List<String> inputFiles, String outputFile) throws IOException {

        Map<String, List<String>> exomeResources = new HashMap<String, List<String>>();
        Map<String, List<String>> lowCovResources = new HashMap<String, List<String>>();
        Set<String> pops = new HashSet<String>();


        BufferedReader br = null;

        for (String inputFile : inputFiles) {
            try {
                br = new BufferedReader(new FileReader(inputFile));

                String line;
                while ((line = br.readLine()) != null) {

                    String[] t1 = Globals.forwardSlashPattern.split(line);
                    String[] t2 = t1[3].split("\\.");

                    String type = t2[5];
                    String name = t2[0] + " " + type;
                    String pop = t2[4];

                    pops.add(pop);

                    Map<String, List<String>> resources = type.equals("exome") ? exomeResources : lowCovResources;

                    List<String> entries = resources.get(pop);

                    if (entries == null) {
                        entries = new ArrayList<String>();
                        resources.put(pop, entries);
                    }

                    String fullPath = prefix + line;
                    if (!(FileUtils.resourceExists(fullPath) && FileUtils.resourceExists(fullPath + ".bai"))) {
                        fullPath = prefix2 + line;
                    }

                    if (!(FileUtils.resourceExists(fullPath) && FileUtils.resourceExists(fullPath + ".bai"))) {
                        System.out.println("Resource not found: " + line);
                    }
                    else {
                        entries.add("<Resource name=\"" + name + "\" path=\"" + fullPath + "\"/>");
                    }
                }

            } finally {
                br.close();
            }

        }

        PrintWriter out = null;

        try {

            out = new PrintWriter(outputFile);

            List<String> popList = new ArrayList(pops);
            Collections.sort(popList);

            for (String pop : popList) {

                out.println("<Category name=\"" + pop + "\">");

                out.println("<Category name=\"exome\">");
                List<String> entries = exomeResources.get(pop);
                for (String e : entries) {
                    out.println(e);
                }
                out.println("</Category>");

                out.println("<Category name=\"low coverage\">");
                entries = lowCovResources.get(pop);
                for (String e : entries) {
                    out.println(e);
                }
                out.println("</Category>");


                out.println("</Category>");
            }

        } finally {
            out.close();
        }
    }


}
