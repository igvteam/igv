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

package org.broad.igv.batch;


import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;

public class TestClient {

    public static void main(String args[]) throws IOException {
        Socket socket = null;
        PrintWriter out = null;
        BufferedReader in = null;
        try {
            socket = new Socket("127.0.0.1", 60151);
            out = new PrintWriter(socket.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
            //testGOTO(out, in);
            //runBatchFile(out, in, "test/data/batch/test_commands.txt");
            //runBatchFile(out, in, "test/data/batch/load_bigwig.txt");
            sendCommand(out, in, "load https://krishna.gs.washington.edu/download/CADD/bigWig/v1.6/GRCh37/CADDv1.6_GRCh37_whole_genome_SNVs.bw");
           // getToolsYaml(out, in, null);
           // manyLoci(out, in);
           // snapshot(out, in);

        } catch (UnknownHostException e) {
            System.err.println("Unknown host exception: " + e.getMessage());
            System.exit(1);
        } catch (IOException e) {
            e.printStackTrace();
            System.err.println("Couldn't get I/O for " + "the connection to IGV");
            System.exit(1);
        } finally {
            in.close();
            out.close();
            socket.close();
        }
    }

    private static void getToolsYaml(PrintWriter out, BufferedReader in, String inputFile) throws IOException {

        out.println("toolsYaml");
        String response = in.readLine();
        System.out.println(response);

    }

    private static void sendCommand(PrintWriter out, BufferedReader in, String command) throws IOException {
        System.out.println("Executing Command: " + command);
        out.println(command);
        String response = in.readLine();
        System.out.println("Response: " + response);
    }

    private static void runBatchFile(PrintWriter out, BufferedReader in, String inputFile) throws IOException {

        String inLine;

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(inputFile);

            while ((inLine = reader.readLine()) != null) {
                if (!(inLine.startsWith("#") || inLine.startsWith("//"))) {
                    System.out.println("Executing Command: " + inLine);
                    out.println(inLine);
                    String response = in.readLine();
                    System.out.println("Response: " + response);
                }
            }


        } catch (IOException ioe) {
            throw new DataLoadException(ioe.getMessage(), inputFile);
        } finally {

            if (reader != null) reader.close();

        }

    }

    private static void testMultiLocus(PrintWriter out, BufferedReader in) throws IOException {

        String cmd = "load https://s3.amazonaws.com/igv.org.demo/GBM-TP.seg.gz";
        out.println(cmd);
        String response = in.readLine();
        System.out.println(cmd + " " + response);

        cmd = "goto EGFR PTEN";
        out.println(cmd);
        response = in.readLine();
        System.out.println(cmd + " " + response);


    }

    private static void testGOTO(PrintWriter out, BufferedReader in) throws IOException {

        out.println("new");
        String response = in.readLine();
        System.out.println(response);

        out.println("load https://1000genomes.s3.amazonaws.com/phase3/data/HG01883/alignment/HG01883.mapped.ILLUMINA.bwa.ACB.low_coverage.20130415.bam");
        response = in.readLine();
        System.out.println(response);

        out.println("load https://1000genomes.s3.amazonaws.com/phase3/data/HG01879/alignment/HG01879.mapped.ILLUMINA.bwa.ACB.low_coverage.20120522.bam");
        response = in.readLine();
        System.out.println(response);

        int cnt = 10;
        while (cnt-- >= 0) {
            int pos = 1000000 + (int) (Math.random() * 1000000);
            String cmd = "goto chr1:" + pos;
            out.println(cmd);
            response = in.readLine();
            System.out.println(cmd + " " + response);

            out.println("snapshot");
            response = in.readLine();
            System.out.println("snapshot " + response);

        }
    }

    private static void manyLoci(PrintWriter out, BufferedReader in) throws IOException {

        int cnt = 100;
        while (cnt-- >= 0) {
            int pos = 1000000 + (int) (Math.random() * 1000000);
            String cmd = "goto chr1:" + pos;
            out.println(cmd);
            String response = in.readLine();
            System.out.println(cmd + " " + response);

            out.println("snapshot");
            response = in.readLine();
            System.out.println("snapshot " + response);

        }
    }

    private static void snapshot(PrintWriter out, BufferedReader in) throws IOException {

        out.println("snapshot");
        String response = in.readLine();
        System.out.println(response);
    }

}
