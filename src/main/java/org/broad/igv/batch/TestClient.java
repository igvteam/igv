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
import org.broad.igv.ui.WaitCursorManager;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.List;

public class TestClient {

    static private String sessionURL = "http://www.broadinstitute.org/mmgp/textReader/IGV/mmrc_session.xml";
    static private String fileURL = "https://data.broadinstitute.org/igvdata/cshcourse/rwpe.washu.merged.bam";

    public static void main(String args[]) throws IOException {
        Socket socket = null;
        PrintWriter out = null;
        BufferedReader in = null;
        try {
            socket = new Socket("127.0.0.1", 60151);
            out = new PrintWriter(socket.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
            testMultiLocus(out, in);
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

        String cmd = "load https://data.broadinstitute.org/igvdata/tcga/gbmsubtypes/Broad.080528.subtypes.seg.gz";
        out.println(cmd);
        String response = in.readLine();
        System.out.println(cmd + " " + response);

        cmd = "goto EGFR PTEN";
        out.println(cmd);
        response = in.readLine();
        System.out.println(cmd + " " + response);


    }

}
