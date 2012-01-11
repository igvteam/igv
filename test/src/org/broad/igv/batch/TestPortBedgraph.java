package org.broad.igv.batch;

import org.broad.igv.util.TestUtils;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;

/**
 * @author Jim Robinson
 * @date 11/30/11
 */
public class TestPortBedgraph {

    private PrintWriter out;
    private BufferedReader in;

    TestPortBedgraph() {
        try {
            Socket socket = new Socket("127.0.0.1", 60151);
            out = new PrintWriter(socket.getOutputStream(), true);
            in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
        } catch (Exception E) {
            System.out.println("IO exception");
        }
    }

    public void importBedGraph(String filename) {
        try {

            out.println("load " + filename);
            String response = in.readLine();
            System.out.println(response);
        } catch (Exception E) {
            System.out.println("IO exception");
        }
    }

    public static void main(String[] args) {
        String testfile = TestUtils.DATA_DIR + "/wig/jira_1409.bedgraph";
        TestPortBedgraph test = new TestPortBedgraph();
        test.importBedGraph(testfile);
    }

}
