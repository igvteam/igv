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

import biz.source_code.base64Coder.Base64Coder;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.HttpURLConnection;
import java.net.Socket;
import java.net.URL;

import static org.junit.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Jul-11
 */
public class CommandListenerTest extends AbstractHeadedTest {

    private static final int port = 60151;

    @BeforeClass
    public static void setUpClass() throws Exception {
        PreferenceManager.getInstance().override(PreferenceManager.PORT_ENABLED, "true");
        PreferenceManager.getInstance().override(PreferenceManager.PORT_NUMBER, "" + port);
        AbstractHeadedTest.setUpClass();
    }

    @Before
    public void setUp() throws Exception{
        super.setUp();
        CommandListener.halt();
        CommandListener.start(port);
        IGV.getInstance().loadGenome(TestUtils.defaultGenome, null, true);
    }

    @After
    public void tearDown() throws Exception{
        super.tearDown();
        CommandListener.halt();
    }

    private static String buildRootURL(){
        return String.format("http://localhost:%d/", port);
    }


    String genId = "mm10";

    @Test
    public void testGenomeSocket() throws Exception{
        String locus = "chr1:1-100";

        Socket socket = new Socket("localhost", port);
        PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
        BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

        out.println("genome " + genId);

        String response = in.readLine();
        System.out.println(response);

        assertEquals(genId, GenomeManager.getInstance().getGenomeId());
    }

    @Test
    public void testGenomeLink() throws Exception{
        String cmd = buildRootURL() + "load?genome=" + genId;
        connect(cmd);

        assertEquals(genId, GenomeManager.getInstance().getGenomeId());
    }

    @Test
    public void testLoadURLLink() throws Exception{

        String urlPath = CommandExecutorTest.urlPathSpaces;
        String name = "mytestfile";
        String cmd = buildRootURL() + "load?file=" + urlPath + "&name=" + name;
        connect(cmd);

        TestUtils.assertTrackLoaded(igv, name);
    }

    @Test
    public void testLoadURLSocket() throws Exception{

        String urlPath = CommandExecutorTest.urlPathSpaces;
        String name = "mytestfile";

        Socket socket = new Socket("localhost", port);
        PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
        BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

        out.println("load " + urlPath + " name=" + name);

        String response = in.readLine();
        System.out.println(response);

        TestUtils.assertTrackLoaded(igv, name);
    }

    @Test
    public void testLoadFileSpacesSocket() throws Exception{
        tstLoadFileSocket(CommandExecutorTest.dirPathSpaces, CommandExecutorTest.fileName01);
    }

    @Test
    public void testLoadFileSpacesPercSocket() throws Exception{
        tstLoadFileSocket(CommandExecutorTest.dirPathSpaces, CommandExecutorTest.fileNamePerc);
    }

    private void tstLoadFileSocket(String fidir, String finame) throws Exception{
        Socket socket = new Socket("localhost", port);
        PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
        BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

        String fileString = String.format("\"%s\"", (new File(fidir, finame)).getPath());
        out.println("load " + fileString);

        String response = in.readLine();
        System.out.println(response);

        TestUtils.assertTrackLoaded(IGV.getInstance(), finame);
    }

    @Test
    public void testLoadFileSpacesLink() throws Exception{
        tstLoadFileLink(CommandExecutorTest.dirPathSpaces, CommandExecutorTest.fileName01);
    }

    @Test
    public void testLoadFileSpacesPercLink() throws Exception{
        tstLoadFileLink(CommandExecutorTest.dirPathSpaces, CommandExecutorTest.fileNamePerc);
    }

    private void tstLoadFileLink(String fidir, String finame) throws Exception{
        String fileString = (new File(fidir, finame)).getPath();
        String urlPath = StringUtils.encodeURL(fileString);
        String name = "mytestfile";

        String cmd = buildRootURL() + "load?file=" + urlPath + "&name=" + name;
        connect(cmd);

        TestUtils.assertTrackLoaded(IGV.getInstance(), name);
    }

    private HttpURLConnection connect(String urlStr) throws Exception{
        URL url = new URL(urlStr);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("GET");
        conn.setRequestProperty("Connection", "Keep-Alive");
        conn.connect();
        System.out.println(conn.getResponseCode() + ":" + conn.getResponseMessage());
        return conn;
    }

    @Test
    public void testSHA1() throws Exception {

        // Example from WebSocket RFC
        String guid = "258EAFA5-E914-47DA-95CA-C5AB0DC85B11";
        String clientSocketKey = "dGhlIHNhbXBsZSBub25jZQ==";
        String serverSocketKey = "s3pPLMBiTxaQ9kYGzzhZRbK+xOo=";

        String generatedKey = CommandListener.computeResponseKey(clientSocketKey + guid);

        assertEquals(new String(generatedKey), serverSocketKey);


    }



}
