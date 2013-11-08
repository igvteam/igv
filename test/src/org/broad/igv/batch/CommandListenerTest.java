/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.batch;

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



}
