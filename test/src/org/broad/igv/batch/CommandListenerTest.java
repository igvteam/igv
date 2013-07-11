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
import org.broad.igv.util.HttpUtils;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
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

    private static String buildRootURL(){
        return String.format("http://localhost:%d/", port);
    }


    String genId = "hg18";

    @Test
    public void testGenomeSocket() throws Exception{
        String locus = "chr1:1-100";
        String gotoCmd = buildRootURL() + "goto?locus=" + locus;

        //This will URLencode

        //InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(gotoCmd));
        //BufferedReader in = new BufferedReader(new InputStreamReader(is));

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

        //This will URLencode
        InputStream is = HttpUtils.getInstance().openConnectionStream(new URL(cmd));
        //BufferedReader in = new BufferedReader(new InputStreamReader(is));

        assertEquals(genId, GenomeManager.getInstance().getGenomeId());
    }



}
