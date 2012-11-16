/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.session;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.net.URL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * @author jrobinso
 */
public class SessionManagerTest {

    public SessionManagerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of getInstance method, of class SessionManager.
     */
    @Test
    public void getInstance() {

        IGVSessionReader.SessionElement elem = IGVSessionReader.SessionElement.findEnum("COLOR_SCALE");
        IGVSessionReader.SessionElement expResult = IGVSessionReader.SessionElement.COLOR_SCALE;
        assertEquals(expResult, elem);
    }

    /**
     * Test of loadSession method, of class SessionManager.
     */
    //@Test
    public void loadSession() {
        System.out.println("loadSession");
        URL url = null;
        Session session = null;
        IGVSessionReader instance = null;
        //instance.loadSession(url, session);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of saveSession method, of class SessionManager.
     */
    //@Test
    public void saveSession() throws Exception {
        System.out.println("saveSession");
        Session session = null;
        File outputFile = null;
        IGVSessionReader instance = null;
        //instance.saveSession(session, outputFile);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createXmlFromSession method, of class SessionManager.
     */
    //@Test
    public void createXmlFromSession() {
        System.out.println("createXmlFromSession");
        Session session = null;
        File outputFile = null;
        IGVSessionReader instance = null;
        String expResult = "";
        //String result = instance.createXmlFromSession(session, outputFile);
        //assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

}