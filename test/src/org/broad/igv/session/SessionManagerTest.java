/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.session;

import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.net.URL;

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

        SessionReader.SessionElement elem = SessionReader.SessionElement.findEnum("COLOR_SCALE");
        SessionReader.SessionElement expResult = SessionReader.SessionElement.COLOR_SCALE;
        assertEquals(expResult, elem);
        // TODO review the generated test code and remove the default call to fail.
    }

    /**
     * Test of loadSession method, of class SessionManager.
     */
    //@Test
    public void loadSession() {
        System.out.println("loadSession");
        URL url = null;
        Session session = null;
        SessionReader instance = null;
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
        SessionReader instance = null;
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
        SessionReader instance = null;
        String expResult = "";
        //String result = instance.createXmlFromSession(session, outputFile);
        //assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

}