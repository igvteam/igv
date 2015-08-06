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