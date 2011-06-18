/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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
 * HDF5UtilsTest.java
 * JUnit based test
 *
 * Created on November 15, 2007, 11:39 AM
 */
package util;

import junit.framework.TestCase;

/**
 * @author jrobinso
 */
public class HDF5UtilsTest extends TestCase {

    public HDF5UtilsTest(String testName) {
        super(testName);
    }

    protected void setUp() throws Exception {
    }

    protected void tearDown() throws Exception {
    }

    /**
     * Test of writeAllData method, of class util.HDF5Utils.
     */
    public void testWriteData() {


        //HDF5Utils.writeAllData(datasetId, type, data);

        // TODO review the generated test code and remove the default call to fail.

    }

    /**
     * Test of writeAttribute method, of class util.HDF5Utils.
     */
    public void testWriteAttribute() {
        System.out.println("writeAttribute");

        int locId = 0;
        String name = "";
        Object value = null;

        //int result = HDF5Utils.writeAttribute(locId, name, value);
    }

    /**
     * Test of readDoubleAttribute method, of class util.HDF5Utils.
     */
    public void testReadDoubleAttribute() {

        int locId = 0;
        String name = "";
        //double result = HDF5Utils.readDoubleAttribute(locId, name);

    }

    /**
     * Test of readIntegerAttribute method, of class util.HDF5Utils.
     */
    public void testReadIntegerAttribute() {
        System.out.println("readIntegerAttribute");

        int locId = 0;
        String name = "";
        //int result = HDF5Utils.readIntegerAttribute(locId, name);
    }

    /**
     * Test of readLongAttribute method, of class util.HDF5Utils.
     */
    public void testReadLongAttribute() {
        System.out.println("readLongAttribute");

        int locId = 0;
        String name = "";
        //HDF5Utils.readLongAttribute(locId, name);

    }

    /**
     * Test of readAllFloats method, of class util.HDF5Utils.
     */
    public void testReadAllFloats() {
        System.out.println("readAllFloats");

        int datasetId = 0;

        float[] expResult = null;
        //float[] result = HDF5Utils.readAllFloats(datasetId);
    }

    /**
     * Test of readAllInts method, of class util.HDF5Utils.
     */
    public void testReadAllInts() {
        System.out.println("readAllInts");

        int locId = 0;
        String dsName = "";

        int[] expResult = null;
        // int[] result = HDF5Utils.readAllInts(locId, dsName);
    }

    /**
     * Test of readFloats method, of class util.HDF5Utils.
     */
    public void testReadFloats() {
        System.out.println("readFloats");

        int fileId = 0;
        String dsName = "";
        int fromIndex = 0;
        int toIndex = 0;

        float[] expResult = null;
        //float[] result = HDF5Utils.readFloats(fileId, dsName, fromIndex, toIndex);
    }


    /**
     * Test of getChildNames method, of class util.HDF5Utils.
     */
    public void testGetChildNames() {
        System.out.println("getChildNames");

        int fileId = 0;
        String groupPath = "";

        String[] expResult = null;
        //String[] result = HDF5Utils.getChildNames(fileId, groupPath);
    }
}
