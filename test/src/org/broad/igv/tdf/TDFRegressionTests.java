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

package org.broad.igv.tdf;

import org.junit.AfterClass;

import static org.junit.Assert.assertEquals;

import org.junit.BeforeClass;
import org.junit.Test;


/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
public class TDFRegressionTests {


     /**
     * IGV-1421 - error reading a version 1 file (error thrown in header)
     */
    @Test
    public void test_IGV1417() {
        String tdfFile = "http://www.broadinstitute.org/igvdata/annotations/hg18/conservation/pi.ewig.tdf";
        TDFReader reader = TDFReader.getReader(tdfFile);
        assertEquals(1, reader.getVersion());
    }

    /**
     * IGV-1421 - error reading a version 2 file (error thrown in header)
     */
    @Test
    public void test_IGV1421() {
        String tdfFile = "http://www.broadinstitute.org/igvdata/1KG/freeze5/NA12878.pilot2.454.bam.tdf";
        TDFReader reader = TDFReader.getReader(tdfFile);
        assertEquals(2, reader.getVersion());
    }
}
