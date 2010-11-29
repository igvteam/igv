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

package org.broad.igv.util;

import org.junit.Test;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jun 9, 2010
 * Time: 8:39:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class EpigeneticsUtilsTest {

    @Test
    /**
     * fixedStep chrom=chr1 start=1 step=25 span=25
-1  0
-1  25
-1  50
1   75     75-125
1   100
2   125    125-200
2   150
2   175
1   200    200-250
1   225
-1  250
-1  275
-1  300
3   325    325-375
3   350
-1  375
-1  400
-1  425
4   450    450-525
4   475
4   500
     */
    public void testCreateCombinedWigs() throws Exception {

        File inputDir = new File("test/data");
        File outputDir = new File("./");
        EpigeneticsUtils.createCombinedWigs(inputDir, outputDir);
    }
}
