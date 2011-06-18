

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

package org.broad.igv.util;

import static org.junit.Assert.assertEquals;
import org.junit.Test;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jan 12, 2010
 * Time: 12:01:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class FileUtilsTest {


    @Test
    public void testGetRelaitvePath() {
        File basePath = new File("src");
        File targetPath = new File("lib/colt.jar");
        String relPath = "../lib/colt.jar";
        assertEquals(relPath, FileUtils.getRelativePath(basePath, targetPath));
    }


    @Test
    public void testLegalFileName() {
        String fn = "?[]/\\=+<>:;\"'*|";
        String legalFN = "_qm__fbr__rbr__fsl__bsl__eq__pl__lt__gt__co__sc__dq__sq__st__pp_";

        System.out.println(fn);

        String conFN = FileUtils.legalFileName(fn);

        System.out.println(conFN);

        assertEquals(legalFN, conFN);
    }

}
