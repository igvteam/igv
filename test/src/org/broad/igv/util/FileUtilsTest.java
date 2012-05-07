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

import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jan 12, 2010
 * Time: 12:01:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class FileUtilsTest {


    @Test
    public void testFindRelativePath() throws IOException {
        String sep = System.getProperty("file.separator");
        File basePath = new File("src");
        File targetPath = new File("lib" + sep + "vcf.jar");
        String relPath = FileUtils.getRelativePath(basePath.getAbsolutePath(), targetPath.getAbsolutePath(), sep);
        assertEquals( ".." + sep + "lib" + sep + "vcf.jar", relPath);
    }

    @Test
    public void testGetRelativePathsUnix() {
        assertEquals("stuff/xyz.dat", FileUtils.getRelativePath("/var/data/", "/var/data/stuff/xyz.dat",  "/"));
        assertEquals("../../b/c", FileUtils.getRelativePath( "/a/x/y/","/a/b/c", "/"));
        assertEquals("../../b/c", FileUtils.getRelativePath( "/m/n/o/a/x/y/", "/m/n/o/a/b/c", "/"));
    }

    @Test
    public void testGetRelativePathWindows() {
        String base = "C:\\Windows\\Boot\\Fonts\\chs_boot.ttf";
        String target = "C:\\Windows\\Speech\\Common\\sapisvr.exe";

        String relPath = FileUtils.getRelativePath(target, base, "\\");
        assertEquals("..\\..\\Boot\\Fonts\\chs_boot.ttf", relPath);
    }

    @Test
    public void testGetRelativePathDirectoryToFile() {
        String base = "C:\\Windows\\Boot\\Fonts\\chs_boot.ttf";
        String target = "C:\\Windows\\Speech\\Common\\";

        String relPath = FileUtils.getRelativePath(target, base, "\\");
        assertEquals("..\\..\\Boot\\Fonts\\chs_boot.ttf", relPath);
    }

    @Test
    public void testGetRelativePathFileToDirectory() {
        String base = "C:\\Windows\\Boot\\Fonts";
        String target = "C:\\Windows\\Speech\\Common\\foo.txt";

        String relPath = FileUtils.getRelativePath(target, base, "\\");
        assertEquals("..\\..\\Boot\\Fonts", relPath);
    }

    @Test
    public void testGetRelativePathDirectoryToDirectory() {
        String base = "C:\\Windows\\Boot\\";
        String target = "C:\\Windows\\Speech\\Common\\";
        String expected = "..\\..\\Boot";

        String relPath = FileUtils.getRelativePath(target, base, "\\");
        assertEquals(expected, relPath);
    }

    @Test
    public void testGetRelativePathDifferentDriveLetters() {
        String base = "D:\\sources\\recovery\\RecEnv.exe";
        String target = "C:\\Java\\workspace\\AcceptanceTests\\Standard test data\\geo\\";
        String expected = target;

        String relPath = FileUtils.getRelativePath(base, target, "\\");
        assertEquals(expected, relPath);
    }



}
