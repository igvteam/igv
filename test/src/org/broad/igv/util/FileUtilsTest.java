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
    public void testFindRelativeHttpPath() throws IOException {
        String basePath = "http://foo.bar.com/baseDir/";
        String targetPath = "http://foo.bar.com/baseDir/dir/test.txt";
        String relPath = FileUtils.getRelativePath(basePath, targetPath, "/");
       System.out.println(relPath);
    }

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

    @Test
    public void testGetParent() throws Exception {

        String windowsPath = "C:\\ path with spaces\\subdir\\foo.txt";
        assertEquals("C:\\ path with spaces\\subdir", FileUtils.getParent(windowsPath));

        String unixPath = "/path with spaces/subdir/foo.txt";
        assertEquals("/path with spaces/subdir", FileUtils.getParent(unixPath));

        String httpPath = "http://host:port/subdir/foo.txt";
        assertEquals("http://host:port/subdir", FileUtils.getParent(httpPath));
    }

    @Test
    public void testGetAbsolutePath() throws Exception {

        String inputPath = "test/mysession.xml";
        String referencePath = "http://foo.bar.com/bob/data/otherdata.xml";
        String absolutePath = "http://foo.bar.com/bob/data/test/mysession.xml";
        assertEquals(absolutePath, FileUtils.getAbsolutePath(inputPath, referencePath));

    }

}
