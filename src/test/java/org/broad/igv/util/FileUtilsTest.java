package org.broad.igv.util;

import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static org.junit.Assert.*;

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
