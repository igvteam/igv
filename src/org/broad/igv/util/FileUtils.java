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

import org.broad.igv.util.ftp.FTPUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.util.*;
import java.util.regex.Pattern;

/**
 * @author jrobinso
 */
public class FileUtils {

    private static final Logger log = Logger.getLogger(FileUtils.class);

    final public static String LINE_SEPARATOR = System.getProperty("line.separator");
    final public static String FILE_SEP = System.getProperty("file.separator");
    final static String[] igvJnlpPrefixes = {"igv", "ichip", "29mammals", "hic"};



    public static boolean resourceExists(String path) {
        try {
            if (isRemote(path)) {
                return HttpUtils.getInstance().resourceAvailable(new URL(path));
            } else {
                return (new File(path)).exists();
            }
        } catch (IOException e) {
            log.error("Error checking existence of: " + path, e);
            return false;
        }
    }


    public static boolean isRemote(String path) {
        if (path == null) {
            return false;
        }
        return path.startsWith("http://") || path.startsWith("https://") ||
                path.startsWith("ftp://");
    }


    public static boolean canWriteTo(File file) {
        FileOutputStream fos = null;
        try {
            file.createNewFile();
            return true;

            // Has permission
        } catch (Exception e) {
            return false;
        } finally {
            file.delete();
        }
    }

    public static String getInstallDirectory() {

        String path = FileUtils.class
                .getProtectionDomain()
                .getCodeSource()
                .getLocation().getPath();
        File f = new File(path);
        if (f.isDirectory()) {
            return f.getAbsolutePath();
        } else {
            return f.getParentFile().getAbsolutePath();
        }
    }

    /**
     * Get the relative path from one file to another, specifying the directory separator.
     * If one of the provided resources does not exist, it is assumed to be a file unless it ends with '/' or
     * '\'.
     *
     * @param targetPath    targetPath is calculated to this file
     * @param basePath      basePath is calculated from this file
     * @param pathSeparator directory separator. The platform default is not assumed so that we can test Unix behaviour when running on Windows (for example)
     * @return
     */
    public static String getRelativePath(String basePath, String targetPath, String pathSeparator) {

        // Normalize the paths
        String normalizedTargetPath = FilenameUtils.normalizeNoEndSeparator(targetPath);
        String normalizedBasePath = FilenameUtils.normalizeNoEndSeparator(basePath);

        // Undo the changes to the separators made by normalization
        if (pathSeparator.equals("/")) {
            normalizedTargetPath = FilenameUtils.separatorsToUnix(normalizedTargetPath);
            normalizedBasePath = FilenameUtils.separatorsToUnix(normalizedBasePath);

        } else if (pathSeparator.equals("\\")) {
            normalizedTargetPath = FilenameUtils.separatorsToWindows(normalizedTargetPath);
            normalizedBasePath = FilenameUtils.separatorsToWindows(normalizedBasePath);

        } else {
            throw new IllegalArgumentException("Unrecognised dir separator '" + pathSeparator + "'");
        }

        String[] base = normalizedBasePath.split(Pattern.quote(pathSeparator));
        String[] target = normalizedTargetPath.split(Pattern.quote(pathSeparator));

        // First get all the common elements. Store them as a string,
        // and also count how many of them there are.
        StringBuffer common = new StringBuffer();

        int commonIndex = 0;
        while (commonIndex < target.length && commonIndex < base.length
                && target[commonIndex].equals(base[commonIndex])) {
            common.append(target[commonIndex] + pathSeparator);
            commonIndex++;
        }

        if (commonIndex == 0) {
            // No single common path element. This most
            // likely indicates differing drive letters, like C: and D:.
            // These paths cannot be relativized.
            return targetPath;
        }

        // The number of directories we have to backtrack depends on whether the base is a file or a dir
        // For example, the relative path from
        //
        // /foo/bar/baz/gg/ff to /foo/bar/baz
        //
        // ".." if ff is a file
        // "../.." if ff is a directory
        //
        // The following is a heuristic to figure out if the base refers to a file or dir. It's not perfect, because
        // the resource referred to by this path may not actually exist, but it's the best I can do
        boolean baseIsFile = true;

        File baseResource = new File(normalizedBasePath);

        if (baseResource.exists()) {
            baseIsFile = baseResource.isFile();

        } else if (basePath.endsWith(pathSeparator)) {
            baseIsFile = false;
        }

        StringBuffer relative = new StringBuffer();

        if (base.length != commonIndex) {
            int numDirsUp = baseIsFile ? base.length - commonIndex - 1 : base.length - commonIndex;

            for (int i = 0; i < numDirsUp; i++) {
                relative.append(".." + pathSeparator);
            }
        }
        relative.append(normalizedTargetPath.substring(common.length()));
        return relative.toString();
    }

    static public String getPlatformIndependentPath(String path) {
        return path.replace('\\', '/');
    }


    // Convenience method
    public static String getRelativePath(String basePath, String targetPath) {
        return getRelativePath(basePath, targetPath, FILE_SEP);
    }

    /**
     * Delete a directory and all subdirectories
     *
     * @param dir
     * @return
     */
    public static boolean deleteDir(File dir) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i = 0; i < children.length; i++) {
                boolean success = deleteDir(new File(dir, children[i]));
                if (!success) {
                    return false;
                }
            }
        }

        // The directory is now empty so delete it
        return dir.delete();
    }


    /**
     * Test to see if the first comment line (first line not starting with #) is tab-delimited with the
     * given number of minimum columns.  Limit the test to the first 1,000 lines.
     */
    public static boolean isTabDelimited(ResourceLocator loc, int minColumnCount) throws IOException {

        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(loc.getPath());
            int nLinesTested = 0;
            String nextLine;
            while ((nextLine = reader.readLine()) != null && nLinesTested < 1000) {
                if (nextLine.startsWith("#")) {
                    continue;
                }
                nLinesTested++;
                String[] tokens = nextLine.split("\t");
                if (tokens.length >= minColumnCount) {
                    return true;
                }
            }
            return false;
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }


    /**
     * Copy a file from one location to another, using buffered writing
     *
     * @param inputFile
     * @param outputFile
     * @throws java.io.IOException
     */

    public static void copyFile(File inputFile, File outputFile) throws IOException {

        OutputStream out = null;
        InputStream in = null;
        try {
            in = new FileInputStream(inputFile);
            out = new FileOutputStream(outputFile);
            byte[] buffer = new byte[64000];
            int bytes_read;
            while ((bytes_read = in.read(buffer)) != -1) {
                out.write(buffer, 0, bytes_read);
            }

        } catch (Exception e) {
            outputFile.delete();
            throw new RuntimeException("<html>Error copying file: " + outputFile.getAbsoluteFile() +
                    "<br/>" + e.toString());

        } finally {
            if (in != null) {
                in.close();
            }
            if (out != null) {
                out.flush();
                out.close();
            }
        }
    }

    private static void cleanup(File dir) {

        if (dir.exists() && dir.isDirectory()) {
            File[] jnlpFiles = dir.listFiles(new FileFilter() {

                public boolean accept(File arg0) {
                    final String name = arg0.getName();
                    for (String pre : igvJnlpPrefixes) {
                        if (name.startsWith(pre) && name.endsWith(".jnlp")) {
                            return true;
                        }
                    }
                    return false;
                }
            });

            // Sort files by ascending version number
            Arrays.sort(jnlpFiles, new Comparator<File>() {

                public int compare(File file1, File file2) {
                    if (org.apache.commons.io.FileUtils.isFileNewer(file1, file2)) {
                        return 1;
                    } else {
                        return -1;
                    }
                }
            });

            // Delete all but the highest version (newest) jnlp file
            for (int i = 0; i < jnlpFiles.length - 1; i++) {
                jnlpFiles[i].delete();
            }

            // Strip the version nuber fro the newest file
            if (jnlpFiles.length > 1) {
                File newestFile = jnlpFiles[jnlpFiles.length - 1];
                String fn = newestFile.getName();
                int dotIndex = fn.indexOf(".jnlp");
                int dashIndex = fn.lastIndexOf("-");
                if (dashIndex > 1) {
                    String newName = fn.substring(0, dashIndex) + fn.substring(dotIndex);
                    newestFile.renameTo(new File(newestFile.getParentFile(), newName));
                }
            }
        }
    }


//    public static void main(String[] args) throws IOException {
//        File inputDirectory = new File(args[0]);
//        File outputDirectory = new File(args[1]);
//        searchAndReplace(inputDirectory, outputDirectory, args[2], args[3]);
//    }

    public static String getFileExtension(String filePath) {
        String extension = null;
        int indexOfExtension = filePath.lastIndexOf(".");
        if (indexOfExtension >= 0) {
            extension = filePath.substring(indexOfExtension, filePath.length());
        }
        return extension;
    }

    /**
     * Returns an absolute path from the parent directory
     * of {@code referencePath} to the sub-element {@code inputPath}.
     * Safe to use with URLs.
     * If {@code inputPath} is an absolute path, it is returned unaltered.
     * <br/>
     * e.g.<br/>
     * String absPath = FileUtils.getAbsolutePath("test/mysession.xml", "/Users/bob/data/otherdata.xml");
     * System.out.println(absPath);<br/>
     * >>>> /Users/bob/data/test/mysession.xml
     *
     * @param inputPath     Relative path element
     * @param referencePath Absolute path root
     * @return
     */
    public static String getAbsolutePath(String inputPath, String referencePath) {
        if (isRemote(inputPath)) {
            return inputPath;
        }
        File inFile = new File(inputPath);
        if (inFile.isAbsolute()) {
            return inFile.getAbsolutePath();
        }

        String absolutePath;

        if (isRemote(referencePath)) {
            int idx = referencePath.lastIndexOf("/");
            String basePath = referencePath.substring(0, idx);
            absolutePath = basePath + "/" + inputPath;
        } else {
            File parent = new File(referencePath).getParentFile();
            File file = new File(parent, inputPath);
            try {
                absolutePath = file.getCanonicalPath();
            } catch (IOException e) {
                absolutePath = file.getAbsolutePath();
            }
        }

        return absolutePath;
    }

    /**
     * Return the path path.  The trailing "/" is not included.
     *
     * @param path
     * @return
     */
    public static String getParent(String path) {
        String piPath = getPlatformIndependentPath(path);
        int lastSlashIdx = piPath.lastIndexOf("/");
        return lastSlashIdx <= 0 ? path : path.substring(0, lastSlashIdx);
    }

    /**
     * Checks the system path for the provided executable.
     * If {@code executable} is a path (contains a path separator)
     * then it is returned unaltered
     *
     * @param executable
     * @return
     */
    public static String findExecutableOnPath(String executable) {
        if (executable.contains(File.separator)) return executable;

        String systemPath = System.getenv("PATH");
        if (systemPath == null) systemPath = System.getenv("path");
        if (systemPath == null || File.pathSeparator == null) return executable;

        String[] pathDirs = systemPath.split(File.pathSeparator);

        String fullPath = executable;
        for (String pathDir : pathDirs) {
            File file = new File(pathDir, executable);
            if (file.isFile()) {
                fullPath = file.getAbsolutePath();
                break;
            }
        }
        return fullPath;
    }



    /**
     * Return the length of the file, which might be remote.
     *
     * @param file
     * @return
     */
    public static long getLength(String file) {

        if (isRemote(file)) {
            try {
                URL url = new URL(file);
                if (file.startsWith("ftp://")) {
                    return FTPUtils.getContentLength(url);
                } else {
                    return HttpUtils.getInstance().getContentLength(url);
                }
            } catch (Exception e) {
                log.error("Error fetching content length for: " + file, e);
                return -1;
            }
        } else {
            File f = new File(file);
            if (f.exists() && f.isFile()) {
                return f.length();
            } else {
                return -1;
            }
        }
    }
}
