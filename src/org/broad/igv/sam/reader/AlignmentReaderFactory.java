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
package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.goby.GobyAlignmentQueryReader;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.*;

/**
 * @author jrobinso
 */
public class AlignmentReaderFactory {
    private static Logger log = Logger.getLogger(AlignmentReaderFactory.class);

    static {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
    }

    public static AlignmentReader getReader(String path, boolean requireIndex) throws IOException {
        return getReader(new ResourceLocator(path), requireIndex);
    }

    public static AlignmentReader getReader(ResourceLocator locator) throws IOException {
        return getReader(locator, true);
    }

    public static AlignmentReader getReader(ResourceLocator locator, boolean requireIndex) throws IOException {

        String pathLowerCase = locator.getPath().toLowerCase();

        AlignmentReader reader = null;

        String samFile = locator.getPath();

        if ("alist".equals(locator.getType())) {
            reader = getMergedReader(locator.getPath(), true);
        } else if (pathLowerCase.startsWith("http") && pathLowerCase.contains("/query.cgi?")) {
            reader = new CGIAlignmentReader(samFile);
        } else if (pathLowerCase.endsWith(".sam")) {
            reader = new SAMReader(samFile, requireIndex);

        } else if (pathLowerCase.endsWith("sorted.txt")
                || pathLowerCase.endsWith(".aligned")
                || pathLowerCase.endsWith(".aligned.txt")
                || pathLowerCase.endsWith("bedz")
                || pathLowerCase.endsWith("bed")
                || pathLowerCase.endsWith("psl")
                || pathLowerCase.endsWith("pslx")) {
            reader = new GeraldReader(samFile, requireIndex);
        } else if (pathLowerCase.endsWith(".bam")) {
            if (locator.isLocal()) {
                reader = new BAMFileReader(new File(samFile));
            } else if (HttpUtils.getInstance().isURL(locator.getPath().toLowerCase())) {
                try {
                    reader = new BAMHttpReader(locator, requireIndex);
                } catch (MalformedURLException e) {
                    log.error("", e);
                    throw new DataLoadException("Error loading BAM file: " + e.toString(), locator.getPath());
                }

            } else {
                if (locator.getServerURL() != null) {
                    reader = new BAMWebserviceReader(locator);
                }
            }
        } else if (pathLowerCase.endsWith(".bam.list")) {
            if (locator.getServerURL() != null) {
                reader = new BAMWebserviceReader(locator);
            } else {
                reader = getBamListReader(locator.getPath(), requireIndex);
            }
        } else if (GobyAlignmentQueryReader.supportsFileType(locator.getPath())) {
            try {
                reader = new GobyAlignmentQueryReader(locator.getPath());
            } catch (IOException e) {
                throw new RuntimeException("Cannot load Goby alignment " + locator.getPath(), e);

            }
        } else {
            throw new RuntimeException("Cannot find reader for aligment file: " + locator.getPath());
        }

        return reader;
    }

    static AlignmentReader getBamListReader(String listFile, boolean requireIndex) {

        List<AlignmentReader> readers = new ArrayList();
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(listFile);
            Map<String, String> replacements = new HashMap();
            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("#replace")) {
                    String[] tokens = nextLine.split("\\s+");
                    if (tokens.length == 2) {
                        String[] kv = tokens[1].split("=");
                        if (kv.length == 2)
                            replacements.put(kv[0], kv[1]);
                    }
                } else {
                    String f = nextLine.trim();
                    for (Map.Entry<String, String> entry : replacements.entrySet()) {
                        f = f.replace(entry.getKey(), entry.getValue());
                    }

                    f = FileUtils.getAbsolutePath(f, listFile);
                    readers.add(AlignmentReaderFactory.getReader(f, requireIndex));
                }
            }
            if (readers.size() == 1) {
                return readers.get(0);
            } else {
                return new MergedAlignmentReader(readers);
            }
        } catch (IOException e) {
            log.error("Error parsing " + listFile, e);
            throw new RuntimeException("Error parsing: " + listFile + " (" + e.toString() + ")");
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
    }

    public static AlignmentReader getMergedReader(String alignmentFileList, boolean requireIndex) {

        String aFile = null;
        try {
            String[] alignmentFiles = ParsingUtils.COMMA_PATTERN.split(alignmentFileList);
            List<AlignmentReader> readers = new ArrayList(alignmentFiles.length);
            for (String f : alignmentFiles) {
                aFile = f;
                readers.add(AlignmentReaderFactory.getReader(aFile, requireIndex));
            }
            if (readers.size() == 1) {
                return readers.get(0);
            } else {
                return new MergedAlignmentReader(readers);
            }
        } catch (IOException e) {
            log.error("Error instantiating reader for: " + aFile, e);
            throw new RuntimeException("Error instantiating reader for : " + aFile + " (" + e.toString() + ")");
        }

    }

    /**
     * @param header
     * @return Return the set of platforms, uppercase. Will be null iff header is null
     */
    public static Set<String> getPlatforms(SAMFileHeader header) {
        Set<String> platforms = null;
        if (header != null) {
            List<SAMReadGroupRecord> readGroups = header.getReadGroups();
            if (readGroups != null) {
                platforms = new HashSet<String>();
                for (SAMReadGroupRecord rg : readGroups) {
                    String platform = rg.getPlatform();
                    if (platform != null) {
                        platforms.add(platform.toUpperCase());
                    }
                }
            }
        }
        return platforms;
    }

}

