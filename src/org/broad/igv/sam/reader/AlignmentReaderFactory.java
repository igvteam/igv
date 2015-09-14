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

package org.broad.igv.sam.reader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.ValidationStringency;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.ga4gh.Ga4ghAPIHelper;
import org.broad.igv.ga4gh.Ga4ghAlignmentReader;
import org.broad.igv.ga4gh.Ga4ghProvider;
import org.broad.igv.goby.GobyAlignmentQueryReader;
import org.broad.igv.util.FileUtils;
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
        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
    }

    public static AlignmentReader getReader(String path, boolean requireIndex) throws IOException {
        return getReader(new ResourceLocator(path), requireIndex);
    }

    public static AlignmentReader getReader(ResourceLocator locator) throws IOException {
        return getReader(locator, true);
    }

    public static AlignmentReader getReader(ResourceLocator locator, boolean requireIndex) throws IOException {
        log.debug("Getting alignment reader for " + locator);
        String pathLowerCase = locator.getPath().toLowerCase();

        AlignmentReader reader = null;

        String samFile = locator.getPath();
        String typeString = locator.getTypeString();

        if ("alist".equals(locator.getType())) {
            reader = getMergedReader(locator.getPath(), true);
        } else if (pathLowerCase.startsWith("http") && pathLowerCase.contains("/query.cgi?")) {
            reader = new CGIAlignmentReader(samFile);
        } else if (typeString.endsWith(".sam")) {
            reader = new SAMReader(samFile, requireIndex);

        } else if (typeString.endsWith(".aligned")
                || typeString.endsWith(".aligned.txt")
                || typeString.endsWith("bedz")
                || typeString.endsWith("bed")
                || typeString.endsWith("psl")
                || typeString.endsWith("pslx")) {
            reader = new GeraldReader(samFile, requireIndex);
        } else if (typeString.endsWith(".bam")) {
            if (locator.isLocal()) {
                reader = new BAMFileReader(new File(samFile));
            } else {
                try {
                    reader = new BAMHttpReader(locator, requireIndex);
                } catch (MalformedURLException e) {
                    log.error(e.getMessage(), e);
                    throw new DataLoadException("Error loading BAM file: " + e.toString(), locator.getPath());
                }

            }
        } else if (typeString.endsWith(".bam.list") || pathLowerCase.endsWith(".sam.list")) {
            reader = getBamListReader(locator.getPath(), requireIndex);
        } else if (GobyAlignmentQueryReader.supportsFileType(locator.getPath())) {
            try {
                reader = new GobyAlignmentQueryReader(locator.getPath());
            } catch (IOException e) {
                throw new RuntimeException("Cannot load Goby alignment " + locator.getPath(), e);

            }
        } else if (Ga4ghAlignmentReader.supportsFileType(locator.getType())) {
            Ga4ghProvider provider = (Ga4ghProvider) locator.getAttribute("provider");
            return new Ga4ghAlignmentReader(provider, locator.getPath());
        }

        else {
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
                    if(f.length() == 0) continue;  // Empty line
                    
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

