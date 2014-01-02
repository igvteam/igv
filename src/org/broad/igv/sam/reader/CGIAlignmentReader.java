/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTextHeaderCodec;
import net.sf.samtools.util.BufferedLineReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.LineReader;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SamAlignment;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Reader for querying a CGI service for SAM (note not BAM) files.
 *
 * @author Jim Robinson
 * @date 12/15/11
 * <p/>
 * Example URLs:
 * http://host/query.cgi?file=input.sam&start=800&end=900&contained=true
 * http://host/samHeader.cgi?file=input.sam
 * http://phhost/getSequenceNames.cgi?file=input.sam
 */

public class CGIAlignmentReader implements AlignmentReader {

    private static Logger log = Logger.getLogger(CGIAlignmentReader.class);

    String baseURL;
    String queryScript = "query.cgi";
    String headerScript = "samHeader.cgi";
    String seqNameScript = "getSequenceNames.cgi";
    String query;
    SAMFileHeader header;
    List<String> sequenceNames;

    public CGIAlignmentReader(String url) throws MalformedURLException {

        URL u = new URL(url);
        int port = u.getPort();
        baseURL = u.getProtocol() + "://" + u.getHost();
        if (port > 0) baseURL += ":" + port;
        baseURL += u.getPath();
        query = u.getQuery();

        loadHeader();
    }


    // The URL methods are package-scope to allow unit testing

    String getHeaderURL() {
        return baseURL.replace(queryScript, headerScript) + "?" + query;
    }

    String getSequenceNamesURL() {
        return baseURL.replace(queryScript, seqNameScript) + "?" + query;
    }

    String getQueryURL() {
        return baseURL + "?" + query;
    }


    public void close() throws IOException {
        //Nothing to do.  Could notify server that we are done, if that is useful
    }

    /**
     * Try to load header, if there are any errors just set to null and continue.  The header is optional
     */
    private void loadHeader() {
        InputStream is = null;
        try {
            URL url = new URL(getHeaderURL());
            is = HttpUtils.getInstance().openConnectionStream(url);

            LineReader reader = new BufferedLineReader(is);
            SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            header = codec.decode(reader, null);

        } catch (Exception ex) {
            log.info("Error loading header : " + ex.getMessage());
            header = null;
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException ex) {
                    log.error("Error closing url stream", ex);
                }
            }
        }
    }

    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            InputStream is = null;
            try {
                URL url = new URL(getSequenceNamesURL());
                is = HttpUtils.getInstance().openConnectionStream(url);
                BufferedReader br = new BufferedReader(new InputStreamReader(is));
                sequenceNames = new ArrayList<String>();
                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    String[] tokens = nextLine.split("\\s+");
                    for (String seq : tokens) {
                        sequenceNames.add(seq);
                    }
                }

            } catch (IOException e) {
                log.error("Error fetching sequence names", e);
                return null;
            } finally {
                if (is != null) try {
                    is.close();
                } catch (IOException e) {
                    log.error(e);
                }
            }
        }
        return sequenceNames;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }

    public SAMFileHeader getFileHeader() {
        if(header == null) {
            loadHeader();
        }
        return header;
    }

    public CloseableIterator<SamAlignment> iterator() {
        try {
            URL url = new URL(getQueryURL());
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);

            SAMFileReader reader = new SAMFileReader(new BufferedInputStream(is, 500000));
            reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            CloseableIterator<SAMRecord> iter = reader.iterator();
            return new SAMQueryIterator(iter);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return null;
        }
    }

    public CloseableIterator<SamAlignment> query(String sequence, int start, int end, boolean contained) throws IOException {
        try {
            //
            final String parameters = "&chr=" + sequence + "&start=" + start + "&end=" + end +
                    "&contained=" + contained;
            //String encodedParameters = URLEncoder.encode(parameters);
            URL url = new URL(getQueryURL() + parameters);
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);

            SAMFileReader reader = new SAMFileReader(new BufferedInputStream(is, 500000));
            reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            CloseableIterator<SAMRecord> iter = reader.iterator();
            return new SAMQueryIterator(sequence, start, end, contained, iter);

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return null;
        }
    }

    public boolean hasIndex() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public static void main(String[] args) throws IOException {
        Globals.setHeadless(true);
        CGIAlignmentReader reader = new CGIAlignmentReader("http://philtest.batcave.net/query.cgi?file=input.sam");
        String chr = "gi|66043271|ref|NC_007005.1|";
        int start = 800;
        int end = 900;
        CloseableIterator<SamAlignment> iter = reader.query(chr, start, end, false);
        while (iter.hasNext()) {
            Alignment a = iter.next();
            System.out.println(a);
        }
    }

}
