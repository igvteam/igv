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

import htsjdk.samtools.*;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.LineReader;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
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

    public CloseableIterator<PicardAlignment> iterator() {
        try {
            URL url = new URL(getQueryURL());
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);

            SAMFileReader reader = new SAMFileReader(new BufferedInputStream(is, 500000));
            reader.setValidationStringency(ValidationStringency.SILENT);
            CloseableIterator<SAMRecord> iter = reader.iterator();
            return new SAMQueryIterator(iter);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            return null;
        }
    }

    public CloseableIterator<PicardAlignment> query(String sequence, int start, int end, boolean contained) throws IOException {
        try {
            //
            final String parameters = "&chr=" + sequence + "&start=" + start + "&end=" + end +
                    "&contained=" + contained;
            //String encodedParameters = URLEncoder.encode(parameters);
            URL url = new URL(getQueryURL() + parameters);
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);

            SAMFileReader reader = new SAMFileReader(new BufferedInputStream(is, 500000));
            reader.setValidationStringency(ValidationStringency.SILENT);
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
        CloseableIterator<PicardAlignment> iter = reader.query(chr, start, end, false);
        while (iter.hasNext()) {
            Alignment a = iter.next();
            System.out.println(a);
        }
    }

}
