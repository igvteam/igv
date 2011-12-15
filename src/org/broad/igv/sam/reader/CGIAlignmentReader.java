package org.broad.igv.sam.reader;

import net.sf.samtools.*;
import net.sf.samtools.util.BufferedLineReader;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.LineReader;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.SamAlignment;
import org.broad.igv.util.HttpUtils;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date 12/15/11
 * <p/>
 * http://philtest.batcave.net/query.cgi?file=input.sam&start=800&end=900&contained=true
 * http://philtest.batcave.net/samHeader.cgi?file=input.sam
 * http://philtest.batcave.net/getSequenceNames.cgi?file=input.sam
 * <p/>
 * host: philtest.batcave.net
 * query: file=input.sam&start=800&end=900&contained=true
 */
public class CGIAlignmentReader implements AlignmentQueryReader {

    private static Logger log = Logger.getLogger(CGIAlignmentReader.class);

    String baseURL;
    String queryPath = "/query.cgi?";
    String headerPath = "/samHeader.cgi?";
    String seqNamePath = "/getSequenceNames.cgi?";
    String file;
    SAMFileHeader header;

    public CGIAlignmentReader(String url) throws MalformedURLException {

        URL u = new URL(url);
        baseURL = u.getProtocol() + "://" + u.getHost();
        file = u.getQuery();
        loadHeader();
    }

    public void close() throws IOException {
        //Nothing to do.  Could notify server that we are done, if that is useful
    }

    // TODO -- nearly exact copy of BAMRemoteQueryReader, only url differs. -- refactor
    private void loadHeader() {
        InputStream is = null;
        try {
            URL url = new URL(baseURL + headerPath + file);
            is = HttpUtils.getInstance().openConnectionStream(url);

            LineReader reader = new BufferedLineReader(is);
            SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            header = codec.decode(reader, null);

        } catch (IOException ex) {
            log.error("Error opening file", ex);
            throw new RuntimeException(ex);
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

    public Set<String> getSequenceNames() {

        InputStream is = null;
        try {
            URL url = new URL(baseURL + seqNamePath + file);
            is = HttpUtils.getInstance().openConnectionStream(url);
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            HashSet<String> seqNames = new HashSet<String>();
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = nextLine.split("\\s+");
                for (String seq : tokens) {
                    seqNames.add(seq);
                }
            }
            return seqNames;

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

    public SAMFileHeader getHeader() throws IOException {
        return header;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public CloseableIterator<Alignment> iterator() {
        try {
            URL url = new URL(baseURL + queryPath + file);
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

    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) throws IOException {
        try {
            //
            final String parameters = file +  "&chr=" + sequence + "&start=" + start + "&end=" + end +
                    "&contained=" + contained;
            String encodedParameters = URLEncoder.encode(parameters);
            URL url = new URL(baseURL + queryPath + encodedParameters);
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
        CloseableIterator<Alignment> iter =  reader.query(chr, start, end, false);
        while(iter.hasNext()) {
            Alignment a = iter.next();
            System.out.println(a);
        }
    }

}
