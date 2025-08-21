package org.broad.igv.htsget;

import com.google.gson.*;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.apache.commons.compress.compressors.gzip.GzipUtils;
import org.broad.igv.util.HttpUtils;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * A low-level class to fetch header and data from htsget servers.  No intrepretation or parsing of the header or
 * data is done here, the class is designed as a helpfer for format specific readers.  Currently this only supports the
 * variant (VCF) format, the htsjdk is used for BAM and CRAM which does not use this helper.
 */
class HtsgetReader {

    private String url;
    private String format;

    /**
     * Factory method to return an htsget reader.
     */
    static HtsgetReader getReader(HtsgetUtils.Metadata metadata) {
        String format = metadata.getFormat();
        return new HtsgetReader(metadata.getUrl(), format);
    }

    HtsgetReader(String url, String format) {
        this.url = url;
        this.format = format.toUpperCase();
    }

    byte[] readHeader() throws IOException {
        URL headerURL = HtsgetUtils.addQueryString(this.url, "class=header&format=" + this.format);
        String ticketString = HttpUtils.getInstance().getContentsAsJSON(headerURL);
        JsonParser parser = new JsonParser();
        JsonObject ticket = parser.parse(ticketString).getAsJsonObject();

        // The UMCCR htsget server returns a bgzipped header for VCF format.  Arguably a server bug, but we

        byte [] bytes =  loadURLs(ticket);
        if (bytes != null && bytes.length >= 2 && bytes[0] == (byte) 0x1F && bytes[1] == (byte) 0x8B) {
            BlockCompressedInputStream bis = new BlockCompressedInputStream(new ByteArrayInputStream(bytes));
            return bis.readAllBytes();
        } else {
            return bytes;
        }
    }

    byte[] readData(String chr, int start, int end) throws IOException {
        final String queryString = "format=" + this.format + "&referenceName=" + chr + "&start=" + start + "&end=" + end;
        URL queryURL = HtsgetUtils.addQueryString(this.url, queryString);
        String ticketString = HttpUtils.getInstance().getContentsAsJSON(queryURL);
        JsonParser parser = new JsonParser();
        JsonObject ticket = parser.parse(ticketString).getAsJsonObject();
        return loadURLs(ticket);
    }

    private byte[] loadURLs(JsonObject ticket) throws IOException {

        JsonObject container = ticket.getAsJsonObject("htsget");
        JsonArray urls = container.get("urls").getAsJsonArray();
        Iterator<JsonElement> iter = urls.iterator();
        byte[] bytes = null;
        while (iter.hasNext()) {

            JsonObject next = iter.next().getAsJsonObject();
            String urlString = next.get("url").getAsString();
            if (urlString.startsWith("data:")) {
                throw new RuntimeException("Data URLs are not currently supported");
            }

            URL url = new URL(urlString);
            Map<String, String> headers = new HashMap<>();
            JsonObject headerObj = next.getAsJsonObject("headers");
            if (headerObj != null) {
                for (String key : headerObj.keySet()) {
                    headers.put(key, headerObj.get(key).getAsString());
                }
            }

            byte[] b = HttpUtils.getInstance().getContentsAsBytes(url, headers);
            if (bytes == null) {
                bytes = b;
            } else {
                byte[] newBytes = new byte[bytes.length + b.length];
                System.arraycopy(bytes, 0, newBytes, 0, bytes.length);
                System.arraycopy(b, 0, newBytes, bytes.length, b.length);
                bytes = newBytes;
            }
        }
        return bytes;
    }
}

