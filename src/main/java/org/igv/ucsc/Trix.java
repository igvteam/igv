package org.igv.ucsc;
//
// This is a port of trix-js from the GMOD repository:  https://github.com/GMOD/trix-js,
// developed by Colin Diesh, Robert Buels, and Matt Morgan.

import htsjdk.samtools.seekablestream.SeekableStream;
import org.igv.Globals;
import org.igv.util.ParsingUtils;
import org.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Reference: https://genome.ucsc.edu/goldenpath/help/trix.html
 */
public class Trix {

    static int ADDRESS_SIZE = 10;

    String ixFile;  // URL to the ix file
    String ixxFile;  // URL to the ixx file
    List<IndexEntry> index;
    Map<Integer, String> bufferCache = new HashMap<>();

    public Trix(String ixxFile, String ixFile) {
        this.ixxFile = ixxFile;
        this.ixFile = ixFile;
    }

    public Map<String, String[]> search(String searchString) throws IOException {

        String[] searchWords = Globals.whitespacePattern.split(searchString);

        // we only support a single search term
        String searchWord = searchWords[0].toLowerCase();
        String str = this._getBuffer(searchWord);
        if (str == null) {
            return null;
        }

        String[] lines = str
                .substring(0, str.lastIndexOf('\n'))
                .split("\n");

        List<String> matches = new ArrayList<>();
        for (String line : lines) {
            if (line.trim().length() == 0) continue;
            String word = Globals.whitespacePattern.split(line)[0];
            boolean match = word.startsWith(searchWord);
            if (match) {
                matches.add(line);
            }
            // we are done scanning if we are lexicographically greater than the search string
            if (word.substring(0, Math.min(word.length(), searchWord.length())).compareTo(searchWord) > 0) {
                break;
            }
        }

        if (matches.size() == 0) {
            return null;
        } else {
            Map<String, String[]> results = new HashMap<>();
            for (String m : matches) {
                String[] parts = Globals.whitespacePattern.split(m);
                String term = parts[0];
                String[] t = new String[parts.length - 1];
                for (int i = 1; i < parts.length; i++) {
                    t[i - 1] = parts[i].split(",")[0];
                }
                results.put(term, t);
            }
            ;
            return results;
        }
    }


    List<IndexEntry> getIndex() throws IOException {
        if (this.index == null) {
            this.index = readIndex();
        }
        return this.readIndex();
    }

    List<IndexEntry> readIndex() throws IOException {
        List<IndexEntry> index = new ArrayList<>();
        try (BufferedReader br = ParsingUtils.openBufferedReader(ixxFile)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().length() == 0) continue;
                int p = line.length() - ADDRESS_SIZE;
                String prefix = line.substring(0, p);
                String posStr = line.substring(p);
                int pos = Integer.parseInt(posStr, 16);
                index.add(new IndexEntry(prefix, pos));
            }
        }
        return index;
    }

    String _getBuffer(String searchWord) throws IOException {

        int start = 0;
        int end = 65536;
        List<IndexEntry> indexes = this.getIndex();
        for (IndexEntry entry : indexes) {
            int trimEnd = Math.min(entry.key.length(), searchWord.length());
            String trimmedKey = entry.key.substring(0, trimEnd);
            if (trimmedKey.compareTo(searchWord) < 0) {
                start = entry.value;
                end = entry.value + 65536;
            } else {
                break;
            }
        }

        // Return the buffer and its end position in the file.
        int len = end - start;
        if (len < 0) {
            return null;
        }

        if (this.bufferCache.containsKey(start)) {
            return this.bufferCache.get(start);
        } else {
            SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(this.ixFile);
            byte[] bytes = new byte[len];
            is.seek(start);
            readFullyOrToEnd(bytes, is);
            String buffer = new String(bytes);
            this.bufferCache.put(start, buffer);
            return buffer;
        }
    }

    public void readFullyOrToEnd(byte b[], SeekableStream is) throws IOException {
        int len = b.length;
        if (len < 0) {
            throw new IndexOutOfBoundsException();
        }
        int n = 0;
        while (n < len) {
            int count = is.read(b, n, len - n);
            if (count < 0) {
                return;  // EOF
            }
            n += count;
        }
    }

    static class IndexEntry {
        String key;
        int value;

        public IndexEntry(String key, int value) {
            this.key = key;
            this.value = value;
        }
    }
}
