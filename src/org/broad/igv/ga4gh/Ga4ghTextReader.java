package org.broad.igv.ga4gh;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.Ga4ghAlignment;
import org.broad.igv.sam.reader.AlignmentReader;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Class for testing GlobalAlliance API.  Loads json from a text file, for development only
 * <p/>
 * Created by jrobinso on 7/18/14.
 */

public class Ga4ghTextReader implements AlignmentReader<Alignment> {

    String path;
    List<Alignment> alignmentList;
    List<String> sequenceNames;

    public Ga4ghTextReader(String path) {
        this.path = path;

        try {
            loadAll();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public List<String> getSequenceNames() {
        return sequenceNames;
    }

    @Override
    public SAMFileHeader getFileHeader() {
        return null;
    }

    @Override
    public Set<String> getPlatforms() {
        return null;
    }

    @Override
    public CloseableIterator<Alignment> iterator() {
        return new MIterator();
    }

    @Override
    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) throws IOException {
        return new MIterator();
    }

    @Override
    public boolean hasIndex() {
        return true;
    }

    private void loadAll() throws IOException {


        alignmentList = new ArrayList<Alignment>();
       HashSet<String> seqNames = new HashSet<String>();

        BufferedReader br;
        StringBuffer sb = new StringBuffer();
        String line;

        br = new BufferedReader(new FileReader(path));

        while ((line = br.readLine()) != null) {
            sb.append(line + "\n");
        }

        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(sb.toString()).getAsJsonObject();

        JsonArray reads = obj.getAsJsonArray("reads");

        Iterator<JsonElement> iter = reads.iterator();
        while (iter.hasNext()) {
            JsonElement next = iter.next();
            Ga4ghAlignment alignment = new Ga4ghAlignment(next.getAsJsonObject());
            seqNames.add(alignment.getChr());
            alignmentList.add(alignment);
        }

        this.sequenceNames = new ArrayList(seqNames);
    }

    public static boolean supportsFileType(String path) {
        return path.contains("ga4gh");
    }

    class MIterator implements CloseableIterator<Alignment> {

        Iterator<Alignment> iter;

        MIterator() {
            iter = alignmentList.iterator();
        }

        @Override
        public void close() {

        }

        @Override
        public boolean hasNext() {
            return iter.hasNext();
        }

        @Override
        public Alignment next() {
            return iter.next();
        }

        @Override
        public void remove() {
            iter.remove();
        }
    }

}
