package org.igv.data.cufflinks;

import org.igv.feature.Range;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Parses various cufflinks (and cuffdiff) output files as described here:
 * <p/>
 * http://cufflinks.cbcb.umd.edu/manual.html
 *
 * @author jrobinso
 *         Date: 3/8/13
 *         Time: 2:30 PM
 */
public class CufflinksParser {

    public static List<? extends Range> parse(String path) throws IOException {

        final String s = path.toLowerCase();
        if (s.endsWith("fpkm_tracking")) {
            AsciiFeatureCodec<FPKMValue> codec = new FPKMTrackingCodec(path);
            return parse(codec, path);
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            AsciiFeatureCodec<ExpDiffValue> codec = new ExpDiffCodec(path);
            return parse(codec, path);
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }

    }

    public static <T extends Range> List<T> parse(AsciiFeatureCodec<T> codec, String path) throws IOException {
        List<T> values = new ArrayList<T>();
        AbstractFeatureReader reader = AbstractFeatureReader.getFeatureReader(path, codec, false);
        Iterator<T> iter = reader.iterator();
        while(iter.hasNext()){
            values.add(iter.next());
        }
        return values;
    }

}
