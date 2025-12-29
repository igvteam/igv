package org.igv.feature.tribble;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.*;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.ParsingUtils;
import org.igv.util.ResourceLocator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 5/8/13
 *         Time: 11:02 AM
 */
public class TribbleReaderWrapper  implements IGVFeatureReader {


    FeatureReader<Feature> wrappedReader;

    public TribbleReaderWrapper(FeatureReader<Feature> wrappedReader) {
        this.wrappedReader = wrappedReader;
    }

    @Override
    public synchronized Iterator<Feature> query(String chr, int start, int end) throws IOException {

        // Tribble iterators must be closed, so we need to copy the features and insure closure before exiting.
        try (CloseableTribbleIterator<Feature> iter = wrappedReader.query(chr, start + 1, end)) {
            List<Feature> featureList = new ArrayList<Feature>();
            while (iter.hasNext()) {
                Feature f = iter.next();
                if (f.getStart() > end) {
                    break;
                } else if (f.getEnd() < start) {
                    continue;
                } else {
                    featureList.add(f);
                }
            }
            return featureList.iterator();
        }
    }

    @Override
    public CloseableIterator<Feature> iterator() throws IOException {
        // Note: Technically this is a file handle leak as the "close" method of the tribble iterator is not called.
        // In practice this is not a problem as the iterator() method is only called by batch programs transversing
        // the entire file.   It is none-the-less a file handle leak that should be addressed at some point.
        return wrappedReader.iterator();
    }

    @Override
    public List<String> getSequenceNames() {
        return wrappedReader.getSequenceNames();
    }

    @Override
    public Object getHeader() {
        return wrappedReader.getHeader();
    }


}
