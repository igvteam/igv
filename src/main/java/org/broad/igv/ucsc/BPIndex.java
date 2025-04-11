package org.broad.igv.ucsc;

import java.io.IOException;

public interface BPIndex {
    long searchForOffset(String term) throws IOException;
}
