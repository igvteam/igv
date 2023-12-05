package org.broad.igv.ucsc;

import java.io.IOException;

public interface BPIndex {
    long [] search(String term) throws IOException;
}
