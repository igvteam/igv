package org.broad.igv.feature.genome;

import java.io.IOException;
import java.util.*;

/**
 * Implementation of Sequence backed by an indexed fasta file
 *
 * @author Jim Robinson
 * @date 3/24/12
 */

public class FastaDirectorySequence implements Sequence {

    Map<String, FastaIndexedSequence> sequenceMap;

    public FastaDirectorySequence(String directoryPath, String [] fastaFiles) throws IOException {
        readIndeces(directoryPath, fastaFiles);
    }

    private void readIndeces(String directoryPath, String [] fastaFiles) throws IOException {
        sequenceMap = new LinkedHashMap<String, FastaIndexedSequence>();
        for (String file : fastaFiles) {
            String fastaPath = directoryPath + "/" + file;
            FastaIndexedSequence fastaSequence = new FastaIndexedSequence(fastaPath);
            for (String chr : fastaSequence.getChromosomeNames()) {
                sequenceMap.put(chr, fastaSequence);
            }
        }
    }

    public Collection<FastaIndexedSequence> getFastaSequences() {
        return sequenceMap.values();
    }

    public byte[] readSequence(String chr, int start, int end) {


        if (!sequenceMap.containsKey(chr)) {
            return null;
        }
        return sequenceMap.get(chr).readSequence(chr, start, end);
    }


}
