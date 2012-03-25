package org.broad.igv.feature.genome;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Implementation of Sequence backed by an indexed fasta file
 *
 * @author Jim Robinson
 * @date 3/24/12
 */

public class FastaDirectorySequence implements Sequence {

    Map<String, FastaSequence> sequenceMap;

    public FastaDirectorySequence(String directoryPath, String [] fastaFiles) throws IOException {
        readIndeces(directoryPath, fastaFiles);
    }

    private void readIndeces(String directoryPath, String [] fastaFiles) throws IOException {
        sequenceMap = new LinkedHashMap<String, FastaSequence>();
        for (String file : fastaFiles) {
            String fastaPath = directoryPath + "/" + file;
            FastaSequence fastaSequence = new FastaSequence(fastaPath);
            for (String chr : fastaSequence.getChromosomeNames()) {
                sequenceMap.put(chr, fastaSequence);
            }
        }
    }

    public Collection<FastaSequence> getFastaSequences() {
        return sequenceMap.values();
    }

    public byte[] readSequence(String chr, int start, int end) {


        if (!sequenceMap.containsKey(chr)) {
            return null;
        }
        return sequenceMap.get(chr).readSequence(chr, start, end);
    }


}
