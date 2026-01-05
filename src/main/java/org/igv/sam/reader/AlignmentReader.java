
package org.igv.sam.reader;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.igv.sam.Alignment;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author jrobinso
 */
public interface AlignmentReader<T extends Alignment> {

    void close() throws IOException;

    /**
     * Return the list of sequence (chromosome) names as defined in the files header or meta-data section.
     */
    List<String> getSequenceNames() throws IOException;

    default Map<String, Long> getSequenceDictionary() {
        return null;
    }

    /**
     * Return the header of the SAM file. May be null
     * @return
     */
    SAMFileHeader getFileHeader();

    /**
     * Return the set of all platforms represented in this file.
     * May return "null"
     */
    Set<String>  getPlatforms();

    CloseableIterator<T> iterator() throws IOException;

    /**
     * Query alignments over a given range. Be careful about start/end,
     * SAMTools uses 1-based, but IGV uses 0-based.
     * This function requires hasIndex() == true.
     *
     *
     * @param sequence
     * @param start 0-based start location
     * @param end 0-based, exclusive-end coordinate
     * @param contained
     * @return
     * @throws IOException
     */
    CloseableIterator<T> query(final String sequence, final int start, final int end, final boolean contained) throws IOException;

    default void cancelQuery() {};

    boolean hasIndex();

}
