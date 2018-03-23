/*
 * The MIT License (MIT)
 *  Copyright (c) 2007-2015 by Institute for Computational Biomedicine,
 *                                          Weill Medical College of Cornell University.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.goby;

import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.lang.MutableString;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Query reader to parse <a href="http://goby.campagnelab.org">Goby</a> alignment files.
 * The compact alignment files must be sorted and indexed.
 * <p/>
 * For further information about Goby, or to obtain sample alignment files, see http://goby.campagnelab.org
 *
 * @author Fabien Campagne
 *         Date: Jun 29, 2010
 *         Time: 11:43:18 AM
 */
public class GobyAlignmentQueryReader implements AlignmentReader {
    private static final Logger LOG = Logger.getLogger(GobyAlignmentQueryReader.class);

    private AlignmentReaderImpl reader = null;
    private final String basename;
    private DoubleIndexedIdentifier targetIdentifiers;
    private boolean isIndexed;
    private List<String> targetSequenceNames;
    private static final CloseableIterator<Alignment> EMPTY_ITERATOR = new CloseableIterator<Alignment>() {
        public void close() {

        }

        public boolean hasNext() {
            return false;
        }

        public Alignment next() {
            return null;
        }

        public void remove() {

        }
    };


    /**
     * Construct a query reader for this filename/basename.
     *
     * @param filename The filename of any file component for a Goby compact alignment, or the alignment basename.
     * @throws IOException If an error occurs reading the alignment.
     */
    public GobyAlignmentQueryReader(String filename) throws IOException {
        basename = filename;
        reader = new AlignmentReaderImpl(filename);
        reader.readHeader();
        if (!reader.isIndexed()) {
            final String errorMessage = "Goby alignment files must be sorted in order to be loaded in IGV. See the IGV tutorial at http://goby.campagnelab.org/ for details.";
            System.err.println(errorMessage);
            throw new UnsupportedOperationException(errorMessage);
        }
        final IndexedIdentifier identifiers = reader.getTargetIdentifiers();
        // add MT as a synonym for M:
        identifiers.put(new MutableString("M"), identifiers.getInt(new MutableString("MT")));
        targetIdentifiers = new DoubleIndexedIdentifier(identifiers);
        isIndexed = reader.isIndexed();
        // reader.close();
        // reader = null;

        targetSequenceNames = new ArrayList();
        for (MutableString ms : identifiers.keySet()) {
            targetSequenceNames.add(ms.toString());
        }


    }


    /**
     * Release resources associated with this query reader. Closes the underlying Goby AlignmentReader.
     *
     * @throws IOException
     */
    public void close() throws IOException {
        if (reader != null) {

            reader.close();
            reader = null;
        }
    }

    public List<String> getSequenceNames() {
        return targetSequenceNames;
    }

    public SAMFileHeader getFileHeader(){
        return null;
    }

    /**
     * @return the list of platforms represented in this alignment? Please define this exactly.
     */
    public Set<String> getPlatforms() {
        return Collections.EMPTY_SET;
    }


    /**
     * Obtain an iterator over the entire file.
     *
     * @return An alignment iterator.
     */
    public CloseableIterator<Alignment> iterator() {
        return new GobyAlignmentIterator(getNewLocalReader(), targetIdentifiers);
    }

    /**
     * Obtain an iterator over a genomic window. The window on reference 'sequence' extends from
     * start to end. The attribute 'contained' is ignored (documentation required).
     *
     * @return An alignment iterator restricted to the sequence [start end] interval.
     */
    public final CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) {
        LOG.debug(String.format("query %s %d %d %b%n", sequence, start, end, contained));

        final MutableString id = new MutableString(sequence);
        int referenceIndex = targetIdentifiers.getIndex(id);
        if (referenceIndex == -1) {
            // try again removing a chr prefix:
            referenceIndex = targetIdentifiers.getIndex(id.replace("chr", ""));
        }
        try {
            if (referenceIndex == -1) {
                // the reference could not be found in the Goby alignment, probably a wrong reference choice. Not sure how
                // to inform the end user, but we send no results here:
                return EMPTY_ITERATOR;
            }
            return new GobyAlignmentIterator(getNewLocalReader(), targetIdentifiers, referenceIndex, sequence, start, end);
        } catch (IOException e) {
            LOG.error(e);
            return null;

        }

    }

    /**
     * Determines whether the file is indexed.
     *
     * @return True. Goby files must be indexed.
     */
    public final boolean hasIndex() {
        return isIndexed;
    }

    /**
     * Determines whether filename can be loaded by this QueryReader.
     * <p/>
     *
     * @param filename Name of a file component or alignment basename.
     * @return True if this implementation can load the alignment corresponding to this filename.
     */
    public static boolean supportsFileType(String filename) {
        if(!(filename.endsWith(".entries") || filename.endsWith(".header") || filename.endsWith(".index"))) {
            return false;
        }
        final boolean result = AlignmentReaderImpl.canRead(filename);
        LOG.debug(String.format("supportsFileType %s result=%b", filename, result));
        return result;
    }

    /**
     * get an unrestricted alignment reader.
     */
    private AlignmentReaderImpl getNewLocalReader() {
        if (reader != null) return reader;
        try {

            return new AlignmentReaderImpl(basename);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * get an alignment reader restricted to a window.
     */
    private AlignmentReaderImpl getNewLocalReader(int referenceIndex, int start, int end) {
        if (reader != null) return reader;
        try {

            return new AlignmentReaderImpl(basename, referenceIndex, start, referenceIndex, end);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

}
