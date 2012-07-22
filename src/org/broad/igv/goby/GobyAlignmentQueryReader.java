/*
 * Copyright (c) 2007-2011 by Institute for Computational Biomedicine,
 *                                          Weill Medical College of Cornell University.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.goby;

import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import edu.cornell.med.icb.identifier.IndexedIdentifier;
import it.unimi.dsi.lang.MutableString;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentReader;

import java.io.IOException;
import java.util.*;

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
