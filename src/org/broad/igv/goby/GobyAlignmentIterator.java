/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

import org.broad.igv.sam.Alignment;
import org.broad.igv.goby.GobyAlignment;
import org.apache.log4j.Logger;
import net.sf.samtools.util.CloseableIterator;
import edu.cornell.med.icb.goby.alignments.AlignmentReader;
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.alignments.AlignmentReaderImpl;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;

import java.util.NoSuchElementException;
import java.io.IOException;

/**
 * An iterator over <a href="http://goby.campagnelab.org">Goby</a> alignment entries.
 * IGV iterator return entries that fall within a window
 * of genomic location. Implement this logic, leveraging sorted/indexed
 * <a href="http://goby.campagnelab.org">Goby</a> alignments for speed.
 * <p/>
 * For further information about Goby, or to obtain sample alignment files, see http://goby.campagnelab.org
 *
 * @author Fabien Campagne
 *         Date: Jun 29, 2010
 *         Time: 12:01:40 PM
 */
public class GobyAlignmentIterator implements CloseableIterator<Alignment> {
    private static final Logger LOG = Logger.getLogger(GobyAlignmentIterator.class);

    private int targetIndex;
    private int startReferencePosition;
    private int endReferencePosition;
    private final AlignmentReader reader;
    protected final DoubleIndexedIdentifier indexToReferenceId;
    private int previousPosition = Integer.MIN_VALUE;
    private int previousReferenceIndex = -1;
    private boolean useWindow;
    // We need to record the sequence (chr), its redundant with the index but I'm not sure how to do the reverse lookup (JTR)
    private String reference;

    /**
     * Construct an iterator over the entire alignment file.
     *
     * @param reader            Goby alignment reader.
     * @param targetIdentifiers Bidirectional map from target index to target identifier/chromosome names.
     */
    public GobyAlignmentIterator(AlignmentReaderImpl reader, DoubleIndexedIdentifier targetIdentifiers) {
        this.reader = reader;
        this.indexToReferenceId = targetIdentifiers;
        startReferencePosition = Integer.MAX_VALUE;
        previousPosition = Integer.MIN_VALUE;
        useWindow = false;

    }

    /**
     * Construct an iterator over a window of the alignment file.
     *
     * @param reader            Goby alignment reader.
     * @param targetIdentifiers Bidirectional map from target index to target identifier/chromosome names.
     * @param referenceIndex    Index of the reference sequence/chromosome.
     * @param start             Minimum genomic location on the reference sequence that alignment entries must have to be returned.
     * @param end               Maximum genomic location on the reference sequence that alignment entries must have to be returned.
     * @throws java.io.IOException If an error occurs reading the Goby alignment.
     */
    public GobyAlignmentIterator(final AlignmentReaderImpl reader, final DoubleIndexedIdentifier targetIdentifiers,
                                 int referenceIndex, String chr, int start, int end) throws IOException {
        this(reader, targetIdentifiers);
        this.useWindow = true;
        if (referenceIndex != -1) {

            this.targetIndex = referenceIndex;
            this.reference = chr;
            this.startReferencePosition = start;
            this.endReferencePosition = end;
            this.previousReferenceIndex = referenceIndex;
            //  LOG.debug(String.format("reposition %d %d", referenceIndex, start));

            reader.reposition(referenceIndex, start);
        }

    }

    /**
     * Release resources used by this iterator. Closes the underlying Goby AlignmentReader.
     */
    public void close() {
        reader.close();
        //  LOG.info("closing " + this);
    }

    /**
     * Will hold the next entry, if it matches the window location criteria:
     */
    private Alignments.AlignmentEntry nextEntry = null;

    /**
     * Determine if this iterator has more alignment entries in the given window.
     *
     * @return True if next() will return an alignment, False otherwise.
     */
    public boolean hasNext() {
        /*    LOG.debug(String.format("previousPosition: %d endReferencePosition %d previousReferenceIndex %d targetIndex %d",
               previousPosition, endReferencePosition, previousReferenceIndex, targetIndex)
       ); */
        // Fetch the next entry with skipTo
        if (nextEntry != null) return true;

        try {
            if (!useWindow) {
                // all results are returned
                nextEntry = reader.next();
            } else {
                // we return only within a window
                nextEntry = reader.skipTo(targetIndex, startReferencePosition);

                if (nextEntry == null ||
                        (nextEntry.getTargetIndex() != targetIndex ||
                                nextEntry.getPosition() < startReferencePosition ||
                                nextEntry.getPosition() > endReferencePosition)) {
                    // No next entry, on a different target sequence, or before the position of interest:
                    nextEntry = null;
                }
            }
        } catch (IOException e) {
            nextEntry = null;
            LOG.error(e);
            throw new RuntimeException("IO error reading next Goby alignment entry", e);
        }

        final boolean result = nextEntry != null;
        // LOG.debug("hasNext returning :" + result);
        return result;
    }

    /**
     * Return the next alignment.
     *
     * @return the next alignment within the window.
     */
    public Alignment next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        Alignments.AlignmentEntry entry;

        entry = nextEntry;

        nextEntry = null;
        //    LOG.debug(String.format("next targetIndex: %d position: %d", targetIndex, entry.getPosition()));

        return new GobyAlignment(this, entry);
    }

    /**
     * This operation is not supported.
     */
    public void remove() {
        throw new UnsupportedOperationException();

    }

    /**
     * Return the reference sequence for this iterator
     *
     * @return The reference sequence
     */
    String getReference() {
        return reference;
    }

}
