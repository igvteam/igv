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
import edu.cornell.med.icb.goby.alignments.Alignments;
import edu.cornell.med.icb.goby.exception.GobyRuntimeException;
import edu.cornell.med.icb.identifier.DoubleIndexedIdentifier;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.lang.MutableString;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;

import java.io.IOException;
import java.util.NoSuchElementException;

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
    private final AlignmentReaderImpl reader;
    protected final DoubleIndexedIdentifier indexToReferenceId;
    private int previousPosition = Integer.MIN_VALUE;
    private int previousReferenceIndex = -1;
    private boolean useWindow;
    // We need to record the sequence (chr), its redundant with the index but I'm not sure how to do the reverse lookup (JTR)
    private String reference;
    private int currentPosition;

    /**
     * Construct an iterator over the entire alignment file.
     *
     * @param reader            Goby alignment reader.
     * @param targetIdentifiers Bidirectional map from target index to target identifier/chromosome names.
     */
    public GobyAlignmentIterator(AlignmentReaderImpl reader, DoubleIndexedIdentifier targetIdentifiers) {
        this.reader = reader;
        this.indexToReferenceId = targetIdentifiers;
        startReferencePosition = 0;
        endReferencePosition=Integer.MAX_VALUE;
        previousPosition = Integer.MIN_VALUE;
        useWindow = false;

    }

    Int2ObjectMap<ObjectArrayList<GobyAlignment>> spliceCache = new Int2ObjectOpenHashMap<ObjectArrayList<GobyAlignment>>();

    /**
     * Cache alignments that we will need to refer to later. This is useful to keep alignments that are linked
     * across genomic positions by splice sites. The iterator only keeps alignments when configured with a window.
     * This means that alignments are not cached when calling the iterator() method to get all entries.
     *
     * @param alignment The alignment to be cached.
     * @return The list of gobyAlignments associated that have the same queryIndex as the cached alignment.
     */
    public ObjectArrayList<GobyAlignment> cacheSpliceComponent(GobyAlignment alignment) {
        if (useWindow) {
            final int queryIndex = alignment.entry.getQueryIndex();
            ObjectArrayList<GobyAlignment> list = spliceCache.get(queryIndex);
            if (list == null) {
                list = new ObjectArrayList<GobyAlignment>();
                spliceCache.put(queryIndex, list);
            }

            list.add(alignment);
            return list;
        }
        return null;

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
     * A constructor useful only for testing. Do not use for production.
     *
     * @param targetIndex Index of the reference sequence over which the window is defined.
     * @param start       Start position
     * @param end         End position
     * @throws IOException
     */
    protected GobyAlignmentIterator(int targetIndex, int start, int end) throws IOException {
        this.useWindow = true;
        this.targetIndex = targetIndex;
        this.startReferencePosition = start;
        this.endReferencePosition = end;
        this.reader = null;
        this.indexToReferenceId = null;
    }

    /**
     * Release resources used by this iterator. Do not close anything because we reused the Goby AlignmentReader (since Goby 1.9.6).
     */
    public void close() {

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
                if (!reader.hasNext()) return false;
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
            //  throw new RuntimeException("IO error reading next Goby alignment entry", e);
            return false;
        } catch (GobyRuntimeException e) {
            nextEntry = null;
            LOG.error(e);
            //  throw new RuntimeException("IO error reading next Goby alignment entry", e);
            return false;
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
        currentPosition = entry.getPosition();
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

    public MutableString getId(int targetIndex) {
        return indexToReferenceId.getId(targetIndex);
    }
}
