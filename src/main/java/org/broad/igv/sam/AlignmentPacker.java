/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import htsjdk.samtools.SAMTag;
import org.broad.igv.feature.Range;
import org.broad.igv.feature.Strand;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.Track;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.*;

/**
 * Packs alignments tightly with no overlap
 *
 * @author jrobinso
 */
public class AlignmentPacker {

    private static final Logger log = LogManager.getLogger(AlignmentPacker.class);

    /**
     * Minimum gap between the end of one alignment and start of another.
     */
    public static final int MIN_ALIGNMENT_SPACING = 2;

    /**
     * Maximum number of rows to search for gaps when tight packing.
     * This limits the worst-case complexity to O(n * MAX_GAP_SEARCH_DEPTH).
     * Set to a reasonable value that balances tight packing with performance.
     */
    private static final int MAX_GAP_SEARCH_DEPTH = 100;

    private static final String NULL_GROUP_VALUE = "";


    /**
     * Allocates each alignment to row such that there is no overlap.
     */
    public PackedAlignments packAlignments(
            AlignmentInterval interval,
            AlignmentTrack.RenderOptions renderOptions,
            ReferenceFrame referenceFrame,
            Track.DisplayMode displayMode) {

        LinkedHashMap<String, List<Row>> packedAlignments = new LinkedHashMap<>();

        List<Alignment> alList = interval.getAlignments();

        // TODO -- means to undo this
        if (renderOptions.isLinkedReads()) {
            alList = linkByTag(alList, renderOptions.getLinkByTag());
        }

        if (renderOptions.getGroupByOption() == AlignmentTrack.GroupOption.NONE) {
            List<Row> alignmentRows = new ArrayList<>(10000);
            if (displayMode == Track.DisplayMode.FULL) {
                packFull(alList, renderOptions, alignmentRows, referenceFrame);
            } else {
                packDense(alList, renderOptions, alignmentRows);
            }
            packedAlignments.put("", alignmentRows);
        } else {

            // Separate alignments into groups.
            Map<Object, List<Alignment>> groupedAlignments = new HashMap<>();
            for (final Alignment alignment : alList) {
                Object groupKey = getGroupValue(alignment, renderOptions);
                if (groupKey == null) {
                    groupKey = NULL_GROUP_VALUE;
                }
                List<Alignment> groupList = groupedAlignments.computeIfAbsent(groupKey, k -> new ArrayList<>(1000));
                groupList.add(alignment);
            }


            // Now sort the groups by their name and pack each group individually
            List<Object> keys = new ArrayList<Object>(groupedAlignments.keySet());
            Comparator<Object> groupComparator = getGroupComparator(renderOptions.getGroupByOption());

            // Certain group options sort descending by default, as indicated by the "reverse" property
            if (renderOptions.getGroupByOption().reverse) {
                groupComparator = groupComparator.reversed();
            }

            if (renderOptions.isInvertGroupSorting()) {
                groupComparator = groupComparator.reversed();
            }
            keys.sort(groupComparator);

            for (Object key : keys) {
                List<Row> alignmentRows = new ArrayList<>(10000);
                List<Alignment> group = groupedAlignments.get(key);

                if (displayMode == Track.DisplayMode.FULL) {
                    packFull(group, renderOptions, alignmentRows, referenceFrame);
                } else {
                    packDense(group, renderOptions, alignmentRows);
                }
                packedAlignments.put(key.toString(), alignmentRows);
            }
        }

        List<AlignmentInterval> tmp = new ArrayList<AlignmentInterval>();
        tmp.add(interval);
        return new PackedAlignments(tmp, packedAlignments);
    }


    private void packDense(List<Alignment> alList,
                           AlignmentTrack.RenderOptions renderOptions,
                           List<Row> alignmentRows) {

        if (alList == null || alList.isEmpty()) return;

        boolean isPairedAlignments = renderOptions.isViewPairs();

        // Process pairing
        Map<String, PairedAlignment> pairs = isPairedAlignments ? new HashMap<>(1000) : null;
        List<Alignment> alignmentsToPack = new ArrayList<>();

        for (Alignment al : alList) {
            if (al.isMapped()) {
                Alignment alignment = al;

                // Pair alignments -- do not pair secondary alignments
                if (isPairedAlignments && isPairable(al)) {
                    String readName = al.getReadName();
                    PairedAlignment pair = pairs.get(readName);
                    if (pair == null) {
                        pair = new PairedAlignment(al);
                        pairs.put(readName, pair);
                        alignment = pair;
                    } else {
                        // Add second alignment to pair.
                        pair.setSecondAlignment(al);
                        pairs.remove(readName);
                        continue;
                    }
                }

                alignmentsToPack.add(alignment);
            }
        }

        // Pack alignments into rows, ensuring tight packing with no overlaps.
        // Optimization: Use a two-phase approach to balance tight packing with performance.
        // Phase 1 (fast-path): Check if alignment fits at the end of any row (common case)
        // Phase 2 (gap-filling): Search a limited number of rows for gaps (bounded complexity)
        for (Alignment alignment : alignmentsToPack) {
            boolean placed = false;
            int alignmentStart = alignment.getStart();

            // Phase 1: Fast-path check for appending to row ends
            // This handles the common case (genomically ordered alignments) in O(1) per alignment
            for (Row row : alignmentRows) {
                if (alignmentStart >= row.getLastEnd() + MIN_ALIGNMENT_SPACING) {
                    row.addAlignment(alignment);
                    placed = true;
                    break;
                }
            }

            // Phase 2: Gap-filling check if alignment couldn't be appended
            // Limit search depth to prevent O(n*m) complexity in worst case
            if (!placed) {
                int searchDepth = Math.min(alignmentRows.size(), MAX_GAP_SEARCH_DEPTH);
                for (int i = 0; i < searchDepth; i++) {
                    Row row = alignmentRows.get(i);
                    if (canFitInRowGap(row, alignment)) {
                        row.addAlignment(alignment);
                        placed = true;
                        break;
                    }
                }
            }

            // If no existing row works, create a new row
            if (!placed) {
                Row newRow = new Row();
                newRow.addAlignment(alignment);
                alignmentRows.add(newRow);
            }
        }

        if (log.isDebugEnabled()) {
            log.debug("Packed " + alList.size() + " alignments into " + alignmentRows.size() + " rows");
        }
    }


    private void packFull(List<Alignment> alList,
                          AlignmentTrack.RenderOptions renderOptions,
                          List<Row> alignmentRows,
                          ReferenceFrame referenceFrame) {

        Map<String, PairedAlignment> pairs = null;

        boolean isPairedAlignments = renderOptions.isViewPairs();

        if (isPairedAlignments) {
            pairs = new HashMap<>(1000);
        }

        Range genomicRange = referenceFrame.getCurrentRange();

        for (Alignment al : alList) {

            if (al.isMapped()) {
                Alignment alignment = al;

                // Pair alignments -- do not pair secondaryalignments
                if (isPairedAlignments && isPairable(al)) {
                    String readName = al.getReadName();
                    PairedAlignment pair = pairs.get(readName);
                    if (pair == null) {
                        pair = new PairedAlignment(al);
                        pairs.put(readName, pair);
                        alignment = pair;
                    } else {
                        // Add second alignment to pair.
                        pair.setSecondAlignment(al);
                        pairs.remove(readName);
                        continue;
                    }
                }

                if (alignment.getEnd() > genomicRange.start && alignment.getStart() < genomicRange.end) {
                    Row row = new Row();
                    row.addAlignment(alignment);
                    alignmentRows.add(row);
                }
            }
        }
    }

    /**
     * Check if an alignment can fit in a gap within a row without overlapping existing alignments.
     * This method is only called when the alignment cannot be appended to the end of the row,
     * so it needs to check for gaps between existing alignments.
     *
     * Performance note: This is the slow path, only used when filling gaps in rows.
     * The common case (appending to row end) is handled by the fast-path check in packDense.
     */
    private boolean canFitInRowGap(Row row, Alignment alignment) {
        int alignmentStart = alignment.getStart();
        int alignmentEnd = alignment.getEnd();

        // Check against each alignment in the row to see if there's a gap
        List<Alignment> existingAlignments = row.getAlignments();
        for (Alignment existing : existingAlignments) {
            int existingStart = existing.getStart();
            int existingEnd = existing.getEnd();

            // Check for overlap with spacing
            if (!(alignmentEnd + MIN_ALIGNMENT_SPACING <= existingStart ||
                    alignmentStart >= existingEnd + MIN_ALIGNMENT_SPACING)) {
                return false;  // Overlaps with this alignment
            }
        }
        return true;
    }

    private boolean isPairable(Alignment al) {
        return al.isPrimary() &&
                al.isPaired() &&
                al.getMate().isMapped() &&
                al.getMate().getChr().equals(al.getChr());
    }

    private List<Alignment> linkByTag(List<Alignment> alList, String tag) {

        List<Alignment> bcList = new ArrayList<>(alList.size() / 10);
        Map<Object, LinkedAlignment> map = new HashMap<>(bcList.size() * 2);

        for (Alignment a : alList) {

            if (a.isPrimary()) {
                Object bc;
                if ("READNAME".equals(tag)) {
                    bc = a.getReadName();
                } else {
                    bc = a.getAttribute(tag);
                }

                if (bc == null) {
                    bcList.add(a);
                } else {
                    LinkedAlignment linkedAlignment = map.get(bc);
                    if (linkedAlignment == null) {
                        linkedAlignment = new LinkedAlignment(tag, bc.toString());
                        map.put(bc, linkedAlignment);
                        bcList.add(linkedAlignment);
                    }
                    linkedAlignment.addAlignment(a);
                }
            } else {
                // Don't link secondary (i.e alternative) alignments
                bcList.add(a);
            }
        }

        // Now copy list, de-linking orphaned alignments (alignments with no linked mates)
        List<Alignment> delinkedList = new ArrayList<>(alList.size());
        for (Alignment a : bcList) {
            if (a instanceof LinkedAlignment) {
                final List<Alignment> alignments = ((LinkedAlignment) a).alignments;
                if (alignments.size() == 1) {
                    delinkedList.add(alignments.get(0));
                } else {
                    a.finish();
                    delinkedList.add(a);
                }
            } else {
                delinkedList.add(a);
            }
        }

        return delinkedList;
    }


    /**
     * Gets the range over which alignmentsList spans. Asssumes all on same chr, and sorted
     *
     * @param alignmentsList
     * @return
     */
    private Range getAlignmentListRange(List<Alignment> alignmentsList) {
        if (alignmentsList == null || alignmentsList.size() == 0) return null;
        Alignment firstAlignment = alignmentsList.get(0);

        int minStart = firstAlignment.getStart();
        int maxEnd = firstAlignment.getEnd();
        for (Alignment alignment : alignmentsList) {
            maxEnd = Math.max(maxEnd, alignment.getEnd());
        }
        return new Range(firstAlignment.getChr(), minStart,
                maxEnd);
    }

    private Comparator<Object> getGroupComparator(AlignmentTrack.GroupOption groupByOption) {
        switch (groupByOption) {
            case PAIR_ORIENTATION:
                return new PairOrientationComparator();
            default:
                //Sort null values towards the end
                return new Comparator<Object>() {
                    @Override
                    public int compare(Object o1, Object o2) {
                        if (o1 == null && o2 == null) {
                            return 0;
                        } else if (o1 == null) {
                            return 1;
                        } else if (o2 == null) {
                            return -1;
                        } else {
                            // no nulls
                            if (o1.equals(o2)) {
                                return 0;
                            } else if (o1 instanceof String && NULL_GROUP_VALUE.equals(o1)) {
                                return 1;
                            } else if (o2 instanceof String && NULL_GROUP_VALUE.equals(o2)) {
                                return -1;
                            } else {
                                if (o1 instanceof Integer && o2 instanceof Integer) {
                                    Integer i1 = (Integer) o1, i2 = (Integer) o2;
                                    return i1.compareTo(i2);
                                } else if (o1 instanceof Float && o2 instanceof Float) {
                                    Float f1 = (Float) o1, f2 = (Float) o2;
                                    return f1.compareTo(f2);
                                } else if (o1 instanceof Double && o2 instanceof Double) {
                                    Double d1 = (Double) o1, d2 = (Double) o2;
                                    return d1.compareTo(d2);
                                } else {
                                    String s1 = o1.toString(), s2 = o2.toString();
                                    return s1.compareToIgnoreCase(s2);
                                }
                            }
                        }
                    }
                };
        }
    }

    private Object getGroupValue(Alignment al, AlignmentTrack.RenderOptions renderOptions) {

        AlignmentTrack.GroupOption groupBy = renderOptions.getGroupByOption();
        String tag = renderOptions.getGroupByTag();
        Range pos = renderOptions.getGroupByPos();
        String readNameParts[], movieName, zmw;

        return switch (groupBy) {
            case CLUSTER -> al.getClusterName();
            case STRAND -> al.isNegativeStrand() ? "-" : "+";
            case SAMPLE -> al.getSample();
            case LIBRARY -> al.getLibrary();
            case READ_GROUP -> al.getReadGroup();
            case LINKED -> (al instanceof LinkedAlignment) ? "Linked" : "";
            case PHASE -> al.getAttribute("HP");
            case TAG -> {
                Object tagValue = tag == null ? null : al.getAttribute(tag);
                if (tagValue == null) {
                    yield null;
                } else if (tagValue instanceof Integer || tagValue instanceof Float || tagValue instanceof Double) {
                    yield tagValue;
                } else {
                    yield tagValue.toString();
                }
            }
            case FIRST_OF_PAIR_STRAND -> {
                Strand strand = al.getFirstOfPairStrand();
                yield strand == Strand.NONE ? null : strand.toString();
            }
            case READ_ORDER -> {
                if (al.isPaired() && al.isFirstOfPair()) {
                    yield "FIRST";
                } else if (al.isPaired() && al.isSecondOfPair()) {
                    yield "SECOND";
                } else {
                    yield "";
                }
            }
            case PAIR_ORIENTATION -> {
                PEStats peStats = AlignmentRenderer.getPEStats(al, renderOptions);
                AlignmentTrack.OrientationType type = AlignmentRenderer.getOrientationType(al, peStats);
                if (type == null) {
                    yield AlignmentTrack.OrientationType.UNKNOWN.name();
                }
                yield type.name();
            }
            case MATE_CHROMOSOME -> {
                ReadMate mate = al.getMate();
                if (mate == null) {
                    yield null;
                }
                if (!mate.isMapped()) {
                    yield "UNMAPPED";
                } else {
                    yield mate.getChr();
                }
            }
            case NONE -> null;
            case CHIMERIC -> al.getAttribute(SAMTag.SA.name()) != null ? "CHIMERIC" : "";
            case SUPPLEMENTARY -> al.isSupplementary() ? "SUPPLEMENTARY" : "";
            case REFERENCE_CONCORDANCE -> !al.isProperPair() ||
                    al.getCigarString().toUpperCase().contains("S") ||
                    al.isSupplementary() ?
                    "DISCORDANT" : "";
            case BASE_AT_POS -> {
                // Use a string prefix to enforce grouping rules:
                //    1: alignments with a base at the position
                //    2: alignments with a gap at the position
                //    3: alignment that do not overlap the position (or are on a different chromosome)

                if (pos != null &&
                        al.getChr().equals(pos.getChr()) &&
                        al.getAlignmentStart() <= pos.getStart() &&
                        al.getAlignmentEnd() > pos.getStart()) {

                    byte[] baseAtPos = new byte[]{al.getBase(pos.getStart())};
                    if (baseAtPos[0] == 0) { // gap at position
                        yield "2:";
                    } else { // base at position
                        yield "1:" + new String(baseAtPos);
                    }
                } else { // does not overlap position
                    yield "3:";
                }
            }
            case INSERTION_AT_POS -> {
                // Use a string prefix to enforce grouping rules:
                //    1: alignments with a base at the position
                //    2: alignments with a gap at the position
                //    3: alignment that do not overlap the position (or are on a different chromosome)
                if (pos != null &&
                        al.getChr().equals(pos.getChr()) &&
                        al.getAlignmentStart() <= pos.getStart() &&
                        al.getAlignmentEnd() > pos.getStart()) {
                    int insertionBaseCount = 0;
                    AlignmentBlock leftInsertion = al.getInsertionAt(pos.getStart() + 1);
                    if (leftInsertion != null) {
                        insertionBaseCount += leftInsertion.getLength();
                    }
                    AlignmentBlock rightInsertion = al.getInsertionAt(pos.getStart());
                    if (rightInsertion != null) {
                        insertionBaseCount += rightInsertion.getLength();
                    }
                    yield insertionBaseCount;

                } else {
                    yield 0;
                }
            }
            case MOVIE -> {
                readNameParts = al.getReadName().split("/");
                if (readNameParts.length < 3) {
                    yield "";
                }
                movieName = readNameParts[0];
                yield movieName; // group PacBio reads by movie
            }
            case ZMW -> {
                readNameParts = al.getReadName().split("/");
                if (readNameParts.length < 3) {
                    yield "";
                }
                movieName = readNameParts[0];
                zmw = readNameParts[1];
                yield movieName + "/" + zmw; // group PacBio reads by ZMW
            }
            case MAPPING_QUALITY -> al.getMappingQuality();
            case DUPLICATE -> al.isDuplicate() ? "duplicate" : "non-duplicate";
            case SELECTED -> renderOptions.getSelectedReadNames().containsKey(al.getReadName()) ? "SELECTED" : "";
        };
    }


    private static class PairOrientationComparator implements Comparator<Object> {
        private final List<AlignmentTrack.OrientationType> orientationTypes;
        //private final Set<String> orientationNames = new HashSet<String>(AlignmentTrack.OrientationType.values().length);

        public PairOrientationComparator() {
            orientationTypes = Arrays.asList(AlignmentTrack.OrientationType.values());
//            for(AlignmentTrack.OrientationType type: orientationTypes){
//                orientationNames.add(type.name());
//            }
        }

        @Override
        public int compare(Object o0, Object o1) {
            String s0 = o0.toString();
            String s1 = o1.toString();
            if (s0 != null && s1 != null) {
                AlignmentTrack.OrientationType t0 = AlignmentTrack.OrientationType.valueOf(s0);
                AlignmentTrack.OrientationType t1 = AlignmentTrack.OrientationType.valueOf(s1);
                return orientationTypes.indexOf(t0) - orientationTypes.indexOf(t1);
            } else if (s0 == null ^ s1 == null) {
                //exactly one is null
                return s0 == null ? 1 : -1;
            } else {
                //both null
                return 0;
            }

        }
    }

}
