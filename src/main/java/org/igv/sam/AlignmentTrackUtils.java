package org.igv.sam;

import org.igv.feature.Range;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.IGV;

import java.util.Comparator;
import java.util.List;


/**
 * Utility methods for alignment tracks (sorting, grouping, etc).  These methods perform sorting, grouping, and
 * coloring operations on all alignment tracks in the current IGV session.
 *
 */
public class AlignmentTrackUtils {

    private static final Logger log = LogManager.getLogger(AlignmentTrackUtils.class);

    /**
     * Sort all alignment tracks according to the specified option
     *
     * @param option
     * @api
     */
    public static void sortAlignmentTracks(SortOption option, String tag, final boolean invertSort) {
       List<AlignmentTrack> alignmentTracks = IGV.getInstance().getAlignmentTracks();
        for (AlignmentTrack track : alignmentTracks) {
            track.sortRows(option, tag, invertSort);
        }
        IGV.getInstance().repaint(alignmentTracks);
    }

    /**
     * Group all alignment tracks by the specified option.
     *
     * @param option
     * @api
     */
    public static void groupAlignmentTracks(AlignmentTrack.GroupOption option, String tag, Range pos) {

        List<AlignmentTrack> alignmentTracks = IGV.getInstance().getAlignmentTracks();
        for (AlignmentTrack t : alignmentTracks) {
            t.groupAlignments(option, tag, pos);
        }
        IGV.getInstance().repaint(alignmentTracks);
    }

    /**
     * Color all alignment tracks by the specified option.
     *
     * @param option
     * @api
     */
    public static void colorAlignmentTracks(AlignmentTrack.ColorOption option, String tag) {

        List<AlignmentTrack> alignmentTracks = IGV.getInstance().getAlignmentTracks();
        for (AlignmentTrack t : alignmentTracks) {
            final AlignmentTrack alignmentTrack = t;
            alignmentTrack.setColorOption(option);
            if (option == AlignmentTrack.ColorOption.BISULFITE && tag != null) {
                try {
                    AlignmentTrack.BisulfiteContext context = AlignmentTrack.BisulfiteContext.valueOf(tag);
                    alignmentTrack.setBisulfiteContext(context);
                } catch (IllegalArgumentException e) {
                    log.error("Error setting bisulfite context for: " + tag, e);
                }
            } else if (tag != null) {
                alignmentTrack.setColorByTag(tag);
            }
        }
        IGV.getInstance().repaint(alignmentTracks);
    }

    public static void packAlignmentTracks() {
        for (AlignmentTrack t : IGV.getInstance().getAlignmentTracks()) {
            t.packAlignments();
        }
    }

    public static void sortAlignments(SortOption option,
                                      String chr,
                                      double location,
                                      String tag,
                                      boolean invertSort,
                                      List<Alignment> alignments) {

        final int center = (int) location;
        byte referenceBase = getReference(chr, center);
        Comparator<Alignment> alignmentComparator = option.getAlignmentComparator(center, tag, referenceBase);
        if (invertSort) alignmentComparator = alignmentComparator.reversed();
        alignments.sort(alignmentComparator);
    }

    private static byte getReference(String chr, int pos) {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome == null) {
            return 0;
        }
        return genome.getReference(chr, pos);
    }
}
