package org.igv.circview;

import htsjdk.samtools.SAMTag;
import htsjdk.tribble.Feature;
import org.igv.bedpe.BedPE;
import org.igv.circview.model.Assembly;
import org.igv.circview.model.Chord;
import org.igv.circview.model.Chromosome;
import org.igv.circview.model.Mate;
import org.igv.circview.ui.CircularView;
import org.igv.circview.ui.CircularViewConfig;
import org.igv.circview.ui.CircularViewPanel;
import org.igv.circview.util.ChrColors;
import org.igv.circview.util.ColorUtils;
import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.sam.Alignment;
import org.igv.ui.IGV;
import org.igv.util.Downsampler;
import org.igv.variant.Variant;
import org.igv.variant.vcf.MateVariant;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Bridge between IGV and the in-process {@link CircularView} widget.
 *
 * <p>This replaces the JBrowse/Electron based {@code org.igv.jbrowse.CircularViewUtilities},
 * which spoke to an external app over a socket. The view is now a singleton Swing
 * {@link JFrame} owned by this class; tracks add chords through the static
 * {@code add*} methods, and the window is created and shown lazily on the first add.
 *
 * <p>Chord features are built from IGV model objects via the {@code Chord.from*}
 * factory methods (the same conversions the old code used), then downsampled to
 * {@link #MAX_CHORDS} to keep rendering responsive.
 */
public class CircularViewUtilities {

    /**
     * Maximum number of chords to render in a single set.
     */
    static int MAX_CHORDS = 10000;

    /**
     * Flanking bases added on each side when navigating to a clicked chord's regions.
     */
    private static final int CLICK_FLANKING = 2000;

    // Singleton view + window, lazily created on the EDT.
    private static CircularView view;
    private static CircularViewPanel panel;
    private static JFrame frame;

    /**
     * Id of the genome currently set as the assembly, to avoid needless resets.
     */
    private static String currentGenomeId;

    private CircularViewUtilities() {
    }

    // ---- Public API ---------------------------------------------------------

    /**
     * True if the circular view window currently exists and is showing.
     */
    public static boolean isOpen() {
        return frame != null && frame.isVisible();
    }

    public static void addBedPE(List<? extends BedPE> features, String trackName, Color color) {
        List<Chord> chords = new ArrayList<>(features.size());
        for (BedPE f : features) {
            chords.add(Chord.fromBedPE(f));
        }
        addChords(chords, trackName, color, 0.5f);
    }

    public static void addAlignments(List<Alignment> alignments, String trackName, Color color) {
        List<Chord> chords = new ArrayList<>();
        for (Alignment a : alignments) {
            if (a.isPaired() && a.getMate().isMapped()) {
                chords.add(Chord.fromPEAlignment(a));
            }
            if (a.getAttribute(SAMTag.SA.name()) != null) {
                chords.addAll(Chord.fromSAString(a));
            }
        }
        addChords(chords, trackName, color, 0.1f);
    }

    public static void addVariants(List<Feature> variants, String trackName, Color color) {
        List<Chord> chords = new ArrayList<>(variants.size());
        for (Feature f : variants) {
            if (f instanceof Variant) {
                Variant v = f instanceof MateVariant ? ((MateVariant) f).mate : (Variant) f;
                chords.add(Chord.fromVariant(v));
            }
        }
        addChords(chords, trackName, color, 0.5f);
    }

    /**
     * Add a set of chords to the view, opening (and creating) the window if needed.
     *
     * @param chords    the chord features
     * @param trackName name for this chord set / track row
     * @param color     base color; {@code alpha} is applied to it
     * @param alpha     opacity fraction in [0, 1]
     */
    public static void addChords(List<Chord> chords, String trackName, Color color, float alpha) {
        Chord[] arr = chords.toArray(new Chord[0]);
        if (arr.length > MAX_CHORDS) {
            arr = new Downsampler<Chord>().sample(arr, MAX_CHORDS);
        }
        final List<Chord> sampled = new ArrayList<>(Arrays.asList(arr));
        final Color c = ColorUtils.setAlpha(color, alpha);
        runOnEdt(() -> {
            CircularView v = getInstance();
            ensureAssembly(v);
            v.addChords(sampled, trackName, c);
            open();
        });
    }

    /**
     * Update the assembly to the given genome. No-op unless the view already
     * exists; on a genome switch the previously added chords no longer apply, so
     * this resets the view (matching {@link CircularView#setAssembly}).
     */
    public static void changeGenome(Genome genome) {
        if (view == null || genome == null) {
            return;
        }
        runOnEdt(() -> {
            view.setAssembly(toAssembly(genome));
            currentGenomeId = genome.getId();
        });
    }

    public static void clearChords() {
        if (view != null) {
            runOnEdt(view::clearChords);
        }
    }

    // ---- Internals ----------------------------------------------------------

    /**
     * Lazily create the singleton view and its window. Must run on the EDT.
     */
    private static CircularView getInstance() {
        if (view == null) {
            CircularViewConfig config = new CircularViewConfig();
            config.onChordClick = CircularViewUtilities::onChordClick;

            view = new CircularView(config);
            panel = new CircularViewPanel(view);

            frame = new JFrame("Circular View");
            frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
            frame.setContentPane(panel);
            frame.pack();
            frame.setLocationRelativeTo(IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null);
        }
        return view;
    }

    /**
     * Must run on the event thread
     */
    public static void open() {
        if (frame == null) {
            getInstance();
        }
        if (!frame.isVisible()) {
            frame.setVisible(true);
        }
        frame.toFront();

    }

    /**
     * Set the assembly from the current genome if it differs from what's shown.
     */
    private static void ensureAssembly(CircularView v) {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome != null && !genome.getId().equals(currentGenomeId)) {
            v.setAssembly(toAssembly(genome));
            currentGenomeId = genome.getId();
        }
    }

    /**
     * Convert an IGV genome to a circular-view {@link Assembly}. Only the "long"
     * (whole-genome) chromosomes are drawn; names are shortened so they match the
     * shortened {@code refName}s on the chords.
     */
    private static Assembly toAssembly(Genome genome) {
        List<Chromosome> chromosomes = new ArrayList<>();
        for (String chr : genome.getLongChromosomeNames()) {
            org.igv.feature.Chromosome c = genome.getChromosome(chr);
            String shortName = Chromosome.shortChrName(chr);
            chromosomes.add(new Chromosome(shortName, c.getLength(), ChrColors.getChrColor(shortName)));
        }
        return new Assembly(genome.getDisplayName(), genome.getId(), chromosomes);
    }

    /**
     * Navigate IGV to a clicked chord's two regions (with flanking), shown
     * side-by-side. Port of the onChordClick callback in circularView.js.
     */
    private static void onChordClick(Chord feature) {
        if (!IGV.hasInstance()) {
            return;
        }
        Mate mate = feature.getMate();
        if (mate == null) {
            return;
        }
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        String locus1 = locusString(genome, feature.getRefName(), feature.getStart(), feature.getEnd());
        String locus2 = locusString(genome, mate.getRefName(), mate.getStart(), mate.getEnd());
        IGV.getInstance().goToLocus(locus1 + " " + locus2);
    }

    private static String locusString(Genome genome, String refName, long start, long end) {
        String chr = (genome != null) ? genome.getCanonicalChrName(refName) : refName;
        long s = Math.max(0, start - CLICK_FLANKING);
        long e = end + CLICK_FLANKING;
        return chr + ":" + s + "-" + e;
    }

    private static void runOnEdt(Runnable r) {
        if (SwingUtilities.isEventDispatchThread()) {
            r.run();
        } else {
            SwingUtilities.invokeLater(r);
        }
    }
}
