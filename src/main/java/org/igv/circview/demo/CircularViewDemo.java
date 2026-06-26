package org.igv.circview.demo;

import org.igv.circview.model.Assembly;
import org.igv.circview.model.Chord;
import org.igv.circview.model.Chromosome;
import org.igv.circview.model.Mate;
import org.igv.circview.ui.CircularView;
import org.igv.circview.ui.CircularViewConfig;
import org.igv.circview.ui.CircularViewPanel;
import org.igv.circview.util.ChrColors;
import org.igv.circview.util.ColorUtils;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

/**
 * Standalone demo: shows the circular view for hg19 with a few sample chords.
 * Clicking a chord prints its feature to the console.
 */
public final class  CircularViewDemo {

    public static void main(String[] args) {
        SwingUtilities.invokeLater(CircularViewDemo::createAndShow);
    }

    private static void createAndShow() {
        CircularViewConfig config = new CircularViewConfig();
        config.onChordClick = feature -> System.out.println("Chord clicked: " + feature);

        CircularView view = new CircularView(config);
        view.setAssembly(hg19());
        view.addChords(sampleChords(), "Structural variants",
                ColorUtils.parseColor("rgba(0, 0, 255, 0.35)"));
        view.addChords(moreChords(), "Translocations",
                ColorUtils.parseColor("rgba(220, 0, 0, 0.45)"));

        CircularViewPanel panel = new CircularViewPanel(view);

        JFrame frame = new JFrame("Circular View (Java) — hg19 demo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setContentPane(panel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    /** hg19 assembly, ported from test/hg19.js (names already short). */
    private static Assembly hg19() {
        long[][] data = {
                {1, 248956422}, {2, 242193529}, {3, 198295559}, {4, 190214555},
                {5, 181538259}, {6, 170805979}, {7, 159345973}, {8, 145138636},
                {9, 138394717}, {10, 133797422}, {11, 135086622}, {12, 133275309},
                {13, 114364328}, {14, 107043718}, {15, 101991189}, {16, 90338345},
                {17, 83257441}, {18, 80373285}, {19, 58617616}, {20, 64444167},
                {21, 46709983}, {22, 50818468},
        };
        List<Chromosome> chromosomes = new ArrayList<>();
        for (long[] d : data) {
            String name = String.valueOf(d[0]);
            chromosomes.add(new Chromosome(name, d[1], ChrColors.getChrColor(name)));
        }
        chromosomes.add(new Chromosome("X", 156040895, ChrColors.getChrColor("X")));
        chromosomes.add(new Chromosome("Y", 57227415, ChrColors.getChrColor("Y")));
        return new Assembly("Human hg19", "hg19", chromosomes);
    }

    /** A handful of intra- and inter-chromosomal chords. */
    private static List<Chord> sampleChords() {
        List<Chord> chords = new ArrayList<>();
        chords.add(chord("1", 30_000_000, 30_200_000, "8", 95_000_000, 95_200_000));
        chords.add(chord("1", 130_000_000, 130_050_000, "1", 200_000_000, 200_050_000));
        chords.add(chord("3", 50_000_000, 50_100_000, "17", 40_000_000, 40_100_000));
        chords.add(chord("5", 12_000_000, 12_100_000, "12", 60_000_000, 60_100_000));
        chords.add(chord("X", 70_000_000, 70_300_000, "7", 100_000_000, 100_300_000));
        chords.add(chord("22", 20_000_000, 20_100_000, "9", 21_000_000, 21_100_000));
        chords.add(chord("2", 90_000_000, 90_200_000, "2", 180_000_000, 180_200_000));
        return chords;
    }

    /** A second chord set, so the controls show more than one row. */
    private static List<Chord> moreChords() {
        List<Chord> chords = new ArrayList<>();
        chords.add(chord("4", 60_000_000, 60_100_000, "11", 70_000_000, 70_100_000));
        chords.add(chord("6", 40_000_000, 40_100_000, "20", 30_000_000, 30_100_000));
        chords.add(chord("10", 50_000_000, 50_100_000, "13", 80_000_000, 80_100_000));
        chords.add(chord("X", 20_000_000, 20_100_000, "16", 10_000_000, 10_100_000));
        return chords;
    }

    private static Chord chord(String ref, long start, long end,
                               String mateRef, long mateStart, long mateEnd) {
        String id = ref + ":" + start + "-" + end + "_" + mateRef + ":" + mateStart + "-" + mateEnd;
        Color color = null; // inherit the chord set color
        return new Chord(id, ref, start, end, new Mate(mateRef, mateStart, mateEnd), color);
    }
}
