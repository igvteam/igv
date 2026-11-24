package figwig;

import java.io.Closeable;
import java.io.IOException;
import java.util.*;

/**
 * Reads windows of a figwig file into arrays, after the FigWigReader of the figwig Python package
 * (docs/figwig.rst of https://github.com/jmschrei/figwig): many windows of one width in one call, from
 * any subset of the file's tracks, into an array of shape (windows, tracks, width).  It is built on the
 * public methods of {@link FigWigFile} alone, and can be run from the command line ({@link #main}).
 *
 * Where it differs from the Python reader: values are returned as double, which holds every figwig type
 * exactly, whatever type they were written as; a base past the end of its chromosome is NaN; and sparse
 * reads, sums per label or bin, GPU output and remote files are left out.
 */
public class FigWigReader implements Closeable {

    private final FigWigFile file;
    private final int[] held;   // the file's index of each of this reader's tracks

    /** A reader of every track of a local file. */
    public FigWigReader(String path) throws IOException {
        this(FigWigFile.open(path), null);
    }

    /** A reader of some tracks of a local file, by name, in the order given. */
    public FigWigReader(String path, String... trackNames) throws IOException {
        this(FigWigFile.open(path), trackNames);
    }

    /**
     * A reader of some tracks of an open file, by name, in the order given; of every track when
     * {@code trackNames} is null.  The reader then acts as if the file held only those.  Closing the
     * reader releases the file, which closes it when nothing else has retained it.
     */
    public FigWigReader(FigWigFile file, String[] trackNames) throws IOException {
        this.file = file;
        if (trackNames == null) {
            held = new int[file.nTracks()];
            for (int t = 0; t < held.length; t++) {
                held[t] = t;
            }
        } else {
            Map<String, Integer> byName = new HashMap<>();
            for (int t = file.nTracks() - 1; t >= 0; t--) {
                byName.put(file.trackName(t), t);
            }
            held = new int[trackNames.length];
            for (int i = 0; i < held.length; i++) {
                Integer t = byName.get(trackNames[i]);
                if (t == null) {
                    file.close();
                    throw new IllegalArgumentException("figwig: " + file.path() + " has no track named " + trackNames[i]);
                }
                held[i] = t;
            }
        }
        for (int t : held) {
            file.retain(t);
        }
    }

    // ------------------------------------------------------------------------------------------
    // What the file holds

    /** The reference genome the file names, or null. */
    public String genome() {
        return file.genome();
    }

    /** Bases per tile, the bin of {@code tileSums} at level 0. */
    public int tile() {
        return file.summaryBases();
    }

    /** The length of each chromosome, in the file's order. */
    public Map<String, Long> chromSizes() {
        Map<String, Long> sizes = new LinkedHashMap<>();
        for (int c = 0; c < file.nChroms(); c++) {
            sizes.put(file.chromName(c), file.chromLength(c));
        }
        return sizes;
    }

    /** The number of tracks this reader holds. */
    public int nTracks() {
        return held.length;
    }

    /** A record per track: its name, its dtype and its attributes. */
    public List<Map<String, String>> tracks() {
        List<Map<String, String>> records = new ArrayList<>(held.length);
        for (int t : held) {
            Map<String, String> record = new LinkedHashMap<>();
            record.put("name", file.trackName(t));
            record.put("dtype", file.trackDtype(t));
            record.putAll(file.trackAttributes(t));
            records.add(record);
        }
        return records;
    }

    /** Per track, in the order of {@code tracks()}: its sum, nonzero bases, least and greatest value. */
    public double[][] stats() {
        double[][] stats = new double[held.length][];
        for (int i = 0; i < held.length; i++) {
            int t = held[i];
            stats[i] = new double[]{file.trackSum(t), file.trackNonzero(t), file.trackMin(t), file.trackMax(t)};
        }
        return stats;
    }

    /** The zoom levels the file keeps besides its tile sums, in tiles per bin: 4, 16, 64 for three. */
    public int[] zooms() {
        int[] zooms = new int[file.levels() - 1];
        for (int k = 1; k <= zooms.length; k++) {
            zooms[k - 1] = (int) (file.levelBases(k) / file.summaryBases());
        }
        return zooms;
    }

    /**
     * The indices of the tracks that have every one of the attributes given as key, value, key, value...,
     * in the reader's order, to pass as {@code tracks}.  "name" and "dtype" can be matched too.
     */
    public int[] select(String... keysAndValues) {
        if (keysAndValues.length % 2 != 0) {
            throw new IllegalArgumentException("figwig: select takes keys and values in pairs");
        }
        List<Map<String, String>> records = tracks();
        int[] matched = new int[records.size()];
        int n = 0;
        for (int i = 0; i < matched.length; i++) {
            boolean all = true;
            for (int k = 0; k < keysAndValues.length; k += 2) {
                all &= keysAndValues[k + 1].equals(records.get(i).get(keysAndValues[k]));
            }
            if (all) {
                matched[n++] = i;
            }
        }
        return Arrays.copyOf(matched, n);
    }

    // ------------------------------------------------------------------------------------------
    // Reading

    /** Every track over each window: an array of shape (windows, tracks, width). */
    public double[][][] read(String[] chroms, int[] starts, int width) throws IOException {
        return read(chroms, starts, width, null);
    }

    /**
     * The tracks given, as indices into {@code tracks()}, over windows of {@code width} bases: window j
     * covers bases starts[j] to starts[j] + width - 1 of chroms[j], 0-based.  The result has shape
     * (windows, tracks, width), with the tracks in the order given.  A base that was never written is 0,
     * and a base past the end of its chromosome is NaN.  The windows may overlap, repeat and come in any
     * order.  A window that starts before 0, or lies on a chromosome the file does not have, is refused.
     */
    public double[][][] read(String[] chroms, int[] starts, int width, int[] tracks) throws IOException {
        if (chroms.length != starts.length) {
            throw new IllegalArgumentException("figwig: " + chroms.length + " chromosomes for " + starts.length + " starts");
        }
        if (width <= 0) {
            throw new IllegalArgumentException("figwig: width " + width);
        }
        if (tracks == null) {
            tracks = new int[held.length];
            for (int i = 0; i < tracks.length; i++) {
                tracks[i] = i;
            }
        }
        double[][][] out = new double[starts.length][tracks.length][width];
        for (int j = 0; j < starts.length; j++) {
            Integer c = file.chromId(chroms[j]);
            if (c == null) {
                throw new IllegalArgumentException("figwig: " + file.path() + " has no chromosome " + chroms[j]);
            }
            if (starts[j] < 0) {
                throw new IllegalArgumentException("figwig: a window starts at " + starts[j]);
            }
            long length = file.chromLength(c);
            for (int i = 0; i < tracks.length; i++) {
                double[] row = out[j][i];
                FigWigFile.Runs runs = file.read(held[tracks[i]], c, starts[j], starts[j] + width);
                for (int r = 0; r < runs.n; r++) {
                    Arrays.fill(row, runs.starts[r] - starts[j], runs.ends[r] - starts[j], runs.values[r]);
                }
                if (starts[j] + (long) width > length) {
                    Arrays.fill(row, (int) Math.max(0, length - starts[j]), width, Double.NaN);
                }
            }
        }
        return out;
    }

    /** One track over each window: an array of shape (windows, width). */
    public double[][] read(String[] chroms, int[] starts, int width, int track) throws IOException {
        double[][][] y = read(chroms, starts, width, new int[]{track});
        double[][] out = new double[y.length][];
        for (int j = 0; j < y.length; j++) {
            out[j] = y[j][0];
        }
        return out;
    }

    /**
     * A track's sum over every bin of a chromosome, as the file keeps it, without reading its values:
     * level 0 is the sum per tile, and level k the sum per bin of zooms()[k - 1] tiles.  Fails when the
     * file keeps no sums for the track, as for a group of thousands of cells.
     */
    public double[] tileSums(int track, String chrom, int level) {
        Integer c = file.chromId(chrom);
        if (c == null) {
            throw new IllegalArgumentException("figwig: " + file.path() + " has no chromosome " + chrom);
        }
        if (level < 0 || level >= file.levels()) {
            throw new IllegalArgumentException("figwig: " + file.path() + " has no level " + level);
        }
        double[] sums = file.levelSums(held[track], level, c);
        if (sums == null) {
            throw new IllegalStateException("figwig: " + file.path() + " keeps no sums for track " + file.trackName(held[track]));
        }
        return sums;
    }

    /** Release this reader's tracks; the file closes when nothing else has retained it. */
    @Override
    public void close() throws IOException {
        for (int t : held) {
            file.release(t);
        }
    }

    // ------------------------------------------------------------------------------------------
    // Command line

    private static final String USAGE = """
            usage: figwig.FigWigReader info FILE
                   figwig.FigWigReader read FILE CHROM START WIDTH [TRACK ...]
            info prints the file's genome, tile, chromosomes and tracks.
            read prints a line per base of the window: its position, then its value in each track
            named (every track when none is)."""
            ;

    public static void main(String[] args) {
        try {
            run(args);
        } catch (java.nio.file.NoSuchFileException e) {
            System.err.println("figwig: no such file: " + e.getMessage());
            System.exit(1);
        } catch (IOException | java.io.UncheckedIOException | IllegalArgumentException | IllegalStateException e) {
            System.err.println(e instanceof NumberFormatException ? "figwig: START and WIDTH must be whole numbers" : e.getMessage());
            System.exit(1);
        }
    }

    private static void run(String[] args) throws IOException {
        if (args.length == 2 && args[0].equals("info")) {
            try (FigWigReader fw = new FigWigReader(args[1])) {
                System.out.println("genome\t" + fw.genome());
                System.out.println("tile\t" + fw.tile());
                System.out.println("zooms\t" + Arrays.toString(fw.zooms()));
                System.out.println("chromosomes\t" + fw.chromSizes().size());
                List<Map<String, String>> tracks = fw.tracks();
                double[][] stats = fw.stats();
                System.out.println("track\tsum\tnonzero\tmin\tmax\trecord");
                for (int i = 0; i < tracks.size(); i++) {
                    System.out.println(i + "\t" + stats[i][0] + "\t" + (long) stats[i][1] + "\t" + stats[i][2]
                            + "\t" + stats[i][3] + "\t" + tracks.get(i));
                }
            }
        } else if (args.length >= 5 && args[0].equals("read")) {
            String[] names = args.length > 5 ? Arrays.copyOfRange(args, 5, args.length) : null;
            int start = Integer.parseInt(args[3]), width = Integer.parseInt(args[4]);
            try (FigWigReader fw = new FigWigReader(FigWigFile.open(args[1]), names)) {
                double[][] y = fw.read(new String[]{args[2]}, new int[]{start}, width)[0];
                StringBuilder line = new StringBuilder("position");
                for (Map<String, String> track : fw.tracks()) {
                    line.append('\t').append(track.get("name"));
                }
                System.out.println(line);
                for (int p = 0; p < width; p++) {
                    line = new StringBuilder().append(start + p);
                    for (double[] track : y) {
                        line.append('\t').append(track[p]);
                    }
                    System.out.println(line);
                }
            }
        } else {
            System.err.println(USAGE);
            System.exit(2);
        }
    }
}
