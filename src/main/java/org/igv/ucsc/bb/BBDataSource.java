package org.igv.ucsc.bb;

import htsjdk.samtools.util.Locatable;
import org.igv.Globals;
import org.igv.data.AbstractDataSource;
import org.igv.data.BasicScore;
import org.igv.data.DataSource;
import org.igv.data.DataTile;
import org.igv.feature.Chromosome;
import org.igv.feature.LocusScore;
import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.track.TrackType;
import org.igv.track.WindowFunction;

import java.io.IOException;
import java.util.*;


public class BBDataSource extends AbstractDataSource implements DataSource {

    private static Logger log = LogManager.getLogger(BBDataSource.class);

    private final Genome genome;
    Collection<WindowFunction> availableWindowFunctions =
            Arrays.asList(WindowFunction.min, WindowFunction.mean, WindowFunction.max, WindowFunction.none);

    BBFile reader;

    private Map<WindowFunction, List<LocusScore>> wholeGenomeScores;

    private double dataMin = 0;
    private double dataMax = 100;

    List<LocusScore> wgScores = null;

    public BBDataSource(BBFile reader, Genome genome) throws IOException {
        super(genome);
        this.reader = reader;
        this.genome = genome;
        this.wholeGenomeScores = new HashMap<>();
    }

    public double getDataMax() {
        return dataMax;
    }

    public double getDataMin() {
        return dataMin;
    }

    @Override
    protected DataTile getRawData(String chr, int start, int end) {

        try {
            long rTreeOffset = reader.getHeader().fullIndexOffset;
            Integer chrIdx = reader.getIdForChr(chr);
            if (chrIdx == null) {
                return null;
            }
            List<byte[]> chunks = this.reader.getLeafChunks(chrIdx, start, chrIdx, end, rTreeOffset);


            List<LocusScore> features = new ArrayList<>();
            for (byte[] c : chunks) {
                reader.decodeWigData(c, chrIdx, start, end, features);
            }

            final int size = features.size();
            int[] starts = new int[size];
            int[] ends = new int[size];
            float[] values = new float[size];
            for (int i = 0; i < size; i++) {
                final LocusScore locusScore = features.get(i);
                starts[i] = locusScore.getStart();
                ends[i] = locusScore.getEnd();
                values[i] = locusScore.getScore();
            }
            return new DataTile(starts, ends, values, null);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }

    }

    /**
     * Return bigwig "zoom data" if available for the resolution encoded by "zoom"
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     */
    @Override
    protected List<LocusScore> getPrecomputedSummaryScores(String chr, int start, int end, int zoom) {
        // Translate IGV zoom number to bpperpixel.
        try {
            if (Globals.CHR_ALL.equals(chr)) {
                if (genome.getHomeChromosome().equals(Globals.CHR_ALL) && windowFunction != WindowFunction.none) {
                    return this.getWholeGenomeScores();
                } else {
                    return null;
                }
            }


            Chromosome chromosome = genome.getChromosome(chr);
            if (chromosome == null) {
                throw new RuntimeException("Unexpected chromosome name: " + chr);
            }


            double nBins = Math.pow(2, zoom);
            double scale = chromosome.getLength() / (nBins * 700);
            BBZoomHeader zlHeader = reader.zoomLevelForScale(scale);

            if (zlHeader == null) {
                return null;
            } else {
                long rTreeOffset = zlHeader.indexOffset;
                Integer chrIdx = reader.getIdForChr(chr);
                if (chrIdx == null) {
                    return Collections.EMPTY_LIST;
                }
                List<byte[]> chunks = this.reader.getLeafChunks(chrIdx, start, chrIdx, end, rTreeOffset);
                List<LocusScore> features = new ArrayList<>();
                for (byte[] chunk : chunks) {
                    reader.decodeZoomData(chr, chunk, chrIdx, start, end, windowFunction, features);
                }
                return features;
            }
        } catch (IOException e) {
            // Wrap the IOException for now, to avoid refactoring of hierarchy
            throw new RuntimeException(e);
        }
    }

    @Override
    public int getLongestFeature(String chr) {
        return 0;
    }


    public TrackType getTrackType() {
        return TrackType.OTHER;
    }


    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return availableWindowFunctions;
    }

    @Override
    public void dispose() {

    }


    private List<LocusScore> getWholeGenomeScores() {

        try {
            if (genome.getHomeChromosome().equals(Globals.CHR_ALL) && windowFunction != WindowFunction.none) {

                if (!wholeGenomeScores.containsKey(windowFunction)) {

                    int screenWidth = 1000;  // nominal
                    double scale = genome.getWGLength() / screenWidth;

                    List<String> longChromosomeNames = genome.getLongChromosomeNames();
                    if (longChromosomeNames.isEmpty()) {
                        return null;
                    }

                    // Find the min and max chromosome ids for whole genome view
                    int minID = Integer.MAX_VALUE;
                    int maxID = -1;
                    for (String chr : longChromosomeNames) {
                        Integer id = reader.getIdForChr(chr);
                        if (id == null) {
                            continue;
                        }
                        if (id < minID) {
                            minID = id;
                        }
                        if (id > maxID) {
                            maxID = id;
                        }
                    }

                    ArrayList<LocusScore> scores = new ArrayList<LocusScore>();
                    wholeGenomeScores.put(windowFunction, scores);

                    BBZoomHeader lowestResHeader = reader.zoomLevelForScale(scale, 1000);
                    if (lowestResHeader == null) return null;

                    Set<String> wgChrNames = new HashSet<>(genome.getLongChromosomeNames());

                    long rTreeOffset = lowestResHeader.indexOffset;
                    List<byte[]> chunks = this.reader.getLeafChunks(minID, 0, maxID, Integer.MAX_VALUE, rTreeOffset);

                    List<LocusScore> features = new ArrayList<>();
                    for (byte[] chunk : chunks) {
                        reader.decodeZoomData(null, chunk, -1, -1, -1, windowFunction, features);
                    }

                    for (LocusScore rec : features) {
                        String chr = genome.getCanonicalChrName(rec.getChr());
                        if (wgChrNames.contains(chr)) {
                            int genomeStart = genome.getGenomeCoordinate(chr, rec.getStart());
                            int genomeEnd = genome.getGenomeCoordinate(chr, rec.getEnd());
                            scores.add(new BasicScore(genomeStart, genomeEnd, rec.getScore()));
                        }
                    }
                    scores.sort(Comparator.comparingInt(Locatable::getStart));

                }
                return wholeGenomeScores.get(windowFunction);
            } else {
                return null;
            }
        } catch (IOException e) {
            log.error("Error getting whole genome scores", e);
            return null;
        }
    }

}
