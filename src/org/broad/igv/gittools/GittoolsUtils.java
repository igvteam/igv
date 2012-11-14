package org.broad.igv.gittools;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.collections.DoubleArrayList;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author jrobinso
 *         Date: 11/13/12
 *         Time: 9:26 PM
 */
public class GittoolsUtils {


    public static void exportTDM(List<String> lociStrings) {

        // Convert the loci strings to a list of lodi, if the loci represents multiple features (e.g. isoforms) use the largest
        int averageFeatureSize = 0;
        List<Feature> loci = new ArrayList<Feature>(lociStrings.size());
        for (String l : lociStrings) {
            Feature feature = FeatureDB.getLongestFeatureNamed(l);
            if (feature == null) {
                feature = Locus.fromString(l);
            }
            if (feature != null) {
                loci.add(feature);
                averageFeatureSize += (feature.getEnd() - feature.getStart());
            }
        }
        if (loci.size() > 0) averageFeatureSize /= loci.size();

        // Determine data types -- all data tracks + mutation, and samples
        LinkedHashSet<TrackType> loadedTypes = new LinkedHashSet<TrackType>();
        LinkedHashSet<String> samples = new LinkedHashSet<String>();

        List<Track> tracks = IGV.getInstance().getAllTracks();
        for (Track t : tracks) {
            if (t instanceof DataTrack || t.getTrackType() == TrackType.MUTATION) {
                loadedTypes.add(t.getTrackType());
                samples.add(t.getSample());
            }
        }

        // set an appropriate zoom level for this feature set.  Very approximate
        int zoom = 0;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome != null) {
            double averageChrLength = genome.getLength() / genome.getChromosomes().size();
            zoom = (int) (Math.log(averageChrLength / averageFeatureSize) / Globals.log2) + 1;
        }


        // Loop though tracks and loci and gather data by sample
        Map<String, SampleData> sampleDataMap = new LinkedHashMap<String, SampleData>();
        for (Track t : tracks) {
            if (t instanceof DataTrack) {
                DataTrack dataTrack = (DataTrack) t;

                for (Feature locus : loci) {
                    double regionScore = dataTrack.getAverageScore(locus.getChr(), locus.getStart(), locus.getEnd(), zoom);
                    if (!Double.isNaN(regionScore)) {
                        String sample = t.getSample();
                        SampleData sd = sampleDataMap.get(sample);
                        if (sd == null) {
                            sd = new SampleData(sample);
                            sampleDataMap.put(sample, sd);
                        }
                        sd.addValue(t.getTrackType(), regionScore);
                    }
                }
            }
        }

        // Finally output data
        System.out.print("Sample");
        for (TrackType tt : loadedTypes) {
            System.out.print("\t" + tt.name());
        }
        System.out.println();

        for (Map.Entry<String, SampleData> entry : sampleDataMap.entrySet()) {
            String sample = entry.getKey();
            SampleData sd = entry.getValue();
            System.out.print(sample);
            for (TrackType tt : loadedTypes) {
                System.out.print('\t');
                double[] values = sd.getValues(tt);
                if (values == null) {
                    // Print nothing, or "null" indicator?
                } else {
                    double avg = StatUtils.max(values);
                    System.out.print(avg);
                }
            }
        }
    }

    static class SampleData {
        String name;
        Map<TrackType, DoubleArrayList> valueMap;

        SampleData(String name) {
            this.name = name;
            this.valueMap = new HashMap<TrackType, DoubleArrayList>();
        }

        void addValue(TrackType tt, double value) {

            DoubleArrayList vlist = valueMap.get(tt);
            if (vlist == null) {
                vlist = new DoubleArrayList();
                valueMap.put(tt, vlist);
            }
            vlist.add(value);
        }

        public double[] getValues(TrackType tt) {
            DoubleArrayList dal = valueMap.get(tt);
            return dal == null ? null : dal.toArray();
        }
    }
}
