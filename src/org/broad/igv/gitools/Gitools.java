package org.broad.igv.gitools;

import org.apache.commons.math.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.collections.DoubleArrayList;

import java.io.*;
import java.util.*;

/**
 * @author jrobinso
 *         Date: 11/13/12
 *         Time: 9:26 PM
 */
public class Gitools {


    public static void exportTDM(List<String> lociStrings, File file) throws IOException {

        // Convert the loci strings to a list of lodi, if the loci represents multiple features (e.g. isoforms) use the largest
        int averageFeatureSize = 0;
        List<NamedFeature> loci = new ArrayList<NamedFeature>(lociStrings.size());
        for (String l : lociStrings) {
            NamedFeature feature = FeatureDB.getFeature(l);
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
            if ((t instanceof DataTrack || t.getTrackType() == TrackType.MUTATION) && t.isVisible()) {
                loadedTypes.add(t.getTrackType());
                samples.add(t.getSample());
            }
        }

        // set an appropriate zoom level for this feature set.  Very approximate
        int zoom = 0;
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome != null) {
            double averageChrLength = genome.getTotalLength() / genome.getChromosomes().size();
            zoom = (int) (Math.log(averageChrLength / averageFeatureSize) / Globals.log2) + 1;
        }


        // Loop though tracks and loci and gather data by sample & gene
        Map<String, SampleData> sampleDataMap = new LinkedHashMap<String, SampleData>();
        for (Track t : tracks) {
            if(!t.isVisible()) continue;

            if (t instanceof DataTrack) {
                DataTrack dataTrack = (DataTrack) t;

                for (NamedFeature locus : loci) {
                    double regionScore = dataTrack.getAverageScore(locus.getChr(), locus.getStart(), locus.getEnd(), zoom);
                    if (!Double.isNaN(regionScore)) {
                        String sample = t.getSample();
                        String locusString =  locus.getName();
                        String key = sample + "_" + locusString;

                        SampleData sd = sampleDataMap.get(key);
                        if (sd == null) {

                            sd = new SampleData(sample, locusString);
                            sampleDataMap.put(key, sd);
                        }
                        sd.addValue(t.getTrackType(), regionScore);
                    }
                }
            }
        }



        // Finally output data
        PrintWriter pw = null;

        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));

            pw.print("Sample\tLocus");
            for (TrackType tt : loadedTypes) {
                pw.print("\t" + tt.name());
            }
            pw.println();

            for (SampleData sd : sampleDataMap.values()) {
                pw.print(sd.sample + "\t" + sd.locus);
                for (TrackType tt : loadedTypes) {
                    pw.print('\t');
                    double[] values = sd.getValues(tt);
                    if (values == null) {
                        pw.print("-");
                    } else {
                        double avg = StatUtils.max(values);
                        pw.print(avg);
                    }
                }
                pw.println();
            }
        } finally {
            if(pw != null) pw.close();
        }
    }

    static class SampleData {

        String sample;
        String locus;
        Map<TrackType, DoubleArrayList> valueMap;

        SampleData(String sample, String locus) {
            this.sample = sample;
            this.locus = locus;
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
