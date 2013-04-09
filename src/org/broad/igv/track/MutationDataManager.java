package org.broad.igv.track;

import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.genome.Genome;
import org.broad.tribble.Feature;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 8:49 AM
 */
public class MutationDataManager {

    String path;
    Interval currentInterval;
    Map<String, List<Mutation>> mutationMap = Collections.synchronizedMap(new HashMap());
    TribbleFeatureSource tribbleFeatureSource;

    public MutationDataManager(String path, Genome genome) throws IOException {
        this.path = path;

        tribbleFeatureSource = new TribbleFeatureSource(path, genome);
    }

    public synchronized Iterator<Mutation> getFeatures(String sample, String chr, int start, int end) throws IOException {
        if (currentInterval == null || !currentInterval.contains(chr, start, end)) {
            Iterator<Feature> mutations = tribbleFeatureSource.getFeatures(chr, start, end);

            while(mutations.hasNext()) {
                Mutation m = (Mutation) mutations.next();
                String sampleId = m.getSampleId();
                List<Mutation> sampleMutations = mutationMap.get(sampleId);
                if(sampleMutations == null) {
                    sampleMutations = new ArrayList<Mutation>();
                    mutationMap.put(sampleId, sampleMutations);
                }
                sampleMutations.add(m);
                currentInterval = new Interval(chr, start, end);
            }

        }
        List<Mutation> mutationList = mutationMap.get(sample);
        return mutationList == null ? Collections.EMPTY_LIST.iterator() : mutationList.iterator();

    }


    // TODO -- variants of this class exist elsewhere, centralize

    static class Interval {

        String chr;
        int start;
        int end;

        Interval(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
        }

        boolean contains(String chr, int start, int end) {
            return this.chr.equals(chr) &&
                    this.start <= start &&
                    this.end >= end;
        }
    }
}
