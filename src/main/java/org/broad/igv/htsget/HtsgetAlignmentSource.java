package org.broad.igv.htsget;

import htsjdk.samtools.SAMFileHeader;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureSource;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class HtsgetAlignmentSource implements FeatureSource {
    HtsgetReader htsgetReader;
    SAMFileHeader header;
    Map<String, String> chrAliasMap;
    int featureWindowSize;
    Genome genome;

    public HtsgetAlignmentSource(HtsgetReader htsgetReader, Genome genome) {
        this.htsgetReader = htsgetReader;
        this.featureWindowSize = 1000;  // default
        this.chrAliasMap = new HashMap<>();
        this.genome = genome;
        init(genome);
    }

    private void init(Genome genome) {
        SAMFileHeader header = (SAMFileHeader) getHeader();
        if (genome != null) {
          // TODO: How does chrAliasMap looks like for Alignments?
        }
    }

    @Override
    public Iterator getFeatures(String chr, int start, int end) throws IOException {
        return null;
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    @Override
    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) { this.featureWindowSize = size; }
}
