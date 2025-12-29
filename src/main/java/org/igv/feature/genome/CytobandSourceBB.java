package org.igv.feature.genome;

import org.igv.feature.BasicFeature;
import org.igv.feature.Cytoband;
import org.igv.feature.IGVFeature;
import org.igv.ucsc.bb.BBFeatureSource;
import org.igv.ucsc.bb.BBFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class CytobandSourceBB implements CytobandSource {

    private BBFeatureSource featureSource;

    public CytobandSourceBB(String path, Genome genome) throws IOException {
        BBFile bbfile = new BBFile(path, genome);
        featureSource = new BBFeatureSource(bbfile, genome);
    }

    @Override
    public List<Cytoband> getCytobands(String chr) throws IOException {
        List<Cytoband> cytobands = new ArrayList<>();
        Iterator<IGVFeature> features = featureSource.getFeatures(chr, 0, Integer.MAX_VALUE);
        while (features.hasNext()) {
            IGVFeature f = features.next();
            cytobands.add(new Cytoband(f.getChr(), f.getStart(), f.getEnd(), f.getName(), f.getAttribute("gieStain")));
        }
        return cytobands;
    }

}

