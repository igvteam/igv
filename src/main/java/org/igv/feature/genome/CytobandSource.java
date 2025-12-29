package org.igv.feature.genome;

import org.igv.feature.Cytoband;

import java.io.IOException;
import java.util.List;

public interface CytobandSource {

    List<Cytoband> getCytobands(String chr) throws IOException;
}