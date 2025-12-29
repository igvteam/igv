package org.igv.feature.genome;

import org.igv.feature.CytoBandFileParser;
import org.igv.feature.Cytoband;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

public class CytobandMap implements CytobandSource {

    private LinkedHashMap<String, List<Cytoband>> cytobandMap;

    CytobandMap(String path) {
        try {
            init(path);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public CytobandMap(LinkedHashMap<String, List<Cytoband>> cytobandMap) {
        this.cytobandMap = cytobandMap;
    }

    void init(String path) throws IOException {
        try(BufferedReader br = ParsingUtils.openBufferedReader(path)) {
            this.cytobandMap = CytoBandFileParser.loadData(br);
        }
    }

    @Override
    public List<Cytoband> getCytobands(String chr) {
        return cytobandMap.get(chr);
    }

}
