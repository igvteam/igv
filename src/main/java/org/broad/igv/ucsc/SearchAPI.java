package org.broad.igv.ucsc;

import com.google.gson.Gson;
import org.broad.igv.feature.Locus;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

// e.g. https://api.genome.ucsc.edu/search?search=brca1&genome=hg38
public class SearchAPI {

    public static List<Locus> search(String searchTerm, String genome) throws IOException {
        String url = "https://api.genome.ucsc.edu/search?search=" + searchTerm + "&genome=" + genome;
        String response = HttpUtils.getInstance().getContentsAsString(new URL(url));
        Gson gson = new Gson();
        return mergeOverlaps(reduceSearchResults(gson.fromJson(response, Map.class), searchTerm));
    }

    public static List<String> reduceSearchResults(Map<String, Object> searchResults, String searchTerm) {
        List<String> reducedResults = new ArrayList<>();
        List<Map<String, Object>> results = (List<Map<String, Object>>) searchResults.get("positionMatches");
        if (results != null) {
            for (Map<String, Object> result : results) {
                List<Map<String, Object>> matches = (List<Map<String, Object>>) result.get("matches");
                if (matches != null) {
                    for (Map<String, Object> match : matches) {
                        if (searchTerm.equalsIgnoreCase(match.get("posName").toString())) {
                            reducedResults.add(match.get("position").toString());
                        }
                    }
                }
            }
        }
        return reducedResults;
    }

    private static List<Locus> mergeOverlaps(List<String> positions) {

        List<Locus> loci = positions.stream()
                .map(p -> Locus.fromString(p))
                .filter(locus -> locus != null)
                .collect(Collectors.toList());

        loci.sort((a, b) -> {
            int chrComparison = a.getChr().compareTo(b.getChr());
            return chrComparison != 0 ? chrComparison : Integer.compare(a.getStart(), b.getStart());
        });

        List<Locus> mergedLoci = new ArrayList<>();
        for (Locus locus : loci) {
            if (!mergedLoci.isEmpty() && mergedLoci.get(mergedLoci.size() - 1).overlaps(locus)) {
                mergedLoci.get(mergedLoci.size() - 1).merge(locus);
            } else {
                mergedLoci.add(locus);
            }
        }

        return mergedLoci;

    }

    public static void main(String[] args) throws IOException {
        String searchTerm = "rs121913433";
        String genome = "hg38";
        List<Locus> results = search(searchTerm, genome);
        for (Locus result : results) {
            System.out.println(result);
        }
    }
}
