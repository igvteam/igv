package org.broad.igv.ucsc;

import org.broad.igv.Globals;
import org.broad.igv.feature.Locus;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Search for a feature from the UCSC or IGV web services.  The UCSC service can be slow, for the time being
 * the IGV service is preferred.
 *
 * <p>
 * This class is not thread safe.  It is not intended to be used in a multi-threaded environment.
 */
public class SearchAPI {


    public static List<Locus> search(String searchTerm, String genome) throws IOException {

        List<String> positions = searchIGV(searchTerm, genome);  //searchUCSC(searchTerm, genome);

        return mergeOverlaps(positions);
    }

    private static List<String> searchIGV(String str, String genomeID) throws IOException {

        String tmp = "https://igv.org/genomes/locus.php?genome=$GENOME$&name=$FEATURE$";

        if (genomeID != null && genomeID.indexOf("/") < 0 && genomeID.indexOf("\\") < 0) {   // Filter out file paths
            URL url = new URL(tmp.replace("$GENOME$", genomeID).replace("$FEATURE$", str));
            String r = HttpUtils.getInstance().getContentsAsString(url);
            return List.of(Globals.whitespacePattern.split(r));
        } else {
            return null;
        }
    }


    private static List<String> searchUCSC(String searchTerm, String genome) throws IOException {
        String url = "https://api.genome.ucsc.edu/search?search=" + searchTerm + "&genome=" + genome;
        String response = HttpUtils.getInstance().getContentsAsString(new URL(url));
        org.json.JSONObject json = new org.json.JSONObject(response);
        return reduceSearchResults(json.toMap(), searchTerm);
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
        String searchTerm = "PKD1";
        String genome = "hg38";
        List<Locus> results = search(searchTerm, genome);
        for (Locus result : results) {
            System.out.println(result);
        }
    }
}
