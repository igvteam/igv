package org.broad.igv.feature.sprite;

import org.apache.log4j.Logger;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by jrobinso on 6/30/18.
 */
public class ClusterParser {

    private static Logger log = Logger.getLogger(ClusterParser.class);

    public static ClusterSet parse(String file) throws IOException {

        int binSize = 1;         // Default => no binning
        List<Cluster> features = new ArrayList<>();

        BufferedReader br = null;

        br = ParsingUtils.openBufferedReader(file);

        String nextLine;
        while ((nextLine = br.readLine()) != null) {

            String[] tokens = ParsingUtils.TAB_PATTERN.split(nextLine);

            if(tokens[0].startsWith("#")) {
                if(tokens[0].startsWith("#binSize")) {
                    String [] t = ParsingUtils.EQ_PATTERN.split(tokens[0]);
                    binSize = Integer.parseInt(t[1].trim());
                }
                else {
                    continue;   // Comment or unrecognized directive
                }
            }

            if (tokens.length < 2) {
                log.info("Skipping line: " + nextLine);
                continue;
            }

            String name = tokens[0];
            Map<String, List<Integer>> positions = new HashMap<>();
            for(int i=1; i<tokens.length; i++) {
                String [] l = tokens[i].split(":");
                String chr = l[0];
                Integer position = Integer.parseInt(l[1]);

                List<Integer> posList = positions.get(chr);
                if(posList == null) {
                    posList = new ArrayList<>();
                    positions.put(chr, posList);
                }

                posList.add(position);
            }

            for(List<Integer> posList : positions.values()) {
                Collections.sort(posList);
            }

            features.add(new Cluster(name, positions));
        }

        return new ClusterSet(binSize, features);

    }


    static public class ClusterSet {

        int binSize;
        List<Cluster> clusters;

        public ClusterSet(int binSize, List<Cluster> clusters) {
            this.binSize = binSize;
            this.clusters = clusters;
        }
    }

}
