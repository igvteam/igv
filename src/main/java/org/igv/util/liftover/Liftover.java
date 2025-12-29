package org.igv.util.liftover;

import org.igv.Globals;
import org.igv.feature.Range;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;


// TODO -- it is assumed that there is a single chain for each target sequenc (tName).  This is not neccessarily
// TODO    the case.  To fix this replace the tName -> chains map with a tName -> tree(chains) map.  The appropriate chain
// TODO    for a range would be fetched from the tree

/**
 * Class for mapping coordinates from one sequence to another based on the UCSC "chain" format for pairwise
 * alignments (https://genome.ucsc.edu/goldenPath/help/chain.html).
 */
public class Liftover {

    private static Logger log = LogManager.getLogger(Liftover.class);

    Map<String, Chain> chains;

    /**
     * Initialize with a list of "Chain" objects.  It is assumed there is one chain per target sequence chromosome ("tName")
     * @param chains
     */
    Liftover(List<Chain> chains) {
        this.chains = new HashMap<>();
        for (Chain c : chains) {
            this.chains.put(c.tName, c);
        }
    }

    /**
     * Map a range from target to query.  A single range can be split intol multiple query ranges
     *
     * @param range
     * @return
     */
    public List<Range> map(Range range) {

        if (chains.containsKey(range.getChr())) {
            List<Range> mapped = chains.get(range.getChr()).map(range);

            // Combine contiguous ranges
            List<Range> combined = new ArrayList<>();
            if (mapped.size() > 0) {
                Collections.sort(mapped, (o1, o2) -> {
                    if (o1.getChr().equals(o2.getChr())) {
                        return o1.getStart() - o2.getStart();
                    } else {
                        return o1.getChr().compareTo(o2.getChr());
                    }
                });
                Range lastRange = mapped.get(0);
                combined.add(lastRange);
                for (int i = 1; i < mapped.size(); i++) {
                    Range currentRange = mapped.get(i);
                    if (currentRange.getChr().equals(lastRange.getChr()) &&
                            lastRange.end == mapped.get(i).start) {
                        lastRange.end = currentRange.end;
                    } else {
                        lastRange = currentRange;
                        combined.add(lastRange);
                    }
                }
            }
            return combined;

        } else {
            return Collections.EMPTY_LIST;
        }
    }

    /**
     * Create a liftover from a chain file.  Chain files typically contain multiple chains.
     *
     * @param path
     * @return
     * @throws IOException
     */
    public static Liftover load(String path) throws IOException {
        BufferedReader br = null;
        try {
            br = ParsingUtils.openBufferedReader(path);

            List<Chain> chains = new ArrayList<>();
            Chain chain = null;
            List<String[]> data = new ArrayList<>();


            String line = "";
            while ((line = br.readLine()) != null) {
                line = line.strip();

                if (line.length() == 0) {
                    if (chain != null) {
                        chain.setAlignments(data);
                        chains.add(chain);
                        chain = null;
                        data.clear();
                    }
                } else if (line.startsWith("chain")) {
                    String[] t = Globals.whitespacePattern.split(line);
                    chain = new Chain(
                            t[2],
                            Integer.parseInt(t[3]),
                            Integer.parseInt(t[5]),
                            Integer.parseInt(t[6]),
                            t[7],
                            Integer.parseInt(t[8]),
                            Integer.parseInt(t[10]),
                            Integer.parseInt(t[11]),
                            t[12]);
                } else {
                    data.add(Globals.whitespacePattern.split(line));
                }
            }


            //Last data
            if (chain != null && data.size() > 0) {
                chain.setAlignments(data);
                chains.add(chain);
            }

            return new Liftover(chains);


        } finally {
            try {
                br.close();
            } catch (IOException e) {
                log.error("Error closing chaing file", e);
            }
        }
    }

}
