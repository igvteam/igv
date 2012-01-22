package org.broad.igv.sam;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.readers.AsciiLineReader;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date 11/29/11
 */
abstract public class BaseAlignmentCounts implements AlignmentCounts {

    private static Logger log = Logger.getLogger(BaseAlignmentCounts.class);

    private static char[] nucleotides = {'a', 'c', 'g', 't', 'n'};
    private static Map<String, Set<Integer>> knownSnps;

    public BaseAlignmentCounts() {
        final PreferenceManager prefs = PreferenceManager.getInstance();
         String snpsFile = prefs.get(PreferenceManager.KNOWN_SNPS, null);
        if (snpsFile != null && knownSnps == null) {
            loadKnownSnps(snpsFile);
        }
    }


    public String getValueStringAt(int pos) {

        StringBuffer buf = new StringBuffer();
        int totalCount = getTotalCount(pos);
        buf.append("Total count: " + totalCount + "<br>");
        for (char c : nucleotides) {
            int negCount = getNegCount(pos, (byte) c);
            int posCount = getPosCount(pos, (byte) c);
            int count = negCount + posCount;
            int percent = (int) Math.round(((float) count) * 100 / totalCount);
            char cU = Character.toUpperCase(c);
            buf.append(cU + "      : " + count);
            if (count == 0) {
                buf.append("<br>");
            } else {
                buf.append("  (" + percent + "%,     " + posCount + "+,   " + negCount + "- )<br>");
            }
        }

        int delCount = getDelCount(pos);
        if (delCount > 0) {
            buf.append("DEL: " + delCount);
        }
        int insCount = getInsCount(pos);
        if (insCount > 0) {
            buf.append("INS: " + insCount);
        }

        return buf.toString();

    }

    public boolean isMismatch(int pos, char ref, String chr, float snpThreshold) {

        Set<Integer> filteredSnps = knownSnps == null ? null : knownSnps.get(chr);
        if (filteredSnps == null || !filteredSnps.contains(pos + 1)) {
            float threshold = snpThreshold * getTotalQuality(pos);
            if (ref > 0) {
                for (char c : nucleotides) {
                    if (c != ref && c != 'n' && getQuality(pos, (byte) c) > threshold) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Load the set of known snps from a tab delimited file, format
     * chr < tab> location
     * The location is "1 base"  (first nucleotide is position 1).
     *
     * @param snpFile
     */
    private static synchronized void loadKnownSnps(String snpFile) {

        // This method might get called many times concurrently, but we only want to load these once.
        if (knownSnps != null) {
            return;
        }

        knownSnps = new HashMap();
        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(new ResourceLocator(snpFile));
            String nextLine = "";
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0];
                Set<Integer> snps = knownSnps.get(chr);
                if (snps == null) {
                    snps = new HashSet(10000);
                    knownSnps.put(chr, snps);
                }
                snps.add(new Integer(tokens[1]));
            }
        } catch (Exception e) {
            knownSnps = null;
            log.error("", e);
            MessageUtils.showMessage("Error loading snps file: " + snpFile + " (" + e.toString() + ")");
        } finally {
            reader.close();
        }


    }


}
