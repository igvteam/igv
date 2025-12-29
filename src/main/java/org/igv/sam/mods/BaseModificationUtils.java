package org.igv.sam.mods;

import htsjdk.samtools.util.SequenceUtil;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.sam.*;

import java.util.*;
import java.util.List;

public class BaseModificationUtils {

    private static Logger log = LogManager.getLogger(BaseModificationUtils.class);

    static Map<String, String> codeValues;

    static {
        codeValues = new HashMap<>();
        codeValues.put("m", "5mC");
        codeValues.put("h", "5hmC");
        codeValues.put("f", "5fC");
        codeValues.put("c", "5caC");
        codeValues.put("g", "5hmU");
        codeValues.put("e", "5fU");
        codeValues.put("b", "5caU");
        codeValues.put("a", "6mA");
        codeValues.put("o", "8xoG");
        codeValues.put("n", "Xao");
        codeValues.put("C", "Unknown C");
        codeValues.put("T", "Unknown T");
        codeValues.put("A", "Unknown A");
        codeValues.put("G", "Unknown G");
        codeValues.put("N", "Unknown");
        codeValues.put("NONE_C", "Unmodified C");
        codeValues.put("NONE_T", "Unmodified T");
        codeValues.put("NONE_G", "Unmodified G");
        codeValues.put("NONE_A", "Unmodified A");
    }

    private static int MM_WARNING_COUNT = 0;

    public static String modificationName(String modification) {
        return ((codeValues.containsKey(modification)) ? codeValues.get(modification) : modification);
    }


    /**
     * Parse the mm tag creating a base modification set for each modification listed.
     *
     * @param mm       MM tag value, string, examples
     *                 C+m?,5,12,0; :   A single modification, 1 set is returned
     *                 C+mh,5,12,0; :   Two modifications, 2 sets are returned
     *                 C+m,5,12,0;C+h,5,12,0;   Two modifications, 2 sets are returned
     * @param ml
     * @param sequence
     * @return
     */
    public static List<BaseModificationSet> getBaseModificationSets(String mm, byte[] ml, byte[] sequence, boolean isNegativeStrand) {

        if (isNegativeStrand) {
            sequence = AlignmentUtils.reverseComplementCopy(sequence);
        }

        List<BaseModificationSet> modificationSets = new ArrayList<>();


        String[] mmTokens = mm.split(";");
        int mlIdx = 0;      // likelihood array index

        for (String mmi : mmTokens) {

            String[] tokens = mmi.split(","); //Globals.commaPattern.split(mm);
            char base = tokens[0].charAt(0);
            char strand = tokens[0].charAt(1);
            boolean skippedBasesCalled;
            if (tokens[0].endsWith(".")) {
                skippedBasesCalled = true;
            } else if (tokens[0].endsWith("?")) {
                skippedBasesCalled = false;
            } else {
                skippedBasesCalled = PreferencesManager.getPreferences().getAsBoolean(Constants.BASEMOD_SKIPPED_BASES);
            }

            String context = base == 'C' ? PreferencesManager.getPreferences().get(Constants.BASEMOD_CYTOSINE_CONTEXT) : null;


            if (tokens.length == 1) {
                // Legal but not handled yet, indicates modification is not present.  Perhaps not relevant for visualization
            } else {

                String modificationString = tokens[0].endsWith(".") || tokens[0].endsWith("?") ?
                        tokens[0].substring(2, tokens[0].length() - 1) :
                        tokens[0].substring(2);

                // Parse modifications, this is rather complex, commensurate with the spec.  Unless a chebi code, modifications
                // are restricted to single characters, a multi-character string that is not a chebi code indicates
                // multiple modifications
                String[] modifications;
                if (modificationString.length() > 1) {
                    if (isChEBI(modificationString)) {
                        modifications = new String[]{modificationString};
                    } else {
                        modifications = new String[modificationString.length()];
                        for (int i = 0; i < modificationString.length(); i++) {
                            modifications[i] = modificationString.substring(i, i + 1);
                        }
                    }
                } else {
                    modifications = new String[]{modificationString};
                }


                // Create a positions -> likelihood map for each modification
                Map<String, Map<Integer, Byte>> likelihoodMap = new HashMap<>();
                for (String m : modifications) {
                    likelihoodMap.put(m, new HashMap<>());
                }

                int idx = 1;  // position array index,  positions start at index 1
                int skip = Integer.parseInt(tokens[idx++]);

                int p = 0;
                int matchCount = 0;

                while (p < sequence.length) {

                    if (base == 'N' || sequence[p] == base) {

                        if (base == 'C' && context != null && !BaseModificationSet.contextMatches(sequence, p, context)) {
                            // Context does not match, skip this base
                            p++;
                            continue;
                        }

                        int position = isNegativeStrand ? sequence.length - 1 - p : p;
                        if (matchCount == skip) { // && idx < tokens.length) {
                            for (String modification : modifications) {
                                byte likelihood = ml == null ? (byte) 255 : ml[mlIdx++];
                                likelihoodMap.get(modification).put(position, likelihood);
                            }
                            if (idx < tokens.length) {
                                skip = Integer.parseInt(tokens[idx++]);
                                matchCount = 0;
                            } else {
                                if (skippedBasesCalled) {
                                    // MM tag is exhausted, but continue scanning for skipped bases
                                    skip = -1;
                                } else {
                                    // If skipped bases are not called unmodified we are done.
                                    break;
                                }
                            }
                        } else {
                            if (skippedBasesCalled) {
                                // Skipped bases =>  "modification present with 0% probability"
                                for (String modification : modifications) {
                                    byte likelihood = 0;
                                    likelihoodMap.get(modification).put(position, likelihood);
                                }
                            }
                            matchCount++;
                        }
                    }
                    p++;
                }

                for (String m : modifications) {
                    modificationSets.add(new BaseModificationSet(base, strand, m, likelihoodMap.get(m)));
                }
            }
        }

        return modificationSets;
    }

    /**
     * If a string can be converted to a positive integer assume its a ChEBI ID
     *
     * @param str
     * @return
     */
    public static boolean isChEBI(String str) {
        if (str == null) {
            return false;
        }
        int length = str.length();
        if (length == 0) {
            return false;
        }
        for (int i = 0; i < length; i++) {
            char c = str.charAt(i);
            if (c < '0' || c > '9') {
                return false;
            }
        }
        return true;
    }

    /**
     * Minimally validate an MM tag.  This will not catch all problems, but will many.  Validation proceeds as follows
     * 1. Validate types of MM and ML tags.  This catches missues of the tags, for example in certain 10X files.
     * 2. If available, validate sequence length vs MN tag.
     * 3. If MN tag is not available, validate implied minimum count of base nucleotide vs actual count.
     *
     * @return
     */
    public static boolean validateMMTag(String readName, String mm, byte[] sequence, boolean isNegativeStrand) {


        // Finally, test implied minimum base count vs actual base count in sequence.  The minimum base count is
        // equal to the number of modified bases + the number of skipped bases as codified in the MM tag.
        // e.g. C+m,5,12,0   => at least 20 "Cs" in the read sequence, 3 with modifications and 17 skipped
        if (PreferencesManager.getPreferences().getAsBoolean(Constants.BASEMOD_VALIDATE_BASE_COUNT)) {

            String[] mmTokens = mm.split(";");

            for (String mmi : mmTokens) {
                String[] tokens = mmi.split(","); //Globals.commaPattern.split(mm);
                int baseCount;
                if (tokens[0].charAt(0) == 'N') {
                    baseCount = sequence.length;
                } else {
                    byte base = (byte) tokens[0].charAt(0);  // "Top strand" base seen by sequencing instrument
                    byte readBase = isNegativeStrand ? SequenceUtil.complement(base) : base;  // Base as reported in BAM file
                    baseCount = 0;
                    for (int i = 0; i < sequence.length; i++) {
                        if (readBase == sequence[i]) baseCount++;
                    }
                }

                // Count # of bases implied by tag
                int modified = tokens.length - 1;    // All tokens but the first are "skip" numbers
                int skipped = 0;
                for (int i = 1; i < tokens.length; i++) {
                    skipped += Integer.parseInt(tokens[i]);
                }
                if (modified + skipped > baseCount) {
                    if (++MM_WARNING_COUNT < 21) {
                        log.warn(readName + "  MM base count validation failed: expected " + (modified + skipped) + "'" + (tokens[0].charAt(0) + "'s" + ", actual count = " + baseCount));
                        if (MM_WARNING_COUNT == 20) {
                            log.warn("MM validation warning count exceeded.  Further failures will not be logged.");
                        }
                    }
                    return false;
                }
            }
        }

        // If we get here assume the tag is valid
        return true;
    }
}