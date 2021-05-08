package org.broad.igv.sam;

import java.awt.*;
import java.util.*;
import java.util.List;

public class BaseModification {

    public String modification;
    public int position;
    public byte likelihood;

    static Set<Character> STANDARD_CODES = new HashSet<>(Arrays.asList('C', 'm', 'h', 'f', 'c', 'T', 'g', 'e', 'b', 'U', 'A', 'a', 'G', 'o', 'N', 'n'));

    public BaseModification(String modification, int position, byte likelihood) {
        this.likelihood = likelihood;
        this.modification = modification;
        this.position = position;
    }

//    public static Map<Integer, Mod> getBaseModificationMap(String mm, byte[] ml, byte[] sequence, boolean isNegativeStrand) {
//        Map<Integer, Mod> map = new HashMap<>();
//        List<BaseModification> mods = BaseModification.getBaseModifications(mm, sequence, isNegativeStrand);
//        for (BaseModification m : mods) {
//            int[] positions = m.positions;
//            for (int i = 0; i < positions.length; i++) {
//                Integer pos = positions[i];
//                byte likelihood = ml == null ? (byte) 255 : ml[i];
//                if (!map.containsKey(pos) || (map.containsKey(pos) && Byte.toUnsignedInt(likelihood) > Byte.toUnsignedInt(map.get(pos).likelihood))) {
//                    map.put(positions[i], new BaseModification.Mod(m.base, m.strand, likelihood));
//                }
//            }
//        }
//        return map;
//    }

    public static List<BaseModification> getBaseModifications(String mm, byte[] ml, byte[] sequence, boolean isNegativeStrand) {

        if (isNegativeStrand) {
            sequence = AlignmentUtils.reverseComplementCopy(sequence);
        }

        List<BaseModification> mods = new ArrayList<>();

        String[] mmTokens = mm.split(";");

        int mlIdx = 0;

        for (String mmi : mmTokens) {

            String[] tokens = mmi.split(","); //Globals.commaPattern.split(mm);
            if (tokens.length > 1) {

                char strand = tokens[0].charAt(1);
                char base = tokens[0].charAt(0);

                String modificationString = tokens[0].substring(2);
                String[] modifications = null;
                if (modificationString.length() > 1) {
                    // If all characters are in the standard code set assume this is a multi-modification
                    boolean multiMod = true;
                    for (byte b : modificationString.getBytes()) {
                        if (!STANDARD_CODES.contains((char) b)) {
                            multiMod = false;
                            break;
                        }
                    }
                    if (multiMod) {
                        modifications = new String[modificationString.length()];
                        for (int i = 0; i < modificationString.length(); i++) {
                            modifications[i] = modificationString.substring(i, i + 1);
                        }
                    }
                }
                if (modifications == null) {
                    modifications = new String[]{modificationString};
                }


                int idx = 0;
                int s = 0;

                int[] positions = new int[tokens.length - 1];
                int skip = Integer.parseInt(tokens[idx + 1]);
                int matchCount = 0;

                while (idx < positions.length) {
                    if (s >= sequence.length) {
                        System.err.println("Ran out of sequence");
                        break;
                    }
                    if (base == 'N' || sequence[s] == base) {
                        if (matchCount == skip) {
                            int position = isNegativeStrand ? sequence.length - 1 - s : s;
                            for(String modification : modifications) {
                                byte likelihood = ml == null ? (byte) 255 : ml[mlIdx];
                                mods.add(new BaseModification(modification, position, likelihood));
                                mlIdx++;
                            }
                            if (idx + 1 == positions.length) {
                                break;
                            }
                            idx++;
                            skip = Integer.parseInt(tokens[idx + 1]);
                            matchCount = 0;
                        } else {
                            matchCount++;
                        }
                    }
                    s++;
                }
            }
        }

        return mods;
    }

    public static class Mod {
        char base;
        char strand;
        byte likelihood;

        public Mod(char base, char strand, byte likelihood) {
            this.base = base;
            this.strand = strand;
            this.likelihood = likelihood;
        }
    }

    public static Color getModColor(String modification, byte likelihood) {

        // TODO -- determine base color by modification
        Color baseColor = Color.red;

        int l = Byte.toUnsignedInt(likelihood);
        if (l > 250) {
            return baseColor;
        }

        l = Math.max(25, l);
        if (!modColorMap.containsKey(l)) {
            modColorMap.put(l, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), l));
        }

        return modColorMap.get(l);
    }


    static Map<Integer, Color> modColorMap = new HashMap<>();

}
