package org.broad.igv.sam;

import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.SAM_BASE_QUALITY_MAX;
import static org.broad.igv.prefs.Constants.SAM_BASE_QUALITY_MIN;

public class BaseModifications {

    public char base;
    public char strand;
    public String modification;
    public int[] positions;

    public BaseModifications(char base, char strand, String modification, int[] positions) {
        this.base = base;
        this.strand = strand;
        this.modification = modification;
        this.positions = positions;
    }

    public static Map<Integer, Mod> getBaseModificationMap(String mm, byte[] ml, byte[] sequence, boolean isNegativeStrand) {
        Map<Integer, Mod> map = new HashMap<>();
        List<BaseModifications> mods = BaseModifications.getBaseModifications(mm, sequence, isNegativeStrand);
        for (BaseModifications m : mods) {
            int[] positions = m.positions;
            for (int i = 0; i < positions.length; i++) {
                Integer pos = positions[i];
                byte likelihood = ml == null ? (byte) 255 : ml[i];
                if (!map.containsKey(pos) || (map.containsKey(pos) && Byte.toUnsignedInt(likelihood) > Byte.toUnsignedInt(map.get(pos).likelihood))) {
                    map.put(positions[i], new BaseModifications.Mod(m.base, m.strand, likelihood));
                }
            }
        }
        return map;
    }

    public static List<BaseModifications> getBaseModifications(String mm, byte[] sequence, boolean isNegativeStrand) {

        if (isNegativeStrand) {
            sequence = AlignmentUtils.reverseComplementCopy(sequence);
        }

        List<BaseModifications> mods = new ArrayList<>();

        String[] mmTokens = mm.split(";");

        for (String mmi : mmTokens) {
            String[] tokens = mmi.split(","); //Globals.commaPattern.split(mm);
            if (tokens.length > 3) {
                char base = tokens[0].charAt(0);
                char strand = tokens[0].charAt(1);
                String modification = tokens[0].substring(2);
                int[] positions = new int[tokens.length - 1];
                int idx = 0;
                int s = 0;
                int skip = Integer.parseInt(tokens[idx + 1]);
                int matchCount = 0;
                while (idx < positions.length && s < sequence.length) {
                    if (sequence[s] == base) {
                        if (matchCount == skip) {
                            positions[idx] = s;
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

                if (isNegativeStrand) {
                    int[] reversedPositions = new int[positions.length];
                    for (int i = 0; i < positions.length; i++) {
                        reversedPositions[i] = sequence.length - 1 - positions[i];
                    }
                    mods.add(new BaseModifications(base, strand, modification, reversedPositions));
                } else {
                    mods.add(new BaseModifications(base, strand, modification, positions));
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
