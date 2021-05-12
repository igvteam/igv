package org.broad.igv.sam;

import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.*;
import java.util.List;

public class BaseModification {

    public String modification;
    char strand;
    public int position;
    public byte likelihood;

    static Set<Character> STANDARD_CODES = new HashSet<>(Arrays.asList('C', 'm', 'h', 'f', 'c', 'T', 'g', 'e', 'b', 'U', 'A', 'a', 'G', 'o', 'N', 'n'));

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
    }

    public BaseModification(String modification, char strand, int position, byte likelihood) {
        this.likelihood = likelihood;
        this.modification = modification;
        this.position = position;
    }

    public String valueString() {
        int l = (int) (100.0 * Byte.toUnsignedInt(likelihood) / 255);
        return "Base modification: " +
                ((codeValues.containsKey(modification)) ? codeValues.get(modification) : "Uknown") + " (" + l + "%)";
    }

    public static List<BaseModification> getBaseModifications(String mm, byte[] ml, byte[] sequence, boolean isNegativeStrand) {
        byte[] origSequence = sequence;
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
                    // If string can be converted to a number assume its a ChEBI code
                    boolean isChEBI = true;
                    try {
                        Integer.parseInt(modificationString);
                    } catch (NumberFormatException e) {
                        isChEBI = false;
                    }
                    if (isChEBI) {
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


                int idx = 0;
                int s = 0;

                int[] positions = new int[tokens.length - 1];
                int skip = Integer.parseInt(tokens[idx + 1]);
                int matchCount = 0;

                while (idx < positions.length) {
                    if (s >= sequence.length) {
                        System.err.println("Ran out of sequence");
                        System.out.println(isNegativeStrand);
                        System.out.println(mm);
                        System.out.println(new String(origSequence));
                        System.out.println(new String(sequence));
                        break;
                    }
                    if (base == 'N' || sequence[s] == base) {
                        if (matchCount == skip) {
                            int position = isNegativeStrand ? sequence.length - 1 - s : s;
                            for (String modification : modifications) {
                                byte likelihood = ml == null ? (byte) 255 : ml[mlIdx];
                                mods.add(new BaseModification(modification, strand, position, likelihood));
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
        Color baseColor;

        if (modification.equals("m")) {
            baseColor = Color.red;
        } else {
            if (modColorPallete == null) {
                modColorPallete = new PaletteColorTable(ColorUtilities.getPalette("Set 1"));
            }
            baseColor = modColorPallete.get(modification);
        }

        int l = Byte.toUnsignedInt(likelihood);
        if (l > 250) {
            return baseColor;
        }

        l = Math.max(25, l);
        String key = modification + "--" + l;
        if (!modColorMap.containsKey(key)) {
            modColorMap.put(key, new Color(baseColor.getRed(), baseColor.getGreen(), baseColor.getBlue(), l));
        }

        return modColorMap.get(key);
    }

    static PaletteColorTable modColorPallete;

    static Map<String, Color> modColorMap = new HashMap<>();

}
