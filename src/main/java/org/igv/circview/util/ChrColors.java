package org.igv.circview.util;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Default per-chromosome colors. Port of chrColor.js.
 *
 * <p>Lookups try the bare name, then a "chr"-prefixed name, and otherwise assign
 * (and cache) a random color, matching getChrColor() in the JS source. Two source
 * entries (chrY, chrUn) were missing a closing paren in the original RGB strings;
 * they are corrected here.
 */
public final class ChrColors {

    private static final Map<String, Color> CHR_COLOR_MAP = new HashMap<>();
    private static final Random RANDOM = new Random();

    private ChrColors() {
    }

    /**
     * Color for a chromosome name. Returns a stable color for known names and a
     * cached random color for unknown ones.
     */
    public static synchronized Color getChrColor(String chr) {
        Color c = CHR_COLOR_MAP.get(chr);
        if (c != null) {
            return c;
        }
        Color prefixed = CHR_COLOR_MAP.get("chr" + chr);
        if (prefixed != null) {
            CHR_COLOR_MAP.put(chr, prefixed);
            return prefixed;
        }
        Color random = randomRGB();
        CHR_COLOR_MAP.put(chr, random);
        return random;
    }

    private static Color randomRGB() {
        return new Color(RANDOM.nextInt(256), RANDOM.nextInt(256), RANDOM.nextInt(256));
    }

    static {
        put("chrX", 204, 153, 0);
        put("chrY", 153, 204, 0);   // source bug: missing ')' in chrColor.js
        put("chrUn", 50, 50, 50);
        put("chr1", 80, 80, 255);
        put("chrI", 139, 155, 187);
        put("chr2", 206, 61, 50);
        put("chrII", 206, 61, 50);
        put("chr2a", 216, 71, 60);
        put("chr2b", 226, 81, 70);
        put("chr3", 116, 155, 88);
        put("chrIII", 116, 155, 88);
        put("chr4", 240, 230, 133);
        put("chrIV", 240, 230, 133);
        put("chr5", 70, 105, 131);
        put("chr6", 186, 99, 56);
        put("chr7", 93, 177, 221);
        put("chr8", 128, 34, 104);
        put("chr9", 107, 215, 107);
        put("chr10", 213, 149, 167);
        put("chr11", 146, 72, 34);
        put("chr12", 131, 123, 141);
        put("chr13", 199, 81, 39);
        put("chr14", 213, 143, 92);
        put("chr15", 122, 101, 165);
        put("chr16", 228, 175, 105);
        put("chr17", 59, 27, 83);
        put("chr18", 205, 222, 183);
        put("chr19", 97, 42, 121);
        put("chr20", 174, 31, 99);
        put("chr21", 231, 199, 111);
        put("chr22", 90, 101, 94);
        put("chr23", 204, 153, 0);
        put("chr24", 153, 204, 0);
        put("chr25", 51, 204, 0);
        put("chr26", 0, 204, 51);
        put("chr27", 0, 204, 153);
        put("chr28", 0, 153, 204);
        put("chr29", 10, 71, 255);
        put("chr30", 71, 117, 255);
        put("chr31", 255, 194, 10);
        put("chr32", 255, 209, 71);
        put("chr33", 153, 0, 51);
        put("chr34", 153, 26, 0);
        put("chr35", 153, 102, 0);
        put("chr36", 128, 153, 0);
        put("chr37", 51, 153, 0);
        put("chr38", 0, 153, 26);
        put("chr39", 0, 153, 102);
        put("chr40", 0, 128, 153);
        put("chr41", 0, 51, 153);
        put("chr42", 26, 0, 153);
        put("chr43", 102, 0, 153);
        put("chr44", 153, 0, 128);
        put("chr45", 214, 0, 71);
        put("chr46", 255, 20, 99);
        put("chr47", 0, 214, 143);
        put("chr48", 20, 255, 177);
    }

    private static void put(String name, int r, int g, int b) {
        CHR_COLOR_MAP.put(name, new Color(r, g, b));
    }
}
