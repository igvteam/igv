package org.broad.igv.hic.data;

import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 11, 2010
 */
public class ChromSizes {

    static Map<String, ChromSizes> genomeCache = new HashMap();
    String genomeId;
    int[] sizes;

    public static ChromSizes getChromSizes(String id) {
        return genomeCache.get(id);
    }

    public ChromSizes(String genomeId, int[] sizes) {
        this.genomeId = genomeId;
        this.sizes = sizes;
    }


    public int getSize(int chr) {
        return sizes[chr];
    }


    static {
        int[] sizes = new int[25];
        sizes[1] = 247249719;
        sizes[2] = 242951149;
        sizes[3] = 199501827;
        sizes[4] = 191273063;
        sizes[5] = 180857866;
        sizes[6] = 170899992;
        sizes[7] = 158821424;
        sizes[8] = 146274826;
        sizes[9] = 140273252;
        sizes[10] = 135374737;
        sizes[11] = 134452384;
        sizes[12] = 132349534;
        sizes[13] = 114142980;
        sizes[14] = 106368585;
        sizes[15] = 100338915;
        sizes[16] = 88827254;
        sizes[17] = 78774742;
        sizes[18] = 76117153;
        sizes[19] = 63811651;
        sizes[20] = 62435964;
        sizes[21] = 46944323;
        sizes[22] = 49691432;
        sizes[23] = 154913754;
        sizes[24] = 57772954;
        sizes[0] = 16571;

        ChromSizes hg18 = new ChromSizes("hg18", sizes);
        genomeCache.put("hg18", hg18);

        int [] sizes2 = new int[8];
        sizes2[1] = 23011544;  // 2L
        sizes2[2] = 21146708;  // 2R
        sizes2[3] = 24543557; // 3L
        sizes2[4] = 27905053; //3R
        sizes2[5] = 1351857; // 4
        sizes2[6] = 22422827; // X
        sizes2[7] = 10049037; // U
        genomeCache.put("dmel", new ChromSizes("dmel", sizes2));
    }
}
