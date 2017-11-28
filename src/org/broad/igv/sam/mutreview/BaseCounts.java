package org.broad.igv.sam.mutreview;


import org.broad.igv.sam.AlignmentCounts;

public class BaseCounts {

    public static char [] bases = {'A', 'C', 'T', 'G', 'N'};
    public int totalCount;
    public int [] negCounts;
    public int [] posCounts;
    public int delCount;
    public int insCount;

    public BaseCounts(AlignmentCounts ac, int pos) {

        negCounts = new int[5];
        posCounts = new int[5];

        for(int i=0; i<5; i++) {
            negCounts[i] = ac.getNegCount(pos, (byte) bases[i]);
            posCounts[i] = ac.getPosCount(pos, (byte) bases[i]);
        }
        delCount = ac.getDelCount(pos);
        insCount = ac.getInsCount(pos);
        totalCount = ac.getTotalCount(pos);
    }

}
