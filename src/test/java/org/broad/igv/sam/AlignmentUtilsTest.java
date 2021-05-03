package org.broad.igv.sam;

import org.junit.Test;

import static org.junit.Assert.*;

public class AlignmentUtilsTest {

    @Test
    public void getModificationPositions() {

        //top-fwd	0	*	0	0	*	*	0	0	AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT	*	Mm:Z:C+m,1,3,0;	Ml:B:C,128,153,179
        byte [] sequence = "AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT".getBytes();
        String mm = "C+m,1,3,0";

        int [] expectedPositions = {7, 30, 31};
        int [] positions = AlignmentUtils.getModificationPositions(mm, sequence);
        assertArrayEquals(expectedPositions, positions);

    }
}