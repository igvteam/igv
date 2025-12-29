package org.broad.igv.feature.tribble;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.bedpe.InteractFeature;
import org.broad.igv.feature.FeatureType;
import org.broad.igv.bedpe.InteractCodec;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 * Date: 4/8/13
 * Time: 2:51 PM
 */
public class InteractCodecTest extends AbstractHeadlessTest {


    /**
     * #chrom  chromStart  chromEnd  name  score  value  exp  color  sourceChrom  sourceStart  sourceEnd  sourceName  sourceStrand  targetChrom  targetStart  targetEnd  targetName  targetStrand
     * chr12    40572709      40618813        rs7974522/LRRK2/muscleSkeletal  0       0.624   muscleSkeletal  #7A67EE  chr12 40572709        40572710  rs7974522       .       chr12   40618812        40618813  LRRK2   +
     * chr12    40579899      40618813        rs17461492/LRRK2/muscleSkeletal  0       0.624   muscleSkeletal  #7A67EE chr12 40579899        40579900  rs17461492       .       chr12   40618812        40618813  LRRK2   +
     * chr12    40614433      40618813        rs76904798/LRRK2/nerveTibial  0       0.625   nerveTibial  #FFD700       chr12 40614433        40614434  rs76904798       .       chr12   40618812        40618813  LRRK2   +
     * chr12    40618812      40652520        rs2723264/LRRK2/lung  0       1.839   lung     #9ACD32  chr12      40652519        40652520  rs2723264       .       chr12   40618812        40618813  LRRK2   +
     */
    @Test
    public void testInteract1() throws Exception {

        String path = TestUtils.DATA_DIR + "bedpe/interactExample1.txt";
        InteractCodec codec = new InteractCodec(genome, FeatureType.INTERACT);
        BufferedReader bufferedReader = null;
        try {
            bufferedReader = new BufferedReader(new FileReader(path));
            String nextLine;
            ArrayList<InteractFeature> features = new ArrayList<>();
            while ((nextLine = bufferedReader.readLine()) != null) {
                InteractFeature f = codec.decode(nextLine);
                if (f != null) {
                    assertEquals("chr12", f.getChr());
                    features.add(f);
                }
            }
            assertEquals(4, features.size());
        } finally {
            bufferedReader.close();
        }
    }



}
