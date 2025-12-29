package org.broad.igv.util;


import htsjdk.tribble.util.LittleEndianOutputStream;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;

import static junit.framework.Assert.assertEquals;

import org.junit.Ignore;


/**
 * @author jrobinso
 *         Date: 7/19/12
 *         Time: 8:54 AM
 */
@Ignore
public class LittleEndianInputStreamTest {

    static String testString = "987043\t1004564\tIDUA MONDO:0001586\t0\t+\t987043\t1004564\t39,103,73\tENST00000514224.2\tENSG00000127415.13\tNM_000203.5\tGENCC_000108-HGNC_5391-MONDO_0001586-HP_0000007-GENCC_100001\tHGNC:5391\tIDUA\tMONDO:0001586\tmucopolysaccharidosis type 1\tMONDO:0001586\tmucopolysaccharidosis type 1\tGENCC:100001\tDefinitive\tHP:0000007\tAutosomal recessive\tGENCC:000108\tMyriad Women’s Health\tHGNC:5391\tIDUA\tMONDO:0001586\tMucopolysaccharidosis type 1\tHP:0000007\tAutosomal recessive inheritance\tGENCC:000108\tMyriad Womens Health\tGENCC:100001\tDefinitive\t2018-06-08 20:05:09\thttps://onlinelibrary.wiley.com/doi/full/10.1002/humu.24033\t\t\thttps://clinicalgenome.org/docs/gene-disease-validity-sop-version-5/\t84\t2020-10-13\n";

    @BeforeClass
    public static void setup() throws Exception {
        createTestFile();
    }

    @Test
    public void testRead() throws Exception {
        LittleEndianInputStream lis = new LittleEndianInputStream(new BufferedInputStream(new FileInputStream("les_test.bin")));

        String str = lis.readString().strip();

        for(int i=0; i<testString.strip().length(); i++) {
//            System.out.println("" + testString.charAt(i) + '\t' + str.charAt(i));
            if(str.charAt(i) != testString.charAt(i)) {
                System.out.println("" + i + '\t' + testString.charAt(i) + '\t' + str.charAt(i) + '\t' + testString.getBytes()[i] + "\t" + str.getBytes()[i]);
            }
            assertEquals(str.charAt(i), testString.charAt(i));
        }

        //assertEquals(testString.strip(), lis.readString().strip());
        assertEquals(Float.MAX_VALUE, lis.readFloat());
        assertEquals(Byte.MAX_VALUE, lis.readByte());
        assertEquals(Short.MAX_VALUE, lis.readShort());
        assertEquals(Integer.MAX_VALUE, lis.readInt());
        assertEquals(Long.MAX_VALUE, lis.readLong());
        assertEquals(Double.MAX_VALUE, lis.readDouble());
        lis.close();
    }


    static void createTestFile() throws IOException {

        FileOutputStream fos = new FileOutputStream("les_test.bin");
        fos.write(testString.getBytes());
        fos.write((byte) 0);
        fos.flush();

        LittleEndianOutputStream los = new LittleEndianOutputStream(new BufferedOutputStream(fos));
        los.writeFloat(Float.MAX_VALUE);
        los.writeByte(Byte.MAX_VALUE);
        los.writeShort(Short.MAX_VALUE);
        los.writeInt(Integer.MAX_VALUE);
        los.writeLong(Long.MAX_VALUE);
        los.writeDouble(Double.MAX_VALUE);

        los.close();


    }
}
