package org.broad.igv.hic.tools;

import org.broad.igv.feature.Chromosome;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Utility class for converting an asscii "pair" file to a compact binary format.  Can greatly speed up calculation
 * of "hic" file
 * Input:   D0J8AACXX120130:6:1101:1003:8700/1 15 61559113 0 D0J8AACXX120130:6:1101:1003:8700/2 15 61559309 16
 * Output:  [chr1 index][pos1][chr2 index][pos 2]    (int, int, int, int)/
 *
 * @author Jim Robinson
 * @date 4/7/12
 */
public class AsciiToBinConverter {

    /**
     * @param inputPath
     * @param outputFile
     */
    public static void convert(String inputPath, String outputFile, List<Chromosome> chromosomes) throws IOException {

        Map<String, Integer> chromosomeOrdinals = new HashMap<String, Integer>();
        for(Chromosome c : chromosomes) {
            chromosomeOrdinals.put(c.getName(), c.getIndex());
        }

        AsciiPairIterator iter = null;
        BufferedOutputStream bos = null;
        try {
            bos = new BufferedOutputStream(new FileOutputStream(outputFile));
            LittleEndianOutputStream les = new LittleEndianOutputStream(bos);
            iter = new AsciiPairIterator(inputPath, chromosomeOrdinals);
            while (iter.hasNext()) {
                AlignmentPair pair = iter.next();
                les.writeInt(pair.getChr1());
                les.writeInt(pair.getPos1());
                les.writeInt(pair.getChr2());
                les.writeInt(pair.getPos2());
            }
            les.flush();
            bos.flush();
        } finally {
            if (iter != null) iter.close();
            if(bos != null) bos.close();

        }


    }
    
    public static void convertBack(String inputPath, String outputFile) throws IOException {
        
        PrintWriter pw = null;
        try {
            File f = new File(outputFile);
            FileWriter fw = new FileWriter (f);
            pw = new PrintWriter (fw);
            Map<String, Integer> chromosomeIndexMap = null;  // TODO
            BinPairIterator iter = new BinPairIterator(inputPath, chromosomeIndexMap);
            while (iter.hasNext()) {
                AlignmentPair pair = iter.next();
                pw.println("R1\t" + pair.getChr1() + "\t" + pair.getPos1() + "\tS1\tR2\t" +pair.getChr2() + "\t" + pair.getPos2() +"\t" + "S2");

            }
        }
        finally {
            pw.close();
        }
    }
}
