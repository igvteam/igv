/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.track;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.ChromAliasDefaults;
import org.broad.igv.feature.genome.Genome;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 */
public class TestUtils {

    public static void main(String[] args) throws IOException {
        createFeatureTestFile();
    }

    public static void createFeatureTestFile() throws IOException {
        String chr = "chr1";
        long chrLength = 560000;
        long featureLength = 1;

        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("data/featureTest.txt")));

        long featureStart = 0;
        long featureEnd = featureStart + featureLength;
        int n = 0;
        while (featureEnd <= chrLength) {
            pw.println(chr + "\t" + featureStart + "\t" + featureEnd + "\t" + ("feature_" + n));
            featureStart += 2 * featureLength;
            featureEnd = featureStart + featureLength;
        }

        pw.close();
    }

}
