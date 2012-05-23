/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.RuntimeUtils;
import org.broad.tribble.Feature;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * A feature source which combines results from other feature sources.
 * Currently uses bedtools to combine results
 * User: jacob
 * Date: 2012/05/01
 */
public class CombinedFeatureSource implements FeatureSource {

    private static Logger log = Logger.getLogger(CombinedFeatureSource.class);

    private FeatureSource sourceA;
    private FeatureSource sourceB;
    private Operation operation;

    private int featureWindowSize = -1;

    /**
     * Checks the global bedtools path, to see if bedtools
     * is actually there. Check is 2-fold:
     * First, we check if path exists.
     * If so, we run version command
     *
     * @return
     */
    public static boolean checkBEDToolsPathValid() throws IOException {
        String path = Globals.BEDtoolsPath;
        File bedtoolsFile = new File(path);
        boolean pathValid = bedtoolsFile.isFile();
        if (pathValid && !bedtoolsFile.canExecute()) {
            log.debug(path + " exists but is not executable. ");
            return false;
        }

        String cmd = path + " --version";
        String resp = RuntimeUtils.executeShellCommand(cmd, null, null);
        String line0 = resp.split("\n")[0].toLowerCase();
        pathValid &= line0.contains("bedtools v");
        pathValid &= !line0.contains("command not found");
        return pathValid;
    }

    /**
     * If known, it is recommended that sourceA be the larger of the two. sourceB will
     * be loaded into memory by BEDTools.
     *
     * @param sourceA
     * @param sourceB
     * @param operation How the two sources will be combined
     */
    public CombinedFeatureSource(FeatureSource sourceA, FeatureSource sourceB, Operation operation) {
        this.sourceA = sourceA;
        this.sourceB = sourceB;
        this.operation = operation;
    }

    /**
     * Stream will be closed after data written
     *
     * @param features
     * @param outputStream
     * @return
     */
    private int writeFeaturesToStream(Iterator<Feature> features, OutputStream outputStream) {
        PrintWriter writer = new PrintWriter(new OutputStreamWriter(outputStream));

        int numLines = 0;
        if (features != null) {
            IGVBEDCodec codec = new IGVBEDCodec();
            while (features.hasNext()) {
                String data = codec.encode(features.next());
                writer.println(data);
                numLines++;
            }
        }
        writer.flush();
        writer.close();

        return numLines;
    }

    /**
     * Perform the actual combination operation between the constituent data
     * sources. This implementation re-runs the operation each call.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws IOException
     */
    @Override
    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {
        Iterator<Feature> iterA = sourceA.getFeatures(chr, start, end);
        Iterator<Feature> iterB = sourceB.getFeatures(chr, start, end);

        //First we need to write out the features in A&B to a file
        //There are size limits on pipes so we don't use stdin
        File outFileA = File.createTempFile("featuresA", ".bed", null);
        File outFileB = File.createTempFile("featuresB", ".bed", null);
        outFileA.deleteOnExit();
        outFileB.deleteOnExit();

        int numA = writeFeaturesToStream(iterA, new FileOutputStream(outFileA));
        int numB = writeFeaturesToStream(iterB, new FileOutputStream(outFileB));

        //Start bedtools process
        String cmd = Globals.BEDtoolsPath + " " + this.operation.getCmd() + " -b " + outFileB.getAbsolutePath() + " -a " + outFileA.getAbsolutePath();
        Process pr = RuntimeUtils.startExternalProcess(cmd, null, null);

        //This is un-necessary I believe, and occasionally will hang
//        try {
//            pr.waitFor();
//        } catch (InterruptedException e) {
//            e.printStackTrace();
//            throw new IOException(e);
//        }

        //Read back in the data which bedtools output
        BufferedReader in = new BufferedReader(new InputStreamReader(pr.getInputStream()));
        BufferedReader err = new BufferedReader(new InputStreamReader(pr.getErrorStream()));

        List<Feature> featuresList = new ArrayList<Feature>(numA + numB);
        IGVBEDCodec codec = new IGVBEDCodec();


        String line;
        Feature feat;
        while ((line = in.readLine()) != null) {
            if (operation == Operation.WINDOW || operation == Operation.CLOSEST) {
                String[] closest = splitDualFeatures(line, 3)[1];
                //If not found, bedtools returns -1 for positions
                if (closest[1].trim().equalsIgnoreCase("-1")) {
                    continue;
                }
                feat = codec.decode(closest);
            } else {
                feat = codec.decode(line);
            }
            featuresList.add(feat);
        }

        in.close();


        while ((line = err.readLine()) != null) {
            log.error(line);
        }
        err.close();

        return featuresList.iterator();
    }

    /**
     * Certain bedtools commands output features side by side
     * e.g.
     * chr1 5   10  chr1    8   20
     * <p/>
     * might be one line, the first 3 columns representing data from file A
     * and the second 3 representing data from file B
     *
     * @param input
     * @param colsPerFeat Number of columns each feature should have.
     *                    input must have >= 2x this number. Extra columns
     *                    are ignored.
     * @return A 2-D string array. First index is length 2, second is the number of
     *         columns each feature has. out[0] is the first feature, out[1] is the second
     */
    private String[][] splitDualFeatures(String input, int colsPerFeat) {
        String[] tokens = Globals.singleTabMultiSpacePattern.split(input);

        assert tokens.length == colsPerFeat * 2;

        String[] feat1 = new String[colsPerFeat];
        String[] feat2 = new String[colsPerFeat];
        for (int cc = 0; cc < colsPerFeat; cc++) {
            feat1[cc] = tokens[cc];
            feat2[cc] = tokens[cc + colsPerFeat];
        }

        String[][] out = new String[][]{feat1, feat2};
        return out;

    }

    /**
     * Bedtools reports certain features as:
     * chr  start   end some_number
     * Where some_number might be coverage, overlap, fraction, etc
     *
     * @param input
     * @return
     */
    private String[] convertBedToolsOutToBed(String[] input) {
        if (input.length < 3) {
            throw new IllegalArgumentException("Input array has only " + input.length + " columns, need at least 3");
        }

        //No score data
        if (input.length == 3) {
            return input;
        }

        String[] output = new String[input.length + 1];
        System.arraycopy(input, 0, output, 0, 3);
        output[3] = "";
        System.arraycopy(input, 3, output, 4, input.length - 3);
        return output;
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return null;
    }

    /**
     * If this track has not had it's feature window size set,
     * we use the minimum of the sources
     *
     * @return
     */
    @Override
    public int getFeatureWindowSize() {
        if (featureWindowSize < 0) {
            return Math.min(sourceA.getFeatureWindowSize(), sourceB.getFeatureWindowSize());
        } else {
            return featureWindowSize;
        }
    }

    @Override
    public void setFeatureWindowSize(int size) {
        featureWindowSize = size;
    }

    public enum Operation {
        //We use these bed flags to ensure output will be in bed format, even
        //if input is bam
        //TODO Include -wo, -wb options
        INTERSECT("intersect -bed"),
        SUBTRACT("subtract"),
        //Identify the "closest" feature in file B for each feature in file A
        CLOSEST("closest"),
        //TODO include -d option
        WINDOW("window -bed"),
        COVERAGE("coverage");


        private String cmd;

        private Operation(String cmd) {
            this.cmd = cmd;
        }

        public String getCmd() {
            return cmd;
        }
    }
}
