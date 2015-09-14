/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.igv.util.StringUtils;
import htsjdk.tribble.Feature;

import java.io.*;
import java.util.*;

/**
 * A feature source which combines results from other feature sources.
 * Currently uses bedtools to combine results
 * User: jacob
 * Date: 2012/05/01
 * @deprecated Use {@link org.broad.igv.cli_plugin.PluginFeatureSource}
 */
@Deprecated
public class CombinedFeatureSource implements FeatureSource {

    private static Logger log = Logger.getLogger(CombinedFeatureSource.class);

    private FeatureSource[] sources;
    private Operation operation;

    /**
     * Checks the global bedtools path, to see if bedtools
     * is actually there. Check is 2-fold:
     * First, we check if path exists.
     * If so, we run version command
     *
     * @return
     */
    public static boolean checkBEDToolsPathValid() {
        String path = FileUtils.findExecutableOnPath(Globals.BEDtoolsPath);
        File bedtoolsFile = new File(path);
        boolean pathValid = bedtoolsFile.isFile();
        if (pathValid && !bedtoolsFile.canExecute()) {
            log.debug(path + " exists but is not executable. ");
            return false;
        }

        String[] cmd = new String[]{path, "--version"};
        String resp;
        try {
            resp = RuntimeUtils.executeShellCommand(cmd, null, null);
        } catch (IOException e) {
            log.error(e);
            return false;
        }
        String line0 = resp.split("\n")[0].toLowerCase();
        pathValid &= line0.contains("bedtools v");
        pathValid &= !line0.contains("command not found");
        return pathValid;
    }

    public CombinedFeatureSource(Collection<Track> tracks, Operation operation) {
        List<FeatureSource> sources = new ArrayList<FeatureSource>(tracks.size());
        for (Track t : tracks) {
            if (t instanceof FeatureTrack) {
                sources.add(((FeatureTrack) t).source);
            }
        }

        init(sources.toArray(new FeatureSource[0]), operation);
    }

    public CombinedFeatureSource(FeatureSource[] sources, Operation operation) {
        init(sources, operation);
    }

    /**
     * If known, it is recommended that source[0] be the larger of the two. sources[1] will
     * be loaded into memory by BEDTools.
     *
     * @param sources
     * @param operation How the two sources will be combined
     */
    private void init(FeatureSource[] sources, Operation operation) {
        this.sources = sources;
        this.operation = operation;
        if (sources.length != 2 && operation != Operation.MULTIINTER) {
            throw new IllegalArgumentException("sources must be length 2 for operation " + operation);
        }
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

        int allNumCols = -1;
        if (features != null) {
            IGVBEDCodec codec = new IGVBEDCodec();
            while (features.hasNext()) {
                String data = codec.encode(features.next());

                writer.println(data);

                //We require consistency of output
                int tmpNumCols = data.split("\t").length;
                if (allNumCols < 0) {
                    allNumCols = tmpNumCols;
                } else {
                    assert tmpNumCols == allNumCols;
                }
            }
        }
        writer.flush();
        writer.close();

        return allNumCols;
    }

    /**
     * Write out data from feature sources within the specified interval
     * to temporary files.
     *
     * @param chr
     * @param start
     * @param end
     * @return LinkedHashMap from TempFileName -> number of columns in data file
     *         A LinkedHashMap has a predictable iteration order, which will be the same as the
     *         insertion order, which will be the order of sources
     * @throws IOException
     */
    private LinkedHashMap<String, Integer> createTempFiles(String chr, int start, int end) throws IOException {
        LinkedHashMap<String, Integer> tempFiles = new LinkedHashMap<String, Integer>(sources.length);
        for (FeatureSource source : sources) {
            Iterator<Feature> iter = source.getFeatures(chr, start, end);

            File outFile = File.createTempFile("features", ".bed", null);
            outFile.deleteOnExit();

            int numCols = writeFeaturesToStream(iter, new FileOutputStream(outFile));
            tempFiles.put(outFile.getAbsolutePath(), numCols);
        }
        return tempFiles;
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

        String cmd = Globals.BEDtoolsPath + " " + this.operation.getCmd();
        LinkedHashMap<String, Integer> tempFiles = createTempFiles(chr, start, end);
        String[] fiNames = tempFiles.keySet().toArray(new String[0]);
        if (operation == Operation.MULTIINTER) {
            assert tempFiles.size() >= 2;
            cmd += " -i " + StringUtils.join(tempFiles.keySet().toArray(), " ");
        } else {
            assert tempFiles.size() == 2;
            cmd += " -a " + fiNames[0] + " -b " + fiNames[1];
        }

        //Start bedtools process
        Process pr = RuntimeUtils.startExternalProcess(new String[]{cmd}, null, null);

        //Read back in the data which bedtools output
        BufferedReader in = new BufferedReader(new InputStreamReader(pr.getInputStream()));

        List<Feature> featuresList = new ArrayList<Feature>();
        //TODO This cast is here as a reminder that we want to use AsciiFeatureCodec
        IGVBEDCodec codec = (IGVBEDCodec) CodecFactory.getCodec(".bed", null);

        String line;
        Feature feat;
        int numCols0 = tempFiles.get(fiNames[0]);
        int numCols1 = tempFiles.get(fiNames[1]);
        while ((line = in.readLine()) != null) {
            //System.out.println(line);
            String[] tokens = line.split("\t");
            if (operation.getCmd().contains("-split")) {
                //When we split, the returned feature still has the exons
                //We don't want to plot them all a zillion times
                tokens = Arrays.copyOfRange(tokens, 0, Math.min(6, tokens.length));
            }

            if (operation == Operation.WINDOW || operation == Operation.CLOSEST) {

                String[] closest = Arrays.copyOfRange(tokens, numCols0, numCols0 + numCols1);
                //If not found, bedtools returns -1 for positions
                if (closest[1].trim().equalsIgnoreCase("-1")) {
                    continue;
                }
                feat = codec.decode(closest);
            } else if (operation == Operation.MULTIINTER) {
                //We only look at regions common to ALL inputs
                //Columns: chr \t start \t \end \t # of files which contained this feature \t comma-separated list files +many more
                int numRegions = Integer.parseInt(tokens[3]);
                if (numRegions < sources.length) {
                    continue;
                }
                String[] intersection = Arrays.copyOf(tokens, 3);
                feat = codec.decode(intersection);
            } else {
                feat = codec.decode(tokens);
            }
            featuresList.add(feat);
        }

        in.close();

        return featuresList.iterator();
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
        int featureWindowSize = Integer.MAX_VALUE;
        for (FeatureSource source : sources) {
            featureWindowSize = Math.min(featureWindowSize, source.getFeatureWindowSize());
        }
        return featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        //no-op
    }

    public enum Operation {
        //We use these bed flags to ensure output will be in bed format, even
        //if input is bam
        //TODO Include -wo, -wb options
        INTERSECT("intersect -bed -split"),
        SUBTRACT("subtract"),
        //Identify the "closest" feature in file B for each feature in file A
        CLOSEST("closest"),
        //TODO include -d option
        WINDOW("window -bed"),
        COVERAGE("coverage -split"),
        MULTIINTER("multiinter");


        private String cmd;

        private Operation(String cmd) {
            this.cmd = cmd;
        }

        public String getCmd() {
            return cmd;
        }

    }
}
