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

package org.broad.igv.dev.plugin;

import org.apache.log4j.Logger;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;

import java.io.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;

/**
 * A feature source which derives its information
 * from a command line plugin
 * User: jacob
 * Date: 2012/05/01
 */
public class PluginFeatureSource implements FeatureSource {

    private static Logger log = Logger.getLogger(PluginFeatureSource.class);

    /**
     * Initial command tokens. This will typically include both
     * the executable and command, e.g. {"/usr/bin/bedtools", "intersect"}
     */
    private final List<String> cmd;
    private final LinkedHashMap<Argument, Object> arguments;

    protected boolean strictParsing = true;
    protected String decodingCodec = null;
    protected URL[] decodingLibURLs = new URL[0];
    protected String format = "bed";

    /**
     * For decoding, we may need to know how many columns
     * were written out in the first place
     */
    protected Map<String, Integer> outputCols = new LinkedHashMap<String, Integer>(2);

    public PluginFeatureSource(List<String> cmd, LinkedHashMap<Argument, Object> arguments, Map<String, String> parsingAttrs, String specPath) {
        this.cmd = cmd;
        this.arguments = arguments;

        setParsingAttributes(parsingAttrs, specPath);
    }

    private void setParsingAttributes(Map<String, String> parsingAttrs, String specPath) {
        this.decodingCodec = parsingAttrs.get("decoding_codec");

        String tmpStrict = parsingAttrs.get("strict");
        String fmt = parsingAttrs.get("format");
        this.format = fmt != null ? fmt : this.format;
        if (tmpStrict != null) this.strictParsing = Boolean.parseBoolean(tmpStrict);
        String libs = parsingAttrs.get("libs");
        libs = libs != null ? libs : "";
        decodingLibURLs = FileUtils.getURLsFromString(libs, specPath);
    }


    /**
     * Stream will be closed after data written
     *
     * @param features
     * @param outputStream
     * @return
     */
    private int writeFeaturesToStream(Iterator<Feature> features, OutputStream outputStream, Argument argument) {
        PrintWriter writer = new PrintWriter(new OutputStreamWriter(outputStream));

        int allNumCols = -1;
        if (features != null) {
            FeatureEncoder codec = getEncodingCodec(argument);
            while (features.hasNext()) {
                String line = codec.encode(features.next());
                if (line == null) continue;
                writer.println(line);

                //We require consistency of output
                int tmpNumCols = codec.getNumCols(line);
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

    private String[] genFullCommand(String chr, int start, int end) throws IOException {

        List<String> fullCmd = new ArrayList<String>(cmd);

        outputCols.clear();
        Map<String, String[]> argValsById = new HashMap<String, String[]>(arguments.size());
        for (Map.Entry<Argument, Object> entry : arguments.entrySet()) {
            Argument arg = entry.getKey();

            assert arg.isValidValue(entry.getValue());
            String[] sVal = null;
            String ts = null;
            switch (arg.getType()) {
                case LONGTEXT:
                case TEXT:
                    ts = (String) entry.getValue();
                    if (ts == null || ts.trim().length() == 0) {
                        continue;
                    }
                    sVal = new String[]{ts};
                    if (arg.getType() == Argument.InputType.TEXT) {
                        sVal = ts.split("\\s+");
                    }
                    break;
                case FEATURE_TRACK:
                    sVal = createTempFiles((FeatureTrack) entry.getValue(), arg, chr, start, end);
                    break;
                case MULTI_FEATURE_TRACK:
                    sVal = createTempFiles((List<FeatureTrack>) entry.getValue(), arg, chr, start, end);
                    break;
            }

            if (arg.getId() != null) {
                argValsById.put(arg.getId(), sVal);
            }

            if (arg.isOutput()) {
                String cmdArg = arg.getCmdArg();
                if (cmdArg.trim().length() > 0) {
                    for (String argId : argValsById.keySet()) {
                        cmdArg = cmdArg.replace("$" + argId, argValsById.get(argId)[0]);
                    }

                    fullCmd.add(cmdArg);
                }
                fullCmd.addAll(Arrays.asList(sVal));
            }
        }

        return fullCmd.toArray(new String[0]);
    }

    /**
     * Convenience for calling {@link #createTempFiles(org.broad.igv.track.FeatureTrack, Argument, String, int, int)};
     * on each track
     *
     * @param tracks
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws IOException
     */
    private String[] createTempFiles(List<FeatureTrack> tracks, Argument argument, String chr, int start, int end) throws IOException {
        String[] fileNames = new String[tracks.size()];
        int fi = 0;
        for (FeatureTrack track : tracks) {
            fileNames[fi++] = (createTempFiles(track, argument, chr, start, end))[0];
        }
        return fileNames;
    }

    /**
     * Write out data from feature sources within the specified interval
     * to temporary files.
     *
     * @param track
     * @param chr
     * @param start
     * @param end
     * @return Array of strings with temp file names.
     * @throws java.io.IOException
     */
    private String[] createTempFiles(FeatureTrack track, Argument argument, String chr, int start, int end) throws IOException {

        List<Feature> features = track.getFeatures(chr, start, end);

        File outFile = File.createTempFile("features", ".tmp", null);
        outFile.deleteOnExit();

        int numCols = writeFeaturesToStream(features.iterator(), new FileOutputStream(outFile), argument);
        String path = outFile.getAbsolutePath();
        outputCols.put(path, numCols);
        return new String[]{path};
    }

    /**
     * Perform the actual combination operation between the constituent data
     * sources. This implementation re-runs the operation each call.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws java.io.IOException
     */
    @Override
    public final Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

        String[] fullCmd = genFullCommand(chr, start, end);

        //Start plugin process
        Process pr = RuntimeUtils.startExternalProcess(fullCmd, null, null);

        //Read back in the data which plugin output
        BufferedReader in = new BufferedReader(new InputStreamReader(pr.getInputStream()));

        List<Feature> featuresList = new ArrayList<Feature>();

        FeatureDecoder codec = getDecodingCodec(decodingCodec, decodingLibURLs, cmd, arguments);
        codec.setOutputColumns(outputCols);

        String line;
        Feature feat;
        while ((line = in.readLine()) != null) {
            try {
                feat = codec.decode(line);
                if (feat != null)
                    featuresList.add(feat);
            } catch (Exception e) {
                if (strictParsing) {
                    throw new RuntimeException(e);
                }
            }
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
    public final int getFeatureWindowSize() {
        int featureWindowSize = getMinSizeFromTracks(arguments.values());
        return featureWindowSize;
    }

    private int getMinSizeFromTracks(Iterable tracks) {
        int featureWindowSize = Integer.MAX_VALUE;
        for (Object val : tracks) {
            int tmpSize = Integer.MAX_VALUE;
            if (val instanceof FeatureTrack) {
                FeatureTrack track = ((FeatureTrack) val);
                tmpSize = track.getFeatureWindowSize();
            } else if (val instanceof List) {
                featureWindowSize = getMinSizeFromTracks((List) val);
            }
            featureWindowSize = Math.min(featureWindowSize, tmpSize);
        }
        return featureWindowSize;
    }

    @Override
    public final void setFeatureWindowSize(int size) {
        //no-op
    }

    protected FeatureEncoder getEncodingCodec(Argument argument) {
        String encodingCodec = argument.getEncodingCodec();

        if (encodingCodec == null) return new IGVBEDCodec();

        URL[] libURLs = argument.getLibURLs();

        try {
            ClassLoader loader = URLClassLoader.newInstance(
                    libURLs,
                    getClass().getClassLoader()
            );
            Class clazz = loader.loadClass(encodingCodec);
            Object codec = clazz.newInstance();
            return (FeatureEncoder) codec;
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(e);
        } catch (InstantiationException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

    }

    protected FeatureDecoder getDecodingCodec(String decodingCodec, URL[] libURLs, List<String> cmd, Map<Argument, Object> argumentMap) {
        if (decodingCodec == null) return
                new DecoderWrapper(CodecFactory.getCodec("." + format, GenomeManager.getInstance().getCurrentGenome()));

        try {

            ClassLoader loader = URLClassLoader.newInstance(
                    libURLs,
                    getClass().getClassLoader()
            );

            Class clazz = loader.loadClass(decodingCodec);
            Constructor constructor = clazz.getConstructor(List.class, Map.class);
            return (FeatureDecoder) constructor.newInstance(cmd, argumentMap);

        } catch (ClassNotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(e);
        } catch (NoSuchMethodException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        } catch (InvocationTargetException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        } catch (InstantiationException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        } catch (IllegalAccessException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

    }

    private class DecoderWrapper implements FeatureDecoder {

        private AsciiFeatureCodec wrappedCodec;

        public DecoderWrapper(AsciiFeatureCodec wrappedCodec) {
            this.wrappedCodec = wrappedCodec;
        }

        @Override
        public Feature decode(String line) {
            return wrappedCodec.decode(line);
        }

        @Override
        public void setOutputColumns(Map<String, Integer> outputColumns) {
            //no-op
        }
    }


}
