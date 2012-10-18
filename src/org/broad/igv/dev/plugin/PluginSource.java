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
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;

import java.io.*;
import java.lang.reflect.Constructor;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;

/**
 * A feature source which derives its information
 * from a command line plugin
 * User: jacob
 * Date: 2012/05/01
 */
abstract class PluginSource<E extends Feature, D extends Feature> {

    private static Logger log = Logger.getLogger(PluginSource.class);

    /**
     * Initial command tokens. This will typically include both
     * the executable and command, e.g. {"/usr/bin/bedtools", "intersect"}
     */
    protected final List<String> commands;
    protected final LinkedHashMap<Argument, Object> arguments;

    protected boolean strictParsing = true;
    protected String decodingCodec = null;
    protected URL[] decodingLibURLs = new URL[0];
    protected String format = "bed";

    /**
     * For decoding, we may need to know how many columns
     * were written out in the first place
     */
    protected List<Map<String, Object>> attributes = new ArrayList<Map<String, Object>>(2);

    public PluginSource(List<String> commands, LinkedHashMap<Argument, Object> arguments, Map<String, String> parsingAttrs, String specPath) {
        this.commands = commands;
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
     * Encode features into strings using {@link #getEncodingCodec(Argument)} and write them to the provided stream.
     * Stream will be closed after data written
     *
     * @param outputStream
     * @param features
     * @param argument
     * @return
     */
    protected final Map<String, Object> writeFeaturesToStream(OutputStream outputStream, Iterator<E> features, Argument argument)
            throws IOException {
        PrintWriter writer = new PrintWriter(new OutputStreamWriter(outputStream));

        Map<String, Object> attributes = null;
        if (features != null) {
            FeatureEncoder codec = getEncodingCodec(argument);

            attributes = codec.encodeAll(outputStream, features);

        }
        writer.flush();
        writer.close();

        return attributes;
    }

    protected final String[] genFullCommand(String chr, int start, int end, int zoom) throws IOException {

        List<String> fullCmd = new ArrayList<String>(commands);

        attributes.clear();
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
                case DATA_TRACK:
                    ts = createTempFile((Track) entry.getValue(), arg, chr, start, end, zoom);
                    sVal = new String[]{ts};
                    break;
                case MULTI_FEATURE_TRACK:
                    sVal = createTempFiles((List<FeatureTrack>) entry.getValue(), arg, chr, start, end, zoom);
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
     * Convenience for calling {@link #createTempFile(org.broad.igv.track.Track, org.broad.igv.dev.plugin.Argument, String, int, int, int)};
     * on each track
     *
     * @param tracks
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     * @throws java.io.IOException
     */
    private String[] createTempFiles(List<FeatureTrack> tracks, Argument argument, String chr, int start, int end, int zoom) throws IOException {
        String[] fileNames = new String[tracks.size()];
        int fi = 0;
        for (FeatureTrack track : tracks) {
            fileNames[fi++] = createTempFile(track, argument, chr, start, end, zoom);
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
     * @return String with temp file name.
     * @throws java.io.IOException
     */
    protected abstract String createTempFile(Track track, Argument argument, String chr, int start, int end, int zoom) throws IOException;

    protected final String createTempFile(List<E> features, Argument argument) throws IOException {
        File outFile = File.createTempFile("features", ".tmp", null);
        outFile.deleteOnExit();

        Map<String, Object> attributes = writeFeaturesToStream(new FileOutputStream(outFile), features.iterator(), argument);
        String path = outFile.getAbsolutePath();
        this.attributes.add(attributes);
        return path;
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
    protected final Iterator<D> getFeatures(String chr, int start, int end, int zoom) throws IOException {

        String[] fullCmd = genFullCommand(chr, start, end, zoom);

        //Start plugin process
        Process pr = RuntimeUtils.startExternalProcess(fullCmd, null, null);

        //Read back in the data which plugin output
        FeatureDecoder<D> codec = getDecodingCodec();
        return codec.decodeAll(pr.getInputStream(), strictParsing);
    }

    /**
     * Create ncoding codec, and apply inputs
     *
     * @param argument
     * @return
     */
    protected final FeatureEncoder<E> getEncodingCodec(Argument argument) {
        FeatureEncoder<E> codec = instantiateEncodingCodec(argument);
        codec.setInputs(Collections.unmodifiableList(commands), Collections.unmodifiableMap(arguments));
        return codec;
    }

    /**
     * Get the encoding codec for this argument. Default
     * is IGVBEDCodec, if there was none specified.
     * <p/>
     * Codec will be reflectively instantiated if it was specified
     * in the {@code argument}
     *
     * @param argument
     * @return
     */
    private final FeatureEncoder<E> instantiateEncodingCodec(Argument argument) {
        String encodingCodec = argument.getEncodingCodec();

        if (encodingCodec == null) return new AsciiEncoder(new IGVBEDCodec());

        URL[] libURLs = argument.getLibURLs();

        try {
            ClassLoader loader = URLClassLoader.newInstance(
                    libURLs, getClass().getClassLoader()
            );
            Class clazz = loader.loadClass(encodingCodec);
            Constructor constructor = clazz.getConstructor();
            Object ocodec = constructor.newInstance();
            FeatureEncoder<E> codec;
            if (!(ocodec instanceof FeatureEncoder) && ocodec instanceof LineFeatureEncoder) {
                codec = new AsciiEncoder((LineFeatureEncoder<D>) ocodec);
            } else {
                codec = (FeatureEncoder<E>) ocodec;
            }
            return codec;

        } catch (ClassNotFoundException e) {
            log.error("Could not find class " + encodingCodec, e);
            throw new IllegalArgumentException(e);
        } catch (Exception e) {
            log.error("Exception getting encoding codec", e);
            throw new RuntimeException(e);
        }
    }

    /**
     * Create decoding codec and set inputs and attributes
     *
     * @return
     */
    protected final FeatureDecoder<D> getDecodingCodec() {
        FeatureDecoder<D> codec = instantiateDecodingCodec(decodingCodec, decodingLibURLs);
        codec.setInputs(Collections.unmodifiableList(commands), Collections.unmodifiableMap(arguments));
        codec.setAttributes(Collections.unmodifiableList(attributes));
        return codec;
    }

    /**
     * Instantiate decodingCodec, using {@code libURLs} as classpath. Will throw exceptions
     * if class cannot be instantiated
     *
     * @param decodingCodec
     * @param libURLs
     * @return
     */
    protected final FeatureDecoder<D> instantiateDecodingCodec(String decodingCodec, URL[] libURLs) {
        if (decodingCodec == null) {
            AsciiFeatureCodec<D> asciiCodec = CodecFactory.getCodec("." + format, GenomeManager.getInstance().getCurrentGenome());
            if (asciiCodec == null) {
                throw new IllegalArgumentException("Unable to find codec for format " + format);
            }
            return new AsciiDecoder.DecoderWrapper<D>(asciiCodec);
        }

        try {

            ClassLoader loader = URLClassLoader.newInstance(libURLs,
                    getClass().getClassLoader());

            Class clazz = loader.loadClass(decodingCodec);
            Constructor constructor = clazz.getConstructor();
            Object codec = constructor.newInstance();
            //User can pass in LineFeatureDecoder, we just wrap it
            if (!(codec instanceof FeatureDecoder) && codec instanceof LineFeatureDecoder) {
                return new AsciiDecoder((LineFeatureDecoder<D>) codec);
            }
            return (FeatureDecoder<D>) codec;

        } catch (ClassNotFoundException e) {
            log.error("Could not find class " + decodingCodec, e);
            throw new IllegalArgumentException(e);
        } catch (Exception e) {
            log.error("Exception getting decoding codec", e);
            throw new RuntimeException(e);
        }

    }

}
