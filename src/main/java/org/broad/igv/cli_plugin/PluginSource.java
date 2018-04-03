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

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.IGVBEDCodec;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.AlignmentInterval;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.RuntimeUtils;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.SimpleBEDFeature;
import htsjdk.tribble.readers.PositionalBufferedStream;

import javax.xml.bind.annotation.*;
import javax.xml.bind.annotation.adapters.XmlAdapter;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.io.*;
import java.lang.reflect.Constructor;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.*;

/**
 * A feature source which derives its information
 * from a command line cli_plugin
 * User: jacob
 * Date: 2012/05/01
 */
@XmlAccessorType(XmlAccessType.NONE)
public abstract class PluginSource<E extends Feature, D extends Feature> {

    private static Logger log = Logger.getLogger(PluginSource.class);

    /**
     * Initial command tokens. This will typically include both
     * the executable and command, e.g. {"/usr/bin/bedtools", "intersect"}
     */
    @XmlList
    @XmlAttribute
    protected List<String> commands;

    @XmlJavaTypeAdapter(MyMapAdapter.class)
    protected LinkedHashMap<Argument, Object> arguments;

    /**
     * We store the argument values by id for later retrieval by
     * other arguments, or the parser. For instance, if a parser
     * needs to know what the output file was.
     */

    @XmlElement
    protected PluginSpecReader.Parser parser;

    protected URL[] decodingLibURLs = new URL[0];

    @XmlAttribute
    protected String specPath = null;

    /**
     * For decoding, we may need to know how many columns
     * were written out in the first place
     */
    protected List<Map<String, Object>> attributes = new ArrayList<Map<String, Object>>(2);

    /**
     * Each time we output data, we give it a unique ID.
     * As yet
     */
    protected String lastRunId;
    private static final String RUN_ID_ATTR = "RUN_ID";

    //private QueryTracker queryTracker = new QueryTracker();

    @SubtlyImportant
    protected PluginSource() {
    }

    public PluginSource(List<String> commands, LinkedHashMap<Argument, Object> arguments, PluginSpecReader.Output outputAttrs, String specPath) {
        this.commands = commands;
        this.arguments = arguments;
        this.parser = outputAttrs.parser;
        this.specPath = specPath;

        String[] libs = this.parser.libs;
        libs = libs != null ? libs : new String[]{};
        try {
            decodingLibURLs = PluginSpecReader.getLibURLs(libs, FileUtils.getParent(specPath));
        } catch (MalformedURLException e) {
            log.error("Error parsing library URL", e);
            throw new RuntimeException(e);
        }
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
    protected final Map<String, Object> writeFeaturesToStream(OutputStream outputStream, Iterator features, Argument argument)
            throws IOException {

        Map<String, Object> attributes = null;
        if (features != null) {
            FeatureEncoder codec = getEncodingCodec(argument);
            attributes = codec.encodeAll(outputStream, features);
        }
        outputStream.flush();
        outputStream.close();

        return attributes;
    }

    protected final String[] genFullCommand(String chr, int start, int end, int zoom) throws IOException {

        List<String> fullCmd = new ArrayList<String>(commands);

        attributes.clear();

        String runId = createNewRunId();

        Map<String, String> idVariables = new HashMap<String, String>(arguments.size());
        idVariables.put(RUN_ID_ATTR, runId);

        for (Map.Entry<Argument, Object> entry : arguments.entrySet()) {
            Argument arg = entry.getKey();

            if (!arg.isValidValue(entry.getValue())) {
                String msg = "Type: " + arg.getType() + " value: " + entry.getValue();
                throw new IllegalArgumentException(msg);
            }

            String[] sVal = null;
            String ts = null;
            switch (arg.getType()) {
                case BOOL:
                    boolean selected = (Boolean) entry.getValue();
                    if (selected) {
                        //Output the cmd_arg, but it doesn't take an argument
                    } else {
                        continue;
                    }
                    break;
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
                case ALIGNMENT_TRACK:
                    ts = createTempFile((AlignmentTrack) entry.getValue(), arg, chr, start, end, zoom);
                    sVal = new String[]{ts};
                    break;
                case VARIANT_TRACK:
                case FEATURE_TRACK:
                case DATA_TRACK:
                    ts = createTempFile((Track) entry.getValue(), arg, chr, start, end, zoom);
                    sVal = new String[]{ts};
                    break;
                case MULTI_FEATURE_TRACK:
                    sVal = createTempFiles((List<FeatureTrack>) entry.getValue(), arg, chr, start, end, zoom);
                    break;
                case LOCUS:
                    ts = writeLocus(arg, chr, start, end);
                    sVal = new String[]{ts};
                    break;
            }

            if (arg.getId() != null && sVal != null) {
                idVariables.put(arg.getId(), sVal[0]);
            }

            if (arg.isOutput()) {
                String cmdArg = arg.getCmdArg();
                if (cmdArg.trim().length() > 0) {
                    cmdArg = replaceStringsFromIds(cmdArg, idVariables);
                    fullCmd.add(cmdArg);
                }
                if (sVal != null) {
                    fullCmd.addAll(Arrays.asList(sVal));
                }
            }
        }

        parser.source = replaceStringsFromIds(parser.source, idVariables);

        return fullCmd.toArray(new String[0]);
    }

    protected String writeLocus(Argument arg, String chr, int start, int end) throws IOException{
        Feature feat = new SimpleBEDFeature(start, end, chr);
        return createTempFile(Arrays.asList(feat), arg);
    }

    /**
     * Replace strings of form $"variablename", similar
     * to how the unix shell deals with variables.
     * <p/>
     * This is best demonstrated by example.
     * inputString = "My name is $myname"
     * idVariables = {"myname" -> "bob", "othervariable" -> ted}
     * replaceStringsFromIds(inputString, idVariables) returns
     * "My name is bob"
     *
     * @param inputString String in which to perform replacement
     * @param idVariables May from string -> string, from to.
     * @return
     */
    private String replaceStringsFromIds(String inputString, Map<String, String> idVariables) {
        for (String argId : idVariables.keySet()) {
            inputString = inputString.replace("$" + argId, idVariables.get(argId));
        }
        return inputString;
    }

    protected String createNewRunId() {
        lastRunId = "" + System.currentTimeMillis();
        return lastRunId;
    }

    String getLastRunId() {
        return lastRunId;
    }

    /**
     * Convenience for calling {@link #createTempFile(org.broad.igv.track.Track, org.broad.igv.cli_plugin.Argument, String, int, int, int)};
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

    protected final String createTempFile(List features, Argument argument) throws IOException {
        String ext = ".tmp";
        switch(argument.getType()){
            case ALIGNMENT_TRACK:
                ext += ".sam";
                break;
            case LOCUS:
                ext += ".bed";
                break;
        }
        File outFile = File.createTempFile("features", ext, null);
        outFile.deleteOnExit();

        Map<String, Object> attributes = writeFeaturesToStream(new FileOutputStream(outFile), features.iterator(), argument);
        String path = outFile.getAbsolutePath();
        this.attributes.add(attributes);
        return path;
    }

    protected List<Alignment> getAlignmentsForRange(AlignmentTrack track, String chr, int start, int end, int zoom) throws IOException {

        Collection<AlignmentInterval> loadedIntervals = track.getDataManager().getLoadedIntervals();
        List<Alignment> alignments = new ArrayList<Alignment>();
        for (AlignmentInterval interval : loadedIntervals) {
            if (interval.overlaps(chr, start, end)) {
                Iterator<Alignment> iter = interval.getAlignmentIterator();
                while (iter.hasNext()) {
                    Alignment al = iter.next();
                    if (al.getStart() <= end && al.getEnd() >= start) {
                        alignments.add(al);
                    }
                }
            }
        }

        return alignments;
    }

    /**
     * Perform the actual combination operation between the constituent data
     * sources.
     *
     * @param chr
     * @param start
     * @param end
     * @return
     * @throws java.io.IOException
     */
    protected final Iterator<D> getFeatures(String chr, int start, int end, int zoom) throws IOException {
        if (parser.source == null) {
            throw new IllegalStateException("Null value for source");
        }

        boolean rerun = true;
        //synchronized (this.queryTracker){
        //    rerun = !this.queryTracker.isQuerySame(chr, start, end, zoom);
        //}

        /**
         * A process might generate multiple output files, we only want to run it once.
         * e.g. Cufflinks generates transcripts.gtf, genes.fpkm_tracking, isoforms.fpkm_tracking
         * If the file exists we just read from it
         */

        InputStream dataStream = null;
        if (!rerun) {
            File inFile = new File(parser.source);
            if (inFile.canRead()) {
                //Use the existing file
                dataStream = new FileInputStream(parser.source);
            }
        }

        if (dataStream == null) {
            //synchronized (this.queryTracker){
            String[] fullCmd = genFullCommand(chr, start, end, zoom);
            //log.debug("interval: " + Locus.getFormattedLocusString(chr, start, end));
            //log.debug(StringUtils.join(fullCmd, " "));

            //Start cli_plugin process
            Process pr = RuntimeUtils.startExternalProcess(fullCmd, null, null);

            if (parser.source.equals(PluginSpecReader.Parser.SOURCE_STDOUT)) {
                dataStream = pr.getInputStream();
            } else {
                try {
                    pr.waitFor();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                dataStream = new FileInputStream(parser.source);
            }
            //this.queryTracker.setLastQuery(chr, start, end, zoom);
            //}
        }
        //Read back in the data which cli_plugin output
        FeatureDecoder<D> codec = getDecodingCodec();

        return codec.decodeAll(dataStream, parser.strict);
    }

    /**
     * Create encoding codec, and apply inputs
     *
     * @param argument
     * @return
     */
    protected final FeatureEncoder<E> getEncodingCodec(Argument argument) {
        FeatureEncoder<E> codec = instantiateEncodingCodec(argument);
        codec.setInputs(Collections.unmodifiableList(commands), Collections.unmodifiableMap(arguments), argument);
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
    private final FeatureEncoder instantiateEncodingCodec(Argument argument) {
        String encodingCodec = argument.getEncodingCodec();

        if (encodingCodec == null) {
            switch(argument.getType()){
                case LOCUS:
                case FEATURE_TRACK:
                case MULTI_FEATURE_TRACK:
                    return new AsciiEncoder(new IGVBEDCodec());
                case ALIGNMENT_TRACK:
                    return new SamAlignmentEncoder();
                case VARIANT_TRACK:
                    return new VCFEncoder();
            default:
                throw new IllegalArgumentException("No encoding codec provided and default not available");
            }
        }

        try {
            URL[] libURLs = PluginSpecReader.getLibURLs(argument.getLibPaths(), FileUtils.getParent(specPath));
            if (libURLs == null) libURLs = new URL[0];
            ClassLoader loader = URLClassLoader.newInstance(
                    libURLs, getClass().getClassLoader()
            );
            Class clazz = loader.loadClass(encodingCodec);
            Constructor constructor = clazz.getConstructor();
            Object ocodec = constructor.newInstance();
            FeatureEncoder<E> codec;
            if (!(ocodec instanceof FeatureEncoder) && ocodec instanceof LineFeatureEncoder) {
                codec = new AsciiEncoder((LineFeatureEncoder<E>) ocodec);
            } else {
                codec = (FeatureEncoder<E>) ocodec;
            }
            return codec;

        } catch (ClassNotFoundException e) {
            log.error("Could not find class " + encodingCodec, e);
            throw new IllegalArgumentException(e);
        } catch (MalformedURLException e) {
            log.error("Malformed library URL", e);
            throw new RuntimeException(e);
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
        FeatureDecoder<D> codec = instantiateDecodingCodec(parser.decodingCodec, decodingLibURLs);
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
            FeatureCodec<D, ?> knownCodec = CodecFactory.getCodec("." + parser.format, GenomeManager.getInstance().getCurrentGenome());
            if (knownCodec == null) {
                throw new IllegalArgumentException("Unable to find codec for format " + parser.format);
            } else if (knownCodec instanceof AsciiFeatureCodec) {
                return new AsciiDecoder.DecoderWrapper<D>((AsciiFeatureCodec) knownCodec);
            } else {
                return new FeatureCodecDecoder<D>((FeatureCodec<D, PositionalBufferedStream>) knownCodec);
            }
        }

        try {

            if (libURLs == null) libURLs = new URL[0];

            Class clazz = RuntimeUtils.loadClassForName(decodingCodec, libURLs);
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


    public void getPersistentState() {
        /**
         * The structure here uses similar tags to the plugin XML spec, but is not exactly
         * the same. Because one and only one command was actually used we don't need
         * to nest as much. We also need to store argument values as well as inputs
         */

//        Map<String, String> parentProps = new HashMap<String, String>(5);

//        parentProps.put(DECODING_CODEC, parser.format);
//        //TODO Less ambiguous name to not clash with argument
//        parentProps.put(DECODING_LIBS, StringUtils.join(decodingLibURLs, ","));
//        parentProps.put(FORMAT, parser.format);
//        parentProps.put("specPath", specPath);
//
//        List<RecursiveAttributes> allChildren = new ArrayList<RecursiveAttributes>(4);
//
//
//        List<RecursiveAttributes> persCommands = new ArrayList<RecursiveAttributes>();
//        for(int ii=0; ii < commands.size(); ii++){
//            String command = commands.get(ii);
//            Map<String, String> props = new HashMap<String, String>(1);
//            props.put("value", command);
//            props.put("index", "" + ii);
//            RecursiveAttributes persCommand = new RecursiveAttributes(PluginSpecReader.CMD_ARG, props);
//            persCommands.add(persCommand);
//        }
//        RecursiveAttributes persCommandParent = new RecursiveAttributes(PluginSpecReader.COMMAND, Collections.<String, String>emptyMap(),
//                persCommands);
//
//
//        allChildren.add(persCommandParent);
//
//        for(Map.Entry<Object, Object> entry: arguments.entrySet()){
//            Argument argument = (Argument) entry.getKey();
//            List<String> values = null;
//            String sval;
//            switch (argument.getType()) {
//                case MULTI_FEATURE_TRACK:
//                    List<FeatureTrack> lVal = (List<FeatureTrack>) entry.getValue();
//                    values = new ArrayList<String>(lVal.size());
//                    for(FeatureTrack fTrack: lVal){
//                        values.add(fTrack.getId());
//                    }
//                    break;
//                case LONGTEXT:
//                case TEXT:
//                    sval = (String) entry.getValue();
//                    values = Arrays.asList(sval);
//                    break;
//                case FEATURE_TRACK:
//                case DATA_TRACK:
//                case ALIGNMENT_TRACK:
//                    sval = ((Track) entry.getValue()).getId();
//                    values = Arrays.asList(sval);
//                    break;
//            }
//
//            if(values == null) continue;
//
//            //Generate list of values for the argument
//            List<RecursiveAttributes> persValues = new ArrayList<RecursiveAttributes>();
//            for(String value: values){
//                Map<String, String> cprop = new HashMap<String, String>(1);
//                cprop.put(VALUE, value);
//                persValues.add(new RecursiveAttributes(VALUE, cprop));
//            }
//
//            //Group values together under relevant argument
//            //RecursiveAttributes persArgument = argument.getPersistentState();
//            //persArgument.getChildren().addAll(persValues);
//
//            //Arguments aren't grouped together, just put flat in the children of the highest level
//            //allChildren.add(persArgument);
//        }
//
//        RecursiveAttributes overall = new RecursiveAttributes(getClass().getName(), parentProps, allChildren);
//        return overall;
    }

    public void updateTrackReferences(List<Track> allTracks) {
        MyMapAdapter.updateTrackReferences(arguments, allTracks);
    }

//    public void setQueryTracker(QueryTracker queryTracker){
//        this.queryTracker = queryTracker;
//    }

    static class XmlMap {
        public List<Argument> arg =
                new ArrayList<Argument>();
    }

    public final static class MyMapAdapter extends XmlAdapter<XmlMap, LinkedHashMap<Argument, Object>> {

        @Override
        public LinkedHashMap<Argument, Object> unmarshal(XmlMap v) throws Exception {
            LinkedHashMap<Argument, Object> argumentMap = new LinkedHashMap(v.arg.size());
            for (Argument argument : v.arg) {
                Object oVal = null;
                switch (argument.getType()) {
                    case BOOL:
                        oVal = Boolean.parseBoolean(argument.value.get(0));
                        break;
                    case LONGTEXT:
                    case TEXT:
                        oVal = argument.value.get(0);
                        break;
                    case MULTI_FEATURE_TRACK:
                    case FEATURE_TRACK:
                    case DATA_TRACK:
                    case ALIGNMENT_TRACK:
                        oVal = findTrackReference(argument, null);
                        break;
                }
                argumentMap.put(argument, oVal);
            }
            return argumentMap;
        }

        private static Object findTrackReference(Argument argument, List<Track> allTracks) {
            Object oVal = null;
            switch (argument.getType()) {
                case MULTI_FEATURE_TRACK:
                    List<FeatureTrack> inputTracks = new ArrayList<FeatureTrack>(argument.value.size());

                    for (String trackId : argument.value) {
                        inputTracks.add((FeatureTrack) IGVSessionReader.getMatchingTrack(trackId, allTracks));
                    }
                    oVal = inputTracks;
                    break;
                case FEATURE_TRACK:
                case DATA_TRACK:
                case ALIGNMENT_TRACK:
                    String trackId = argument.value.get(0);
                    oVal = IGVSessionReader.getMatchingTrack(trackId, allTracks);
                    break;
            }
            return oVal;
        }

        public static void updateTrackReferences(Map<Argument, Object> argumentMap, List<Track> allTracks) {
            for (Argument argument : argumentMap.keySet()) {
                //Reference already resolved
                if (argumentMap.get(argument) != null) continue;
                Object oVal = findTrackReference(argument, allTracks);
                argumentMap.put(argument, oVal);
            }
        }

        @Override
        public XmlMap marshal(LinkedHashMap<Argument, Object> v) throws Exception {
            XmlMap outMap = new XmlMap();
            for (Map.Entry<Argument, Object> loopEntry : v.entrySet()) {
                Argument argument = loopEntry.getKey();
                List<String> values = null;
                String sval;
                switch (argument.getType()) {
                    case MULTI_FEATURE_TRACK:
                        List<FeatureTrack> lVal = (List<FeatureTrack>) loopEntry.getValue();
                        values = new ArrayList<String>(lVal.size());
                        for (FeatureTrack fTrack : lVal) {
                            values.add(fTrack.getId());
                        }
                        break;
                    case BOOL:
                        sval = Boolean.toString((Boolean) loopEntry.getValue());
                        values = Arrays.asList(sval);
                        break;
                    case LONGTEXT:
                    case TEXT:
                        sval = (String) loopEntry.getValue();
                        values = Arrays.asList(sval);
                        break;
                    case FEATURE_TRACK:
                    case DATA_TRACK:
                    case ALIGNMENT_TRACK:
                        sval = ((Track) loopEntry.getValue()).getId();
                        values = Arrays.asList(sval);
                        break;
                }

                argument.value = values;
                outMap.arg.add(argument);
            }
            return outMap;
        }

    }


}
