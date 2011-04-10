/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.das.DASFeatureSource;
import org.broad.igv.data.*;
import org.broad.igv.data.expression.GCTDataset;
import org.broad.igv.data.expression.GCTDatasetParser;
import org.broad.igv.data.rnai.RNAIDataSource;
import org.broad.igv.data.rnai.RNAIGCTDatasetParser;
import org.broad.igv.data.rnai.RNAIGeneScoreParser;
import org.broad.igv.data.rnai.RNAIHairpinParser;
import org.broad.igv.data.seg.FreqData;
import org.broad.igv.data.seg.SegmentedAsciiDataSet;
import org.broad.igv.data.seg.SegmentedBinaryDataSet;
import org.broad.igv.data.seg.SegmentedDataSource;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.dranger.DRangerParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.goby.GobyAlignmentQueryReader;
import org.broad.igv.gwas.GWASData;
import org.broad.igv.gwas.GWASParser;
import org.broad.igv.gwas.GWASTrack;
import org.broad.igv.maf.MAFTrack;
import org.broad.igv.maf.conservation.OmegaDataSource;
import org.broad.igv.maf.conservation.OmegaTrack;
import org.broad.igv.sam.reader.IndexNotFoundException;
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.renderer.*;
import org.broad.igv.sam.*;
import org.broad.igv.synteny.BlastMapping;
import org.broad.igv.synteny.BlastParser;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.vcf.VCFTrack;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * User: jrobinso
 * Date: Feb 14, 2010
 */
public class TrackLoader {

    private static Logger log = Logger.getLogger(TrackLoader.class);

    /**
     * Switches on various attributes of locator (mainly locator path extension and whether the locator is indexed)
     * to call the appropriate loading method.
     *
     * @param locator
     * @return
     */
    public List<Track> load(ResourceLocator locator) {

        try {
            String typeString = locator.getType();
            if (typeString == null) {
                typeString = locator.getPath().toLowerCase();
                if (!typeString.endsWith("_sorted.txt") &&
                        (typeString.endsWith(".txt") || typeString.endsWith(
                                ".xls") || typeString.endsWith(".gz"))) {
                    typeString = typeString.substring(0, typeString.lastIndexOf("."));
                }
            }
            typeString = typeString.toLowerCase();

            if (typeString.endsWith(".tbi")) {
                MessageUtils.showMessage("<html><b>Error:</b>File type '.tbi' is not recognized.  If this is a 'tabix' index <br>" +
                        " load the associated gzipped file, which should have an extension of '.gz'");
            }

            //TODO Why is this not inside of isIndexed?
            //dhmay seconding this question -- this appears to be taken care of already in isIndexed()
            // Check for index
            boolean hasIndex = false;
            if (locator.isLocal()) {
                File indexFile = new File(locator.getPath() + ".sai");
                hasIndex = indexFile.exists();
            }

            //This list will hold all new tracks created for this locator
            List<Track> newTracks = new ArrayList<Track>();

            if (typeString.equals("das")) {
                loadDASResource(locator, newTracks);
            } else if (isIndexed(locator.getPath())) {
                loadIndexed(locator, newTracks);
            } else if (typeString.endsWith(".vcf") || typeString.endsWith(".vcf4")) {
                // TODO This is a hack,  vcf files must be indexed.  Fix in next release.
                throw new IndexNotFoundException(locator.getPath());
            } else if (typeString.endsWith("h5") || typeString.endsWith("hbin")) {
                loadH5File(locator, newTracks);
            } else if (typeString.endsWith(".rnai.gct")) {
                loadRnaiGctFile(locator, newTracks);
            } else if (typeString.endsWith(".gct") || typeString.endsWith("res") || typeString.endsWith("tab")) {
                loadGctFile(locator, newTracks);
            } else if (typeString.endsWith(".cn") || typeString.endsWith(".xcn") || typeString.endsWith(".snp") ||
                    typeString.endsWith(".igv") || typeString.endsWith(".loh")) {
                loadIGVFile(locator, newTracks);
            } else if (typeString.endsWith(".mut")) {
                loadMutFile(locator, newTracks);
            } else if (typeString.endsWith(".cbs") || typeString.endsWith(".seg") ||
                    typeString.endsWith("glad") || typeString.endsWith("birdseye_canary_calls")) {
                loadSegFile(locator, newTracks);
            } else if (typeString.endsWith(".seg.zip")) {
                loadBinarySegFile(locator, newTracks);
            } else if (typeString.endsWith(".gistic")) {
                loadGisticFile(locator, newTracks);
            } else if (typeString.endsWith(".gs")) {
                loadRNAiGeneScoreFile(locator, newTracks, RNAIGeneScoreParser.Type.GENE_SCORE);
            } else if (typeString.endsWith(".riger")) {
                loadRNAiGeneScoreFile(locator, newTracks, RNAIGeneScoreParser.Type.POOLED);
            } else if (typeString.endsWith(".hp")) {
                loadRNAiHPScoreFile(locator);
            } else if (typeString.endsWith("gene")) {
                loadGeneFile(locator, newTracks);
            } else if (typeString.contains(".tabblastn") || typeString.endsWith(".orthologs")) {
                loadSyntentyMapping(locator, newTracks);
            } else if (typeString.endsWith(".sam") || typeString.endsWith(".bam") ||
                    typeString.endsWith(".sam.list") || typeString.endsWith(".bam.list") ||
                    typeString.endsWith("_sorted.txt") ||
                    typeString.endsWith(".aligned") || typeString.endsWith(".sai") ||
                    typeString.endsWith(".bai")) {
                loadAlignmentsTrack(locator, newTracks);
            } else if (typeString.endsWith(".bedz") || (typeString.endsWith(".bed") && hasIndex)) {
                loadIndexdBedFile(locator, newTracks);
            } else if (typeString.endsWith(".omega")) {
                loadOmegaTrack(locator, newTracks);
            } else if (typeString.endsWith(".wig") || (typeString.endsWith(".bedgraph")) ||
                    typeString.endsWith("cpg.txt") || typeString.endsWith(".expr")) {
                loadWigFile(locator, newTracks);
            } else if (typeString.endsWith(".list")) {
                loadListFile(locator, newTracks);
            } else if (typeString.contains(".dranger")) {
                loadDRangerFile(locator, newTracks);
            } else if (typeString.endsWith(".ewig.tdf") || (typeString.endsWith(".ewig.ibf"))) {
                loadEwigIBFFile(locator, newTracks);
            } else if (typeString.endsWith(".ibf") || typeString.endsWith(".tdf")) {
                loadTDFFile(locator, newTracks);
            } else if (typeString.endsWith(".psl") || typeString.endsWith(".psl.gz") ||
                    typeString.endsWith(".pslx") || typeString.endsWith(".pslx.gz")) {
                loadPslFile(locator, newTracks);
                //AbstractFeatureParser.getInstanceFor() is called twice.  Wasteful
            } else if (AbstractFeatureParser.getInstanceFor(locator) != null) {
                loadFeatureFile(locator, newTracks);
            } else if (MutationParser.isMutationAnnotationFile(locator)) {
                this.loadMutFile(locator, newTracks);
            } else if (WiggleParser.isWiggle(locator)) {
                loadWigFile(locator, newTracks);
            } else if (locator.getPath().toLowerCase().contains(".maf")) {
                loadMAFTrack(locator, newTracks);
            } else if (GCTDatasetParser.parsableMAGE_TAB(locator)) {
                locator.setDescription("MAGE_TAB");
                loadGctFile(locator, newTracks);
            } else if (IGVDatasetParser.parsableMAGE_TAB(locator)) {
                locator.setDescription("MAGE_TAB");
                loadIGVFile(locator, newTracks);

            }

            // For PLINK GWAS results files
            else if (typeString.endsWith(".logistic") || typeString.endsWith(".linear") || typeString.endsWith(".assoc") || typeString.endsWith(".qassoc") || typeString.endsWith(".gwas")) {
                loadGWASFile(locator, newTracks);

            } else if (GobyAlignmentQueryReader.supportsFileType(locator.getPath())) {
                loadAlignmentsTrack(locator, newTracks);
            } else {
                AttributeManager.getInstance().loadSampleInfo(locator);
            }

            // Track line
            TrackProperties tp = null;
            String trackLine = locator.getTrackLine();
            if (trackLine != null) {
                tp = new TrackProperties();
                ParsingUtils.parseTrackLine(trackLine, tp);
            }

            for (Track track : newTracks) {

                if (locator.getUrl() != null) {
                    track.setUrl(locator.getUrl());
                }
                if (tp != null) {
                    track.setTrackProperties(tp);
                }
                if (locator.getColor() != null) {
                    track.setColor(locator.getColor());
                }
                if (locator.getSampleId() != null) {
                    track.setSampleId(locator.getSampleId());
                }

                IGV.getInstance().getTrackManager().addLoadedType(track.getTrackType());
            }


            return newTracks;
        } catch (DataLoadException dle) {
            throw dle;
        } catch (Exception e) {
            log.error(e);
            throw new DataLoadException(e.getMessage(), locator.getPath());
        }

    }

    private void loadIndexed(ResourceLocator locator, List<Track> newTracks) throws IOException {

        TribbleFeatureSource src = new TribbleFeatureSource(locator.getPath());
        String typeString = locator.getPath();
        //Track t;

        if (typeString.endsWith("vcf") || typeString.endsWith("vcf.gz")) {

            VCFTrack t = new VCFTrack(locator, src);
            // VCF tracks handle their own margin
            t.setMargin(0);
            newTracks.add(t);
        } else {

            // Create feature source and track
            FeatureTrack t = new FeatureTrack(locator, src);
            t.setName(locator.getTrackName());
            //t.setRendererClass(BasicTribbleRenderer.class);

            // Set track properties from header
            Object header = src.getHeader();
            if (header != null && header instanceof FeatureFileHeader) {
                FeatureFileHeader ffh = (FeatureFileHeader) header;
                if (ffh.getTrackType() != null) {
                    t.setTrackType(ffh.getTrackType());
                }
                if (ffh.getTrackProperties() != null) {
                    t.setTrackProperties(ffh.getTrackProperties());
                }

                if (ffh.getTrackType() == TrackType.REPMASK) {
                    t.setHeight(15);
                    t.setPreferredHeight(15);
                }
            }
            newTracks.add(t);
        }

    }

    /**
     * Load the input file as a BED or Attribute (Sample Info) file.  First assume
     * it is a BED file,  if no features are found load as an attribute file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadGeneFile(ResourceLocator locator, List<Track> newTracks) {

        FeatureParser featureParser = AbstractFeatureParser.getInstanceFor(locator);
        if (featureParser != null) {
            List<FeatureTrack> tracks = featureParser.loadTracks(locator);
            newTracks.addAll(tracks);
        }

    }

    private void loadSyntentyMapping(ResourceLocator locator, List<Track> newTracks) {

        List<BlastMapping> mappings = (new BlastParser()).parse(locator.getPath());
        List<org.broad.tribble.Feature> features = new ArrayList<org.broad.tribble.Feature>(mappings.size());
        features.addAll(mappings);

        FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features));
        track.setName(locator.getTrackName());
        // track.setRendererClass(AlignmentBlockRenderer.class);
        newTracks.add(track);
    }

    private void loadDRangerFile(ResourceLocator locator, List<Track> newTracks) {

        DRangerParser parser = new DRangerParser();
        newTracks.addAll(parser.loadTracks(locator));
    }

    /**
     * Load the input file as a feature, muation, or maf (multiple alignment) file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadPslFile(ResourceLocator locator, List<Track> newTracks) throws IOException {

        PSLParser featureParser = new PSLParser();
        List<FeatureTrack> tracks = featureParser.loadTracks(locator);
        newTracks.addAll(tracks);
        for (FeatureTrack t : tracks) {
            t.setMinimumHeight(10);
            t.setHeight(30);
            t.setPreferredHeight(30);
            t.setDisplayMode(Track.DisplayMode.EXPANDED);

        }


    }

    /**
     * Load the input file as a feature, muation, or maf (multiple alignment) file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadFeatureFile(ResourceLocator locator, List<Track> newTracks) throws IOException {

        if (locator.isLocal() && (locator.getPath().endsWith(".bed") ||
                locator.getPath().endsWith(".bed.txt") ||
                locator.getPath().endsWith(".bed.gz"))) {
            //checkSize takes care of warning the user
            if (!checkSize(locator.getPath())) {
                return;
            }
        }

        FeatureParser featureParser = AbstractFeatureParser.getInstanceFor(locator);
        if (featureParser != null) {
            List<FeatureTrack> tracks = featureParser.loadTracks(locator);
            newTracks.addAll(tracks);
        } else if (MutationParser.isMutationAnnotationFile(locator)) {
            this.loadMutFile(locator, newTracks);
        } else if (WiggleParser.isWiggle(locator)) {
            loadWigFile(locator, newTracks);
        } else if (locator.getPath().toLowerCase().contains(".maf")) {
            loadMAFTrack(locator, newTracks);
        }


    }

    /**
     * Load the input file as a feature, muation, or maf (multiple alignment) file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadIndexdBedFile(ResourceLocator locator, List<Track> newTracks) throws IOException {

        File featureFile = new File(locator.getPath());
        File indexFile = new File(locator.getPath() + ".sai");
        FeatureSource src = new IndexedBEDFeatureSource(featureFile, indexFile);
        Track t = new FeatureTrack(locator, src);
        newTracks.add(t);

    }

    /**
     * Load GWAS PLINK result file
     *
     * @param locator
     * @param newTracks
     * @throws IOException
     */


    private void loadGWASFile(ResourceLocator locator, List<Track> newTracks) throws IOException {

        GWASParser gwasParser = new GWASParser(locator);
        GWASData gwasData = gwasParser.parse();

        GWASTrack gwasTrack = new GWASTrack(locator, locator.getPath(), locator.getFileName(), gwasData, gwasParser);
        newTracks.add(gwasTrack);

    }


    private void loadRnaiGctFile(ResourceLocator locator, List<Track> newTracks) {

        RNAIGCTDatasetParser parser = new RNAIGCTDatasetParser(locator);

        Collection<RNAIDataSource> dataSources = parser.parse();
        if (dataSources != null) {
            String path = locator.getPath();
            for (RNAIDataSource ds : dataSources) {
                String trackId = path + "_" + ds.getName();
                DataSourceTrack track = new DataSourceTrack(locator, trackId, ds.getName(), ds);

                // Set attributes.
                track.setAttributeValue("SCREEN", ds.getScreen());
                track.setHeight(80);
                track.setPreferredHeight(80);
                newTracks.add(track);
            }
        }
    }

    private void loadGctFile(ResourceLocator locator, List<Track> newTracks) {

        if (locator.isLocal()) {
            if (!checkSize(locator.getPath())) {
                return;
            }
        }

        GCTDatasetParser parser = null;
        GCTDataset ds = null;

        String fName = locator.getTrackName();

        // TODO -- handle remote resource
        try {
            parser = new GCTDatasetParser(locator, null,
                    GenomeManager.getInstance().getGenomeId());
        } catch (IOException e) {
            log.error("Error creating GCT parser.", e);
            throw new DataLoadException("Error creating GCT parser: " + e, locator.getPath());
        }
        ds = parser.createDataset();
        ds.setName(fName);
        ds.setNormalized(true);
        ds.setLogValues(true);

        /*
         * File outputFile = new File(IGV.DEFAULT_USER_DIRECTORY, file.getName() + ".h5");
         * OverlappingProcessor proc = new OverlappingProcessor(ds);
         * proc.setZoomMax(0);
         * proc.process(outputFile.getAbsolutePath());
         * loadH5File(outputFile, messages, attributeList, group);
         */

        // Counter for generating ID
        TrackProperties trackProperties = ds.getTrackProperties();
        String path = locator.getPath();
        for (String trackName : ds.getTrackNames()) {
            Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
            DatasetDataSource dataSource = new DatasetDataSource(currentGenome, trackName, ds);
            String trackId = path + "_" + trackName;
            Track track = new DataSourceTrack(locator, trackId, trackName, dataSource);
            track.setRendererClass(HeatmapRenderer.class);
            track.setTrackProperties(trackProperties);
            newTracks.add(track);
        }
    }

    private void loadIGVFile(ResourceLocator locator, List<Track> newTracks) {

        if (locator.isLocal()) {
            if (!checkSize(locator.getPath())) {
                return;
            }
        }


        String dsName = locator.getTrackName();
        String currentGenomeId = GenomeManager.getInstance().getGenomeId();

        IGVDataset ds = new IGVDataset(currentGenomeId, locator);
        ds.setName(dsName);

        TrackProperties trackProperties = ds.getTrackProperties();
        String path = locator.getPath();
        TrackType type = ds.getType();
        for (String trackName : ds.getTrackNames()) {

            Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
            DatasetDataSource dataSource = new DatasetDataSource(currentGenome, trackName, ds);
            String trackId = path + "_" + trackName;
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            // track.setRendererClass(HeatmapRenderer.class);
            track.setTrackType(ds.getType());
            track.setTrackProperties(trackProperties);

            if (type == TrackType.ALLELE_FREQUENCY) {
                track.setRendererClass(ScatterplotRenderer.class);
                track.setHeight(40);
                track.setPreferredHeight(40);
            }
            newTracks.add(track);
        }
    }

    private boolean checkSize(String file) {

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SIZE_WARNING)) {
            return true;
        }

        File f = new File(file);
        String tmp = file;
        if (f.exists()) {
            long size = f.length();
            if (file.endsWith(".gz")) {
                size *= 3;
                tmp = file.substring(0, file.length() - 3);
            }

            if (size > 50000000) {
                String message = "";
                if (tmp.endsWith(".bed") || tmp.endsWith(".bed.txt")) {
                    message = "The file " + file + " is large (" + (size / 1000000) + " mb).  It is recommended " +
                            "that large files be indexed using IGVTools or Tabix. Loading un-indexed " +
                            "ascii fies of this size can lead to poor performance or unresponsiveness (freezing).  " +
                            "<br><br>IGVTools can be launched from the <b>Tools</b> menu or separately as a command line program.  " +
                            "See the user guide for more details.<br><br>Click <b>Continue</b> to continue loading, or <b>Cancel</b>" +
                            " to skip this file.";

                } else {

                    message = "The file " + file + " is large (" + (size / 1000000) + " mb).  It is recommended " +
                            "that large files be converted to the binary <i>.tdf</i> format using the IGVTools " +
                            "<b>tile</b> command. Loading  unconverted ascii fies of this size can lead to poor " +
                            "performance or unresponsiveness (freezing).  " +
                            "<br><br>IGVTools can be launched from the <b>Tools</b> menu or separately as a " +
                            "command line program. See the user guide for more details.<br><br>Click <b>Continue</b> " +
                            "to continue loading, or <b>Cancel</b> to skip this file.";
                }

                return ConfirmDialog.optionallyShowConfirmDialog(message, PreferenceManager.SHOW_SIZE_WARNING, true);

            }
        }
        return true;
    }

    private void loadDOTFile(ResourceLocator locator, List<Track> newTracks) {

        //GraphTrack gt = new GraphTrack(locator);
        //gt.setHeight(80);
        //newTracks.add(gt);

    }

    private void loadWigFile(ResourceLocator locator, List<Track> newTracks) {

        if (locator.isLocal()) {
            if (!checkSize(locator.getPath())) {
                return;
            }
        }

        String genome = GenomeManager.getInstance().getGenomeId();

        WiggleDataset ds = (new WiggleParser(locator, genome)).parse();
        TrackProperties props = ds.getTrackProperties();

        // In case of conflict between the resource locator display name and the track properties name,
        // use the resource locator
        String name = props == null ? null : props.getName();
        String label = locator.getName();
        if (name == null) {
            name = locator.getFileName();
        } else if (label != null) {
            props.setName(label);  // erase name rom track properties
        }

        String path = locator.getPath();
        boolean multiTrack = ds.getTrackNames().length > 1;

        for (String heading : ds.getTrackNames()) {

            String trackId = multiTrack ? path + "_" + heading : path;
            String trackName = multiTrack ? heading : name;

            Genome currentGenome = GenomeManager.getInstance().getCurrentGenome();
            DatasetDataSource dataSource = new DatasetDataSource(currentGenome, trackId, ds);

            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            String displayName = (label == null || multiTrack) ? heading : label;
            track.setName(displayName);
            track.setTrackProperties(props);

            track.setTrackType(ds.getType());

            if (ds.getType() == TrackType.EXPR) {
                track.setWindowFunction(WindowFunction.none);
            }


            newTracks.add(track);
        }
    }

    private void loadTDFFile(ResourceLocator locator, List<Track> newTracks) {


        if (log.isDebugEnabled()) {
            log.debug("Loading TDFFile: " + locator.toString());
        }

        TDFReader reader = TDFReader.getReader(locator.getPath());
        TrackType type = reader.getTrackType();

        if (log.isDebugEnabled()) {
            log.debug("Parsing track line ");
        }
        TrackProperties props = null;
        String trackLine = reader.getTrackLine();
        if (trackLine != null && trackLine.length() > 0) {
            props = new TrackProperties();
            ParsingUtils.parseTrackLine(trackLine, props);
        }

        // In case of conflict between the resource locator display name and the track properties name,
        // use the resource locator
        String name = locator.getName();
        if (name != null && props != null) {
            props.setName(name);
        }

        if (name == null) {
            name = props == null ? locator.getTrackName() : props.getName();
        }

        int trackNumber = 0;
        String path = locator.getPath();
        boolean multiTrack = reader.getTrackNames().length > 1;

        for (String heading : reader.getTrackNames()) {

            String trackId = multiTrack ? path + "_" + heading : path;
            String trackName = multiTrack ? heading : name;
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, new TDFDataSource(reader, trackNumber, heading));

            String displayName = (name == null || multiTrack) ? heading : name;
            track.setName(displayName);
            track.setTrackType(type);
            if (props != null) {
                track.setTrackProperties(props);
            }
            newTracks.add(track);
            trackNumber++;
        }
    }

    private void loadEwigIBFFile(ResourceLocator locator, List<Track> newTracks) {

        TDFReader reader = TDFReader.getReader(locator.getPath());
        TrackProperties props = null;
        String trackLine = reader.getTrackLine();
        if (trackLine != null && trackLine.length() > 0) {
            props = new TrackProperties();
            ParsingUtils.parseTrackLine(trackLine, props);
        }

        EWigTrack track = new EWigTrack(locator);
        if (props != null) {
            track.setTrackProperties(props);
        }
        track.setName(locator.getTrackName());
        newTracks.add(track);
    }

    private void loadListFile(ResourceLocator locator, List<Track> newTracks) {
        try {
            FeatureSource source = new FeatureDirSource(locator);
            FeatureTrack track = new FeatureTrack(locator, source);
            track.setName(locator.getTrackName());
            track.setVisibilityWindow(100000);
            newTracks.add(track);
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }

    private void loadGisticFile(ResourceLocator locator, List<Track> newTracks) {

        GisticTrack track = GisticFileParser.loadData(locator);
        track.setName(locator.getTrackName());
        newTracks.add(track);
    }

    /**
     * Load a rnai gene score file and create a datasource and track.
     * <p/>
     * // TODO -- change parser to use resource locator rather than path.
     *
     * @param locator
     * @param newTracks
     */
    private void loadRNAiGeneScoreFile(ResourceLocator locator,
                                       List<Track> newTracks, RNAIGeneScoreParser.Type type) {

        String genomeId = GenomeManager.getInstance().getGenomeId();
        RNAIGeneScoreParser parser = new RNAIGeneScoreParser(locator.getPath(), genomeId, type);

        Collection<RNAIDataSource> dataSources = parser.parse();
        String path = locator.getPath();
        for (RNAIDataSource ds : dataSources) {
            String name = ds.getName();
            String trackId = path + "_" + name;
            DataSourceTrack track = new DataSourceTrack(locator, trackId, name, ds);

            // Set attributes.  This "hack" is neccessary to register these attributes with the
            // attribute manager to get displayed.
            track.setAttributeValue("SCREEN", ds.getScreen());
            if ((ds.getCondition() != null) && (ds.getCondition().length() > 0)) {
                track.setAttributeValue("CONDITION", ds.getCondition());
            }
            track.setHeight(80);
            track.setPreferredHeight(80);
            //track.setDataRange(new DataRange(-3, 0, 3));
            newTracks.add(track);
        }

    }

    /**
     * Load a RNAi haripin score file.  The results of this action are hairpin scores
     * added to the RNAIDataManager.  Currently no tracks are created for hairpin
     * scores, although this could change.
     * <p/>
     * // TODO -- change parser to use resource locator rather than path.
     *
     * @param locator
     */
    private void loadRNAiHPScoreFile(ResourceLocator locator) {
        (new RNAIHairpinParser(locator.getPath())).parse();
    }

    private void loadMAFTrack(ResourceLocator locator, List<Track> newTracks) {
        MAFTrack t = new MAFTrack(locator);
        t.setName("Multiple Alignments");
        newTracks.add(t);
    }

    private void loadOmegaTrack(ResourceLocator locator, List<Track> newTracks) {
        OmegaDataSource ds = new OmegaDataSource();
        OmegaTrack track = new OmegaTrack(locator, ds);
        track.setName("Conservation (Omega)");
        track.setHeight(40);
        track.setPreferredHeight(40);
        newTracks.add(track);
    }

    /**
     * Load a rnai gene score file and create a datasource and track.
     * <p/>
     * // TODO -- change parser to use resource locator rather than path.
     *
     * @param locator
     * @param newTracks
     */
    private void loadAlignmentsTrack(ResourceLocator locator, List<Track> newTracks) {

        try {
            String dsName = locator.getTrackName();
            String fn = locator.getPath().toLowerCase();
            boolean isBed = fn.endsWith(".bedz") || fn.endsWith(".bed") || fn.endsWith(".bed.gz");

            // If the user tried to load the index,  look for the file (this is a common mistake)
            if (locator.getPath().endsWith(".sai") || locator.getPath().endsWith(".bai")) {
                MessageUtils.showMessage("<html><b>ERROR:</b> Loading SAM/BAM index files are not supported:  " + locator.getPath() +
                        "<br>Load the SAM or BAM file directly. ");
                return;
            }
            AlignmentDataManager dataManager = new AlignmentDataManager(locator);

            if (locator.getPath().toLowerCase().endsWith(".bam")) {
                if (!dataManager.hasIndex()) {
                    MessageUtils.showMessage("<html>Could not load index file for: " +
                            locator.getPath() + "<br>  An index file is required for SAM & BAM files.");
                    return;
                }
            }

            AlignmentTrack alignmentTrack = new AlignmentTrack(locator, dataManager);    // parser.loadTrack(locator, dsName);
            alignmentTrack.setName(dsName);
            if (isBed) {
                alignmentTrack.setRenderer(new BedRenderer());
                alignmentTrack.setPreferredHeight(40);
                alignmentTrack.setHeight(40);
            }

            // Create coverage track
            CoverageTrack covTrack = new CoverageTrack(locator.getPath() + "_coverage", alignmentTrack.getName() + " Coverage");
            covTrack.setVisible(PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_COV_TRACK));
            newTracks.add(covTrack);
            alignmentTrack.setCoverageTrack(covTrack);
            if (!isBed) {
                covTrack.setDataManager(dataManager);
                dataManager.setCoverageTrack(covTrack);
            }

            // Search for precalculated coverage data
            String covPath = locator.getCoverage();
            if (covPath == null) {
                String path = locator.getPath();
                covPath = path + ".tdf";
            }
            if (covPath != null) {
                try {
                    if ((new File(covPath)).exists() || (IGVHttpUtils.isURL(covPath) &&
                            IGVHttpUtils.resourceAvailable(new URL(covPath)))) {
                        TDFReader reader = TDFReader.getReader(covPath);
                        TDFDataSource ds = new TDFDataSource(reader, 0, alignmentTrack.getName() + " coverage");
                        covTrack.setDataSource(ds);
                    }
                } catch (MalformedURLException e) {
                    // This is expected if
                    //    log.info("Could not loading coverage data: MalformedURL: " + covPath);
                }
            }

            boolean showSpliceJunctionTrack = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK);
            if (showSpliceJunctionTrack) {
                SpliceJunctionFinderTrack spliceJunctionTrack = new SpliceJunctionFinderTrack(locator.getPath() + "_junctions",
                        alignmentTrack.getName() + " Junctions", dataManager);
//            spliceJunctionTrack.setDataManager(dataManager);
                spliceJunctionTrack.setHeight(60);
                spliceJunctionTrack.setPreferredHeight(60);
                spliceJunctionTrack.setVisible(showSpliceJunctionTrack);
                newTracks.add(spliceJunctionTrack);
                alignmentTrack.setSpliceJunctionTrack(spliceJunctionTrack);
            }

            newTracks.add(alignmentTrack);

        } catch (IndexNotFoundException e) {
            MessageUtils.showMessage("<html>Could not find the index file for  <br><br>&nbsp;&nbsp;" + e.getSamFile() +
                    "<br><br>Note: The index file can be created using igvtools and must be in the same directory as the .sam file.");
        } catch (Exception e) {
            MessageUtils.showMessage("Error loading " + locator.getPath() + ": " + e.getMessage());
        }
    }


    /**
     * Load a ".mut" file (muation file) and create tracks.
     *
     * @param locator
     * @param newTracks
     */
    private void loadMutFile(ResourceLocator locator, List<Track> newTracks) {

        MutationParser parser = new MutationParser();
        List<FeatureTrack> mutationTracks = parser.loadMutationTracks(locator);
        for (FeatureTrack track : mutationTracks) {
            track.setTrackType(TrackType.MUTATION);
            track.setRendererClass(MutationRenderer.class);
            newTracks.add(track);
        }
    }

    private void loadSegFile(ResourceLocator locator, List<Track> newTracks) {

        // TODO - -handle remote resource
        SegmentedAsciiDataSet ds = new SegmentedAsciiDataSet(locator);
        String path = locator.getPath();
        TrackProperties props = ds.getTrackProperties();

        // The "freq" track.  TODO - make this optional
        if (ds.getSampleNames().size() > 1) {
            FreqData fd = new FreqData(ds);
            String freqTrackId = path;
            String freqTrackName = (new File(path)).getName();
            CNFreqTrack freqTrack = new CNFreqTrack(locator, freqTrackId, freqTrackName, fd);
            newTracks.add(freqTrack);
        }

        for (String trackName : ds.getDataHeadings()) {
            String trackId = path + "_" + trackName;
            SegmentedDataSource dataSource = new SegmentedDataSource(trackName, ds);
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);
            track.setRendererClass(HeatmapRenderer.class);
            track.setTrackType(ds.getType());

            if (props != null) {
                track.setTrackProperties(props);
            }

            newTracks.add(track);
        }
    }

    private void loadBinarySegFile(ResourceLocator locator, List<Track> newTracks) {

        SegmentedBinaryDataSet ds = new SegmentedBinaryDataSet(locator);
        String path = locator.getPath();

        // The "freq" track.  Make this optional?
        FreqData fd = new FreqData(ds);
        String freqTrackId = path;
        String freqTrackName = (new File(path)).getName();
        CNFreqTrack freqTrack = new CNFreqTrack(locator, freqTrackId, freqTrackName, fd);
        newTracks.add(freqTrack);


        for (String trackName : ds.getSampleNames()) {
            String trackId = path + "_" + trackName;
            SegmentedDataSource dataSource = new SegmentedDataSource(trackName, ds);
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);
            track.setRendererClass(HeatmapRenderer.class);
            track.setTrackType(ds.getType());
            newTracks.add(track);
        }
    }


    /**
     * Load the data from an HDF5 file and add resulting tracks to the supplied TrackGroup.
     * Error messages are appended to the MessageCollection
     *
     * @param locator
     */
    private void loadH5File(ResourceLocator locator, List<Track> newTracks) {

        TrackSet trackSet = null;

        // TODO -- temporary until "name" property is added to locator
        // TODO -- "label" has been added, how does that affect this?

        String fName = locator.getTrackName();

        HDFDataManager dataManager = HDFDataManagerFactory.getDataManager(locator);
        TrackProperties trackProperties = dataManager.getTrackProperties();
        String[] trackNames = dataManager.getTrackNames();

        List<Track> tracks = new ArrayList();

        for (int trackNumber = 0; trackNumber < trackNames.length; trackNumber++) {
            String name = trackNames.length == 1 ? fName : trackNames[trackNumber];

            Track track = null;
            try {
                track = new HDFDataTrack(dataManager, locator, name, trackNumber);
            } catch (FileNotFoundException fe) {
                throw new RuntimeException(fe);
            }
            if (trackProperties != null) {
                track.setTrackProperties(trackProperties);
            }

            tracks.add(track);
        }

        trackSet = new TrackSet(tracks);

        if (trackSet.isEmpty()) {
            throw new RuntimeException("No data found in file");
        }

        newTracks.addAll(trackSet.getTracks());
    }

    private void loadDASResource(ResourceLocator locator, List<Track> currentTracks) {

        //TODO Connect and get all the attributes of the DAS server, and run the appropriate load statements
        //TODO Currently we are only going to be doing features
        // TODO -- move the source creation to a factory


        DASFeatureSource featureSource = null;
        try {
            featureSource = new DASFeatureSource(locator);
        } catch (MalformedURLException e) {
            log.error("Malformed URL", e);
            throw new DataLoadException("Error: Malformed URL ", locator.getPath());
        }

        FeatureTrack track = new FeatureTrack(locator, featureSource);

        // Try to create a sensible name from the path
        String name = locator.getPath();
        if (locator.getPath().contains("genome.ucsc.edu")) {
            name = featureSource.getType();
        } else {
            name = featureSource.getPath().replace("/das/", "").replace("/features", "");
        }
        track.setName(name);

        // A hack until we can notate this some other way
        if (locator.getPath().contains("cosmic")) {
            track.setRendererClass(CosmicFeatureRenderer.class);
            track.setMinimumHeight(2);
            track.setHeight(20);
            track.setPreferredHeight(20);
            track.setDisplayMode(Track.DisplayMode.EXPANDED);
        } else {
            track.setRendererClass(IGVFeatureRenderer.class);
            track.setMinimumHeight(35);
            track.setHeight(45);
            track.setPreferredHeight(45);
        }
        currentTracks.add(track);
    }


    public static boolean isIndexed(String path) {
        String indexExtension = path.endsWith("gz") ? ".tbi" : ".idx";
        String indexPath = path + indexExtension;
        try {
            if (IGVHttpUtils.isURL(path)) {
                return IGVHttpUtils.resourceAvailable(new URL(indexPath));
            } else {
                File f = new File(path + indexExtension);
                return f.exists();
            }

        } catch (IOException e) {
            return false;
        }

    }
}
