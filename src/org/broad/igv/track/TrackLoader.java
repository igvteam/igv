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
import org.broad.igv.PreferenceManager;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bigwig.BigWigDataSource;
import org.broad.igv.das.DASFeatureSource;
import org.broad.igv.data.*;
import org.broad.igv.data.expression.ExpressionDataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.data.rnai.RNAIDataSource;
import org.broad.igv.data.rnai.RNAIGCTDatasetParser;
import org.broad.igv.data.rnai.RNAIGeneScoreParser;
import org.broad.igv.data.rnai.RNAIHairpinParser;
import org.broad.igv.data.seg.*;
import org.broad.igv.dev.SegmentedReader;
import org.broad.igv.dev.affective.AffectiveAnnotationParser;
import org.broad.igv.dev.affective.AffectiveAnnotationTrack;
import org.broad.igv.dev.affective.AffectiveUtils;
import org.broad.igv.dev.affective.Annotation;
import org.broad.igv.dev.db.DBTable;
import org.broad.igv.dev.db.SQLCodecSource;
import org.broad.igv.dev.db.SampleInfoSQLReader;
import org.broad.igv.dev.db.SegmentedSQLReader;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.*;
import org.broad.igv.feature.dranger.DRangerParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.goby.GobyAlignmentQueryReader;
import org.broad.igv.goby.GobyCountArchiveDataSource;
import org.broad.igv.gwas.GWASData;
import org.broad.igv.gwas.GWASParser;
import org.broad.igv.gwas.GWASTrack;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.lists.VariantListManager;
import org.broad.igv.maf.MultipleAlignmentTrack;
import org.broad.igv.methyl.MethylTrack;
import org.broad.igv.peaks.PeakTrack;
import org.broad.igv.renderer.*;
import org.broad.igv.sam.*;
import org.broad.igv.sam.reader.IndexNotFoundException;
import org.broad.igv.synteny.BlastMapping;
import org.broad.igv.synteny.BlastParser;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.util.PedigreeUtils;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * User: jrobinso
 * Date: Feb 14, 2010
 */
public class TrackLoader {

    private static Logger log = Logger.getLogger(TrackLoader.class);

    /**
     * Calls {@linkplain TrackLoader#load(org.broad.igv.util.ResourceLocator, org.broad.igv.feature.genome.Genome)}
     * with genome from IGV instance (if not null).
     *
     * @param locator
     * @param igv
     * @return
     */
    public List<Track> load(ResourceLocator locator, IGV igv) {
        Genome genome = igv != null ? GenomeManager.getInstance().getCurrentGenome() : null;
        return load(locator, genome);
    }

    /**
     * Switches on various attributes of locator (mainly locator path extension and whether the locator is indexed)
     * to call the appropriate loading method.
     *
     * @param locator
     * @param genome
     * @return
     */
    public List<Track> load(ResourceLocator locator, Genome genome) {

        final String path = locator.getPath();
        try {
            String typeString = locator.getType();
            if (typeString == null) {

                // Genome space hack -- check for explicit type converter
                //https://dmtest.genomespace.org:8444/datamanager/files/users/SAGDemo/Step1/TF.data.tab
                //   ?dataformat=http://www.genomespace.org/datamanager/dataformat/gct/0.0.0
                if (path.contains("?dataformat")) {
                    if (path.contains("dataformat/gct")) {
                        typeString = ".gct";
                    } else if (path.contains("dataformat/bed")) {
                        typeString = ".bed";
                    } else if (path.contains("dataformat/cn")) {
                        typeString = ".cn";
                    }

                } else {
                    typeString = path.toLowerCase();
                    if (!typeString.endsWith("_sorted.txt") &&
                            (typeString.endsWith(".txt") || typeString.endsWith(
                                    ".xls") || typeString.endsWith(".gz"))) {
                        typeString = typeString.substring(0, typeString.lastIndexOf("."));
                    }
                }
            }
            typeString = typeString.toLowerCase();

            if (typeString.endsWith(".tbi")) {
                MessageUtils.showMessage("<html><b>Error:</b>File type '.tbi' is not recognized.  If this is a 'tabix' index <br>" +
                        " load the associated gzipped file, which should have an extension of '.gz'");
            }

            //This list will hold all new tracks created for this locator
            List<Track> newTracks = new ArrayList<Track>();

            String serverURL = locator.getServerURL();
            if (serverURL != null && serverURL.startsWith("jdbc:")) {
                this.loadFromDatabase(locator, newTracks, genome);
            } else if (typeString.endsWith(".dbxml")) {
                loadFromDBProfile(locator, newTracks);
            } else if (typeString.endsWith(".gmt")) {
                loadGMT(locator);
            } else if (typeString.equals("das")) {
                loadDASResource(locator, newTracks);
            } else if (isIndexed(path, genome)) {
                loadIndexed(locator, newTracks, genome);
            } else if (typeString.endsWith(".vcf.list")) {
                loadVCFListFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".vcf") || typeString.endsWith(".vcf4")) {
                // VCF files must be indexed.
                throw new IndexNotFoundException(path);
            } else if (typeString.endsWith(".trio")) {
                loadTrioData(locator);
            } else if (typeString.endsWith("varlist")) {
                VariantListManager.loadVariants(locator);
            } else if (typeString.endsWith("samplepathmap")) {
                VariantListManager.loadSamplePathMap(locator);
            } else if (typeString.endsWith("h5") || typeString.endsWith("hbin")) {
                throw new DataLoadException("HDF5 files are no longer supported", locator.getPath());
            } else if (typeString.endsWith(".rnai.gct")) {
                loadRnaiGctFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".gct") || typeString.endsWith("res") || typeString.endsWith("tab")) {
                loadGctFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".cn") || typeString.endsWith(".xcn") || typeString.endsWith(".snp") ||
                    typeString.endsWith(".igv") || typeString.endsWith(".loh")) {
                loadIGVFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".cbs") || typeString.endsWith(".seg") ||
                    typeString.endsWith("glad") || typeString.endsWith("birdseye_canary_calls")
                    || typeString.endsWith(".seg.zip")) {
                loadSegFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".gistic")) {
                loadGisticFile(locator, newTracks);
            } else if (typeString.endsWith(".gs")) {
                loadRNAiGeneScoreFile(locator, newTracks, RNAIGeneScoreParser.Type.GENE_SCORE, genome);
            } else if (typeString.endsWith(".riger")) {
                loadRNAiGeneScoreFile(locator, newTracks, RNAIGeneScoreParser.Type.POOLED, genome);
            } else if (typeString.endsWith(".hp")) {
                loadRNAiHPScoreFile(locator);
            } else if (typeString.endsWith("gene")) {
                loadGeneFile(locator, newTracks, genome);
            } else if (typeString.contains(".tabblastn") || typeString.endsWith(".orthologs")) {
                loadSyntentyMapping(locator, newTracks);
            } else if (typeString.endsWith(".sam") || typeString.endsWith(".bam") ||
                    typeString.endsWith(".sam.list") || typeString.endsWith(".bam.list") ||
                    typeString.endsWith("_sorted.txt") ||
                    typeString.endsWith(".aligned") || typeString.endsWith(".sai") ||
                    typeString.endsWith(".bai") || typeString.equals("alist")) {
                loadAlignmentsTrack(locator, newTracks, genome);
            } else if (typeString.endsWith(".wig") || (typeString.endsWith(".bedgraph")) ||
                    typeString.endsWith("cpg.txt") || typeString.endsWith(".expr")) {
                loadWigFile(locator, newTracks, genome);
            } else if (typeString.contains(".dranger")) {
                loadDRangerFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".ewig.tdf") || (typeString.endsWith(".ewig.ibf"))) {
                loadEwigIBFFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".bw") || typeString.endsWith(".bb") || typeString.endsWith(".bigwig") ||
                    typeString.endsWith(".bigbed")) {
                loadBWFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".ibf") || typeString.endsWith(".tdf")) {
                loadTDFFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".counts")) {
                loadGobyCountsArchive(locator, newTracks, genome);
            } else if (GFFFeatureSource.isGFF(locator.getPath())) {
                loadGFFfile(locator, newTracks, genome);
            } else if (AbstractFeatureParser.canParse(locator.getPath())) {
                loadFeatureFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".mut")) { //  MutationParser.isMutationAnnotationFile(locator)) {
                this.loadMutFile(locator, newTracks, genome);
            } else if (WiggleParser.isWiggle(locator)) {
                loadWigFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".maf") || typeString.endsWith(".maf.annotated")) {
                if (MutationTrackLoader.isMutationAnnotationFile(locator)) {
                    loadMutFile(locator, newTracks, genome);
                } else {
                    loadMultipleAlignmentTrack(locator, newTracks, genome);
                }
            } else if (typeString.endsWith(".maf.dict")) {
                loadMultipleAlignmentTrack(locator, newTracks, genome);
            } else if (path.toLowerCase().contains(".peak.bin")) {
                loadPeakTrack(locator, newTracks, genome);
            } else if ("mage-tab".equals(locator.getType()) || ExpressionFileParser.parsableMAGE_TAB(locator)) {
                locator.setDescription("MAGE_TAB");
                loadGctFile(locator, newTracks, genome);
            } else if (GWASParser.isGWASFile(typeString)) {
                loadGWASFile(locator, newTracks, genome);
            } else if (GobyAlignmentQueryReader.supportsFileType(path)) {
                loadAlignmentsTrack(locator, newTracks, genome);
            } else if (typeString.endsWith(".list")) {
                // This should be deprecated
                loadListFile(locator, newTracks, genome);
            } else if (path.contains("Participant") && path.endsWith(".csv")) {
                loadAffectiveAnnotationTrack(locator, newTracks, genome);
            } else if (AttributeManager.isSampleInfoFile(locator)) {
                // This might be a sample information file.
                AttributeManager.getInstance().loadSampleInfo(locator);
            } else {
                MessageUtils.showMessage("<html>Unknown file type: " + path + "<br>Check file extenstion");
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
                    track.setProperties(tp);
                }
                if (locator.getColor() != null) {
                    track.setColor(locator.getColor());
                }
                if (locator.getSampleId() != null) {
                    track.setSampleId(locator.getSampleId());
                }
            }


            return newTracks;
        } catch (IOException e) {
            log.error(e.getMessage(), e);
            throw new DataLoadException(e.getMessage(), path);
        }

    }

    private void loadGMT(ResourceLocator locator) throws IOException {
        List<GeneList> lists = GeneListManager.getInstance().importGMTFile(locator.getPath());
        if (lists.size() == 1) {
            GeneList gl = lists.get(0);
            IGV.getInstance().setGeneList(gl.getName(), true);
        } else {
            MessageUtils.showMessage("Loaded " + lists.size() + " gene lists.");
        }
    }

    private void loadIndexed(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        TribbleFeatureSource src = GFFFeatureSource.isGFF(locator.getPath()) ?
                new GFFFeatureSource(locator.getPath(), genome) :
                new TribbleFeatureSource(locator.getPath(), genome);
        String typeString = locator.getPath();
        //Track t;

        if (typeString.endsWith("vcf") || typeString.endsWith("vcf.gz")) {

            VCFHeader header = (VCFHeader) src.getHeader();

            // Test if the input VCF file contains methylation rate data:

            // This is determined by testing for the presence of two sample format fields: MR and GB, used in the
            // rendering of methylation rate.
            // MR is the methylation rate on a scale of 0 to 100% and GB is the number of bases that pass
            // filter for the position. GB is needed to avoid displaying positions for which limited coverage
            // prevents reliable estimation of methylation rate.
            boolean enableMethylationRateSupport = (header.getFormatHeaderLine("MR") != null &&
                    header.getFormatHeaderLine("GB") != null);

            List<String> allSamples = new ArrayList(header.getGenotypeSamples());

            VariantTrack t = new VariantTrack(locator, src, allSamples, enableMethylationRateSupport);

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
                    t.setProperties(ffh.getTrackProperties());
                }

                if (ffh.getTrackType() == TrackType.REPMASK) {
                    t.setHeight(15);
                }
            }
            newTracks.add(t);
        }

    }


    private void loadVCFListFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        TribbleListFeatureSource src = new TribbleListFeatureSource(locator.getPath(), genome);

        VCFHeader header = (VCFHeader) src.getHeader();

        // Test if the input VCF file contains methylation rate data:

        // This is determined by testing for the presence of two sample format fields: MR and GB, used in the
        // rendering of methylation rate.
        // MR is the methylation rate on a scale of 0 to 100% and GB is the number of bases that pass
        // filter for the position. GB is needed to avoid displaying positions for which limited coverage
        // prevents reliable estimation of methylation rate.
        boolean enableMethylationRateSupport = (header.getFormatHeaderLine("MR") != null &&
                header.getFormatHeaderLine("GB") != null);

        List<String> allSamples = new ArrayList(header.getGenotypeSamples());

        VariantTrack t = new VariantTrack(locator, src, allSamples, enableMethylationRateSupport);

        // VCF tracks handle their own margin
        t.setMargin(0);
        newTracks.add(t);

    }

    private void loadGeneFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        FeatureParser featureParser = AbstractFeatureParser.getInstanceFor(locator, genome);
        if (featureParser != null) {
            List<FeatureTrack> tracks = featureParser.loadTracks(locator, genome);
            newTracks.addAll(tracks);
        }

    }

    private void loadSyntentyMapping(ResourceLocator locator, List<Track> newTracks) {

        List<BlastMapping> mappings = (new BlastParser()).parse(locator.getPath());
        List<org.broad.tribble.Feature> features = new ArrayList<org.broad.tribble.Feature>(mappings.size());
        features.addAll(mappings);

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features, genome));
        track.setName(locator.getTrackName());
        // track.setRendererClass(AlignmentBlockRenderer.class);
        newTracks.add(track);
    }

    private void loadDRangerFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        DRangerParser parser = new DRangerParser();
        newTracks.addAll(parser.loadTracks(locator, genome));
    }


    /**
     * Load the input file as a feature, mutation, or maf (multiple alignment) file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadGFFfile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        GFFParser featureParser = new GFFParser(locator.getPath());
        List<FeatureTrack> tracks = featureParser.loadTracks(locator, genome);
        newTracks.addAll(tracks);
    }

    /**
     * Load the input file as a feature, mutation, or maf (multiple alignment) file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadFeatureFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        if (locator.isLocal() && (locator.getPath().endsWith(".bed") ||
                locator.getPath().endsWith(".bed.txt"))) {
            //checkSize takes care of warning the user
            if (!checkSize(locator.getPath())) {
                return;
            }
        }

        FeatureCodec codec = CodecFactory.getCodec(locator.getPath(), genome);
        if (codec != null) {
            AbstractFeatureReader<Feature> bfs = AbstractFeatureReader.getFeatureReader(locator.getPath(), codec, false);
            Iterable<Feature> iter = bfs.iterator();
            Object header = bfs.getHeader();
            TrackProperties trackProperties = getTrackProperties(header);
            List<FeatureTrack> tracks = AbstractFeatureParser.loadTracks(iter, locator, genome, trackProperties);

            if (locator.getPath().contains(".narrowPeak") || locator.getPath().contains(".broadPeak")) {
                for (FeatureTrack t : tracks) {
                    t.setUseScore(true);
                }
            }


            newTracks.addAll(tracks);
        } else if (MutationTrackLoader.isMutationAnnotationFile(locator)) {
            this.loadMutFile(locator, newTracks, genome);
        } else if (WiggleParser.isWiggle(locator)) {
            loadWigFile(locator, newTracks, genome);
        } else if (locator.getPath().toLowerCase().contains(".maf") || locator.getPath().toLowerCase().endsWith(".maf.dict")) {
            loadMultipleAlignmentTrack(locator, newTracks, genome);
        }
    }

    /**
     * Load GWAS PLINK result file
     *
     * @param locator
     * @param newTracks
     * @throws IOException
     */


    private void loadGWASFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        GWASParser gwasParser = new GWASParser(locator, genome);
        GWASData gwasData = gwasParser.parse();

        GWASTrack gwasTrack = new GWASTrack(locator, locator.getPath(), locator.getFileName(), gwasData, gwasParser);
        newTracks.add(gwasTrack);

    }


    private void loadRnaiGctFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        RNAIGCTDatasetParser parser = new RNAIGCTDatasetParser(locator, genome);

        Collection<RNAIDataSource> dataSources = parser.parse();
        if (dataSources != null) {
            String path = locator.getPath();
            for (RNAIDataSource ds : dataSources) {
                String trackId = path + "_" + ds.getName();
                DataSourceTrack track = new DataSourceTrack(locator, trackId, ds.getName(), ds);

                // Set attributes.
                track.setAttributeValue("SCREEN", ds.getScreen());
                track.setHeight(80);
                newTracks.add(track);
            }
        }
    }

    private void loadGctFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        if (locator.isLocal()) {
            if (!checkSize(locator.getPath())) {
                return;
            }
        }

        ExpressionFileParser parser = null;
        ExpressionDataset ds = null;
        parser = new ExpressionFileParser(locator, null, genome);
        ds = parser.createDataset();
        if (ds.isEmpty()) {
            String message = "The probes in the file <br>&nbsp;&nbsp;&nbsp;" + locator.getPath() + "<br>" +
                    "could not be mapped to genomic positions.  This can be corrected by specify a probe mapping<br>" +
                    "file from the Preferences window (Probes tab), or by specifing the genomic positions in the<br>" +
                    "expression data file.  Please see the user guide for more details.";
            MessageUtils.showMessage(message);

        } else {
            ds.setName(locator.getTrackName());
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
                DatasetDataSource dataSource = new DatasetDataSource(trackName, ds, genome);
                String trackId = path + "_" + trackName;
                Track track = new DataSourceTrack(locator, trackId, trackName, dataSource);
                track.setRendererClass(HeatmapRenderer.class);
                track.setProperties(trackProperties);
                newTracks.add(track);
            }
        }


    }

    private void loadIGVFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        if (locator.isLocal()) {
            if (!checkSize(locator.getPath())) {
                return;
            }
        }


        String dsName = locator.getTrackName();
        IGVDataset ds = new IGVDataset(locator, genome);
        ds.setName(dsName);

        TrackProperties trackProperties = ds.getTrackProperties();
        String path = locator.getPath();
        TrackType type = ds.getType();
        for (String trackName : ds.getTrackNames()) {

            DatasetDataSource dataSource = new DatasetDataSource(trackName, ds, genome);
            String trackId = path + "_" + trackName;
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            // track.setRendererClass(HeatmapRenderer.class);
            track.setTrackType(ds.getType());
            track.setProperties(trackProperties);

            if (type == TrackType.ALLELE_FREQUENCY) {
                track.setRendererClass(PointsRenderer.class);
                track.setHeight(40);
            }
            newTracks.add(track);
        }
    }

    private void loadAffectiveAnnotationTrack(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        AffectiveAnnotationParser parser = new AffectiveAnnotationParser();
        Map<String, List<Annotation>> annotMap = parser.parse(locator.getPath());
        AffectiveAnnotationTrack track = new AffectiveAnnotationTrack("id", "Annotations", annotMap);
        newTracks.add(track);

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

    private void loadWigFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        if (locator.isLocal()) {
            if (!checkSize(locator.getPath())) {
                return;
            }
        }

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


            DatasetDataSource dataSource = new DatasetDataSource(trackId, ds, genome);

            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            String displayName = (label == null || multiTrack) ? heading : label;
            track.setName(displayName);
            track.setProperties(props);

            track.setTrackType(ds.getType());

            if (ds.getType() == TrackType.EXPR) {
                track.setWindowFunction(WindowFunction.none);
            }


            newTracks.add(track);
        }
    }

    public void loadTDFFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {


        TDFReader reader = TDFReader.getReader(locator);
        TrackType type = reader.getTrackType();

        if (reader.getTrackType() == TrackType.AFFECTIVE) {
            AffectiveUtils.loadTDFFile(locator, newTracks, genome, reader);
            return;
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
            final DataSource dataSource = locator.getPath().endsWith(".counts") ?
                    new GobyCountArchiveDataSource(locator) :
                    new TDFDataSource(reader, trackNumber, heading, genome);
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            String displayName = (name == null || multiTrack) ? heading : name;
            track.setName(displayName);
            track.setTrackType(type);
            if (props != null) {
                track.setProperties(props);
            }
            newTracks.add(track);
            trackNumber++;
        }


    }

    public void loadBWFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        String trackName = locator.getTrackName();
        String trackId = locator.getPath();

        String path = locator.getPath();
        BBFileReader reader = new BBFileReader(path);
        BigWigDataSource bigwigSource = new BigWigDataSource(reader, genome);

        if (reader.isBigWigFile()) {
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, bigwigSource);
            newTracks.add(track);
        } else if (reader.isBigBedFile()) {

            if (locator.getPath().contains("RRBS_cpgMethylation") ||
                    locator.getPath().contains("BiSeq_cpgMethylation") ||
                    (reader.getAutoSql() != null && reader.getAutoSql().startsWith("table BisulfiteSeq"))) {
                loadMethylTrack(locator, reader, newTracks, genome);
            } else {
                FeatureTrack track = new FeatureTrack(locator, trackId, trackName, bigwigSource);
                newTracks.add(track);
            }
        } else {
            throw new RuntimeException("Unknown BIGWIG type: " + locator.getPath());
        }
    }

    private void loadMethylTrack(ResourceLocator locator, BBFileReader reader, List<Track> newTracks, Genome genome) throws IOException {

        MethylTrack track = new MethylTrack(locator, reader, genome);
        newTracks.add(track);
    }


    private void loadGobyCountsArchive(ResourceLocator locator, List<Track> newTracks, Genome genome) {


        if (log.isDebugEnabled()) {
            log.debug("Loading Goby counts archive: " + locator.toString());
        }


        String trackId = locator.getSampleId() + " coverage";
        String trackName = locator.getFileName();
        final DataSource dataSource = new GobyCountArchiveDataSource(locator);

        DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

        newTracks.add(track);


    }

    private void loadEwigIBFFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        TDFReader reader = TDFReader.getReader(locator.getPath());
        TrackProperties props = null;
        String trackLine = reader.getTrackLine();
        if (trackLine != null && trackLine.length() > 0) {
            props = new TrackProperties();
            ParsingUtils.parseTrackLine(trackLine, props);
        }

        EWigTrack track = new EWigTrack(locator, genome);
        if (props != null) {
            track.setProperties(props);
        }
        track.setName(locator.getTrackName());
        newTracks.add(track);
    }

    private void loadListFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {
        try {
            FeatureSource source = new FeatureDirSource(locator, genome);
            FeatureTrack track = new FeatureTrack(locator, source);
            track.setName(locator.getTrackName());
            track.setVisibilityWindow(0);
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
     *
     * @param locator
     * @param newTracks
     */
    private void loadRNAiGeneScoreFile(ResourceLocator locator,
                                       List<Track> newTracks, RNAIGeneScoreParser.Type type,
                                       Genome genome) {

        RNAIGeneScoreParser parser = new RNAIGeneScoreParser(locator.getPath(), type, genome);

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
            //track.setDataRange(new DataRange(-3, 0, 3));
            newTracks.add(track);
        }

    }

    /**
     * Load a RNAi haripin score file.  The results of this action are hairpin scores
     * added to the RNAIDataManager.  Currently no tracks are created for hairpin
     * scores, although this could change.
     *
     * @param locator
     */
    private void loadRNAiHPScoreFile(ResourceLocator locator) {
        (new RNAIHairpinParser(locator.getPath())).parse();
    }

    private void loadMultipleAlignmentTrack(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {
        MultipleAlignmentTrack t = new MultipleAlignmentTrack(locator, genome);
        t.setName("Multiple Alignments");
        newTracks.add(t);
    }

    private void loadPeakTrack(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {
        PeakTrack t = new PeakTrack(locator, genome);
        newTracks.add(t);
    }


    private void loadAlignmentsTrack(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        try {
            String dsName = locator.getTrackName();

            // If the user tried to load the index,  look for the file (this is a common mistake)
            if (locator.getPath().endsWith(".sai") || locator.getPath().endsWith(".bai")) {
                MessageUtils.showMessage("<html><b>ERROR:</b> Loading SAM/BAM index files are not supported:  " + locator.getPath() +
                        "<br>Load the SAM or BAM file directly. ");
                return;
            }

            AlignmentDataManager dataManager = new AlignmentDataManager(locator, genome);

            // Check that alignments we loaded actually match some data.  Many BAM files will contain some sequences
            // not represented in the genome, buf if there are no matches warn the user.
            List<String> seqNames = dataManager.getSequenceNames();
            if (seqNames != null && seqNames.size() > 0) {
                if (!checkSequenceNames(locator.getPath(), genome, seqNames)) {
                    return;
                }
            }

            if (locator.getPath().toLowerCase().endsWith(".bam")) {
                if (!dataManager.hasIndex()) {
                    MessageUtils.showMessage("<html>Could not load index file for: " +
                            locator.getPath() + "<br>  An index file is required for SAM & BAM files.");
                    return;
                }
            }

            AlignmentTrack alignmentTrack = new AlignmentTrack(locator, dataManager, genome);    // parser.loadTrack(locator, dsName);
            alignmentTrack.setName(dsName);


            // Create coverage track
            CoverageTrack covTrack = new CoverageTrack(locator, alignmentTrack.getName() + " Coverage", genome);
            covTrack.setVisible(PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_COV_TRACK));
            newTracks.add(covTrack);
            alignmentTrack.setCoverageTrack(covTrack);
            covTrack.setDataManager(dataManager);
            dataManager.setCoverageTrack(covTrack);

            // Search for precalculated coverage data   -- DON'T DO THIS FOR CGI READER
            String covPath = locator.getCoverage();
            if (covPath == null) {
                String path = locator.getPath();
                if (!path.contains("/query.cgi?")) {
                    covPath = path + ".tdf";
                }
            }
            if (covPath != null) {
                try {
                    if ((new File(covPath)).exists() || (HttpUtils.isRemoteURL(covPath) &&
                            HttpUtils.getInstance().resourceAvailable(new URL(covPath)))) {
                        TDFReader reader = TDFReader.getReader(covPath);
                        TDFDataSource ds = new TDFDataSource(reader, 0, alignmentTrack.getName() + " coverage", genome);
                        covTrack.setDataSource(ds);
                    }
                } catch (MalformedURLException e) {
                    // This is expected if
                    //    log.info("Could not loading coverage data: MalformedURL: " + covPath);
                }
            }

            boolean showSpliceJunctionTrack = PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_SHOW_JUNCTION_TRACK);
            if (showSpliceJunctionTrack) {
                SpliceJunctionFinderTrack spliceJunctionTrack = new SpliceJunctionFinderTrack(locator,
                        alignmentTrack.getName() + " Junctions", dataManager, genome);
                spliceJunctionTrack.setHeight(60);

                spliceJunctionTrack.setVisible(showSpliceJunctionTrack);
                newTracks.add(spliceJunctionTrack);
                alignmentTrack.setSpliceJunctionTrack(spliceJunctionTrack);
            }

            newTracks.add(alignmentTrack);

        } catch (IndexNotFoundException e) {
            MessageUtils.showMessage("<html>Could not find the index file for  <br><br>&nbsp;&nbsp;" + e.getSamFile() +
                    "<br><br>Note: The index file can be created using igvtools and must be in the same directory as the .sam file.");
        }
    }


    /**
     * Compare the sequence names against sequence (chromosome) names in the genome.  If no matches warn the user.
     *
     * @param filename
     * @param genome
     * @param seqNames
     * @return true if there is at least one sequence match, false otherwise
     */
    private boolean checkSequenceNames(String filename, Genome genome, List<String> seqNames) {
        boolean atLeastOneMatch = false;
        for (String seqName : seqNames) {
            if (genome.getChromosome(seqName) != null) {
                atLeastOneMatch = true;
                break;
            }
        }
        if (!atLeastOneMatch) {
            StringBuffer message = new StringBuffer();
            message.append("<html>File: " + filename +
                    "<br>does not contain any sequence names which match the current genome.");
            message.append("<br><br>File: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");
            int n = 0;
            for (String sn : seqNames) {
                message.append(sn + ", ");
                n++;
                if (n > 3) {
                    message.append(" ...");
                    break;
                }
            }
            message.append("<br>Genome: ");
            n = 0;
            for (String cn : genome.getAllChromosomeNames()) {
                message.append(cn + ", ");
                n++;
                if (n > 3) {
                    message.append(" ...");
                    break;
                }
            }
            MessageUtils.showMessage(message.toString());
        }
        return atLeastOneMatch;
    }


    /**
     * Load a mutation file (".mut" or ".maf").
     *
     * @param locator
     * @param newTracks
     */
    private void loadMutFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        MutationTrackLoader parser = new MutationTrackLoader();
        List<FeatureTrack> mutationTracks = parser.loadMutationTracks(locator, genome);
        for (FeatureTrack track : mutationTracks) {
            track.setTrackType(TrackType.MUTATION);
            track.setRendererClass(MutationRenderer.class);
            newTracks.add(track);
        }
    }

    private void loadFromDBProfile(ResourceLocator profileLocator, List<Track> newTracks) throws IOException {
        List<DBTable> tableList = DBTable.parseProfile(profileLocator.getPath());
        for (DBTable table : tableList) {
            SQLCodecSource source = SQLCodecSource.getFromTable(table);
            if (source != null) {
                CachingFeatureSource cachingReader = new CachingFeatureSource(source);
                FeatureTrack track = new FeatureTrack(profileLocator, cachingReader);
                track.setName(source.getTableName());
                newTracks.add(track);
            } else if (table.getFormat().equals("seg")) {
                Genome genome = GenomeManager.getInstance().getCurrentGenome();
                SegmentedAsciiDataSet ds = (new SegmentedReader(table.getDbLocator(), genome)).loadFromDB(table);
                loadSegTrack(table.getDbLocator(), newTracks, genome, ds);

            } else if (table.getFormat().equals("sample.info")) {
                //TODO sampleIdColumnLabel was previously hardcoded as "SAMPLE_ID_ARRAY"
                //TODO Basically I'm shoehorning this information into a field usually used for something else. Only slightly better
                String sampleIdColumnLabel = table.getBinColName();
                if (sampleIdColumnLabel == null) {
                    throw new IllegalArgumentException("Profile must have binColName specifying the sample id column label");
                }
                (new SampleInfoSQLReader(table, sampleIdColumnLabel)).load();
            }
        }

    }


    /**
     * @param locator
     * @param newTracks
     * @param genome
     * @deprecated See loadFromDBProfile, which loads from an xml file specifying table characteristics
     */
    private void loadFromDatabase(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        if (".seg".equals(locator.getType())) {

            //TODO Don't hardcode table name, this might note even be right for our target case
            SegmentedAsciiDataSet ds = (new SegmentedSQLReader(locator, "CNV", genome)).load();
            loadSegTrack(locator, newTracks, genome, ds);
        } else {
            (new SampleInfoSQLReader(locator, "SAMPLE_INFO", "SAMPLE_ID_ARRAY")).load();
        }
    }

    private void loadSegFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        // TODO - -handle remote resource
        SegmentedDataSet ds;
        String path = locator.getPath().toLowerCase();

        if (path.endsWith("seg.zip")) {
            ds = new SegmentedBinaryDataSet(locator);
        } else {
            SegmentFileParser parser = new SegmentFileParser(locator);
            ds = parser.loadSegments(locator, genome);
        }

        loadSegTrack(locator, newTracks, genome, ds);
    }

    /**
     * Add the provided SegmentedDataSet to the list of tracks,
     * set other relevant properties
     *
     * @param locator
     * @param newTracks
     * @param genome
     * @param ds
     */
    private void loadSegTrack(ResourceLocator locator, List<Track> newTracks, Genome genome, SegmentedDataSet ds) {
        String path = locator.getPath();

        TrackProperties props = null;
        if (ds instanceof SegmentedAsciiDataSet) {
            props = ((SegmentedAsciiDataSet) ds).getTrackProperties();
        }

        // The "freq" track.  TODO - make this optional
        if ((ds.getType() == TrackType.COPY_NUMBER || ds.getType() == TrackType.CNV) &&
                ds.getSampleNames().size() > 4) {
            FreqData fd = new FreqData(ds, genome);
            String freqTrackId = path;
            String freqTrackName = "CNV Summary";
            CNFreqTrack freqTrack = new CNFreqTrack(locator, freqTrackId, freqTrackName, fd);
            newTracks.add(freqTrack);
        }

        for (String trackName : ds.getSampleNames()) {
            String trackId = path + "_" + trackName;
            SegmentedDataSource dataSource = new SegmentedDataSource(trackName, ds);
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);
            track.setRendererClass(HeatmapRenderer.class);
            track.setTrackType(ds.getType());

            if (props != null) {
                track.setProperties(props);
            }

            newTracks.add(track);
        }
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
        String name = locator.getName();
        if (name == null || name.length() == 0) {
            if (locator.getPath().contains("genome.ucsc.edu")) {
                name = featureSource.getType();
            } else {
                name = featureSource.getPath().replace("/das/", "").replace("/features", "");
            }
        }
        track.setName(name);

        // A hack until we can notate this some other way
        if (locator.getPath().contains("cosmic")) {
            track.setRendererClass(CosmicFeatureRenderer.class);
            track.setMinimumHeight(2);
            track.setHeight(20);
            track.setDisplayMode(Track.DisplayMode.EXPANDED);
        } else {
            track.setRendererClass(IGVFeatureRenderer.class);
            track.setMinimumHeight(35);
            track.setHeight(45);
        }
        currentTracks.add(track);
    }


    private void loadTrioData(ResourceLocator locator) throws IOException {
        PedigreeUtils.parseTrioFile(locator.getPath());
    }


    public static boolean isIndexed(String path, Genome genome) {

        // Checking for the index is expensive over HTTP.  First see if this is an indexable format by fetching the codec
        if (!isIndexable(path, genome)) {
            return false;
        }

        String indexExtension = path.endsWith("gz") ? ".tbi" : ".idx";
        String indexPath = path + indexExtension;
        try {
            if (HttpUtils.isRemoteURL(path)) {
                return HttpUtils.getInstance().resourceAvailable(new URL(indexPath));
            } else {
                File f = new File(path + indexExtension);
                return f.exists();
            }

        } catch (IOException e) {
            return false;
        }

    }


    /**
     * Return true if a file represented by "path" is indexable.  This method is an optimization, we could just look
     * for the index but that is expensive to do for remote resources.  All tribble indexable extensions should be
     * listed here.
     *
     * @param path
     * @return
     */
    private static boolean isIndexable(String path, Genome genome) {
        String fn = path.toLowerCase();
        if (fn.endsWith(".gz")) {
            int l = fn.length() - 3;
            fn = fn.substring(0, l);
        }
        // The vcf extension is for performance, it doesn't matter which codec is returned all vcf files
        // are indexable.
        return path.endsWith(".vcf") || CodecFactory.getCodec(path, genome) != null;
    }

    public static TrackProperties getTrackProperties(Object header) {
        try {
            FeatureFileHeader ffHeader = (FeatureFileHeader) header;
            if (ffHeader != null) {
                return ffHeader.getTrackProperties();
            } else {
                return null;
            }
        } catch (ClassCastException e) {
            return null;
        }
    }


}
