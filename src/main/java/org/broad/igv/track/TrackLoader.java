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

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.log4j.Logger;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bigwig.BigWigDataSource;
import org.broad.igv.blast.BlastMapping;
import org.broad.igv.blast.BlastParser;
import org.broad.igv.data.*;
import org.broad.igv.data.cufflinks.*;
import org.broad.igv.data.expression.ExpressionDataset;
import org.broad.igv.data.expression.ExpressionFileParser;
import org.broad.igv.data.seg.*;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.BasePairFileUtils;
import org.broad.igv.feature.GisticFileParser;
import org.broad.igv.feature.MutationTrackLoader;
import org.broad.igv.feature.ShapeFileUtils;
import org.broad.igv.feature.basepair.BasePairTrack;
import org.broad.igv.bedpe.BedPEParser;
import org.broad.igv.bedpe.InteractionTrack;
import org.broad.igv.feature.bionano.SMAPParser;
import org.broad.igv.feature.bionano.SMAPRenderer;
import org.broad.igv.feature.dranger.DRangerParser;
import org.broad.igv.feature.dsi.DSIRenderer;
import org.broad.igv.feature.dsi.DSITrack;
import org.broad.igv.feature.genome.load.GenbankParser;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.gff.GFFFeatureSource;
import org.broad.igv.feature.sprite.ClusterParser;
import org.broad.igv.feature.sprite.ClusterTrack;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.feature.tribble.GFFCodec;
import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.goby.GobyAlignmentQueryReader;
import org.broad.igv.goby.GobyCountArchiveDataSource;
import org.broad.igv.google.Ga4ghAPIHelper;
import org.broad.igv.google.GoogleUtils;
import org.broad.igv.gwas.*;
import org.broad.igv.lists.GeneList;
import org.broad.igv.lists.GeneListManager;
import org.broad.igv.maf.MultipleAlignmentTrack;
import org.broad.igv.methyl.MethylTrack;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.HeatmapRenderer;
import org.broad.igv.renderer.MutationRenderer;
import org.broad.igv.renderer.PointsRenderer;
import org.broad.igv.sam.*;
import org.broad.igv.sam.reader.IndexNotFoundException;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.ConfirmDialog;
import org.broad.igv.ui.util.ConvertFileDialog;
import org.broad.igv.ui.util.ConvertOptions;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.util.PedigreeUtils;

import java.io.IOException;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * User: jrobinso
 * Date: Feb 14, 2010
 */
public class TrackLoader {

    private static Logger log = Logger.getLogger(TrackLoader.class);

    private static Collection<? extends Class> NOLogExceptions = Arrays.asList(TribbleIndexNotFoundException.class);

    /**
     * Switches on various attributes of locator (mainly locator path extension and whether the locator is indexed)
     * to call the appropriate loading method.
     *
     * @param locator
     * @param genome
     * @return
     */
    public List<Track> load(ResourceLocator locator, Genome genome) throws DataLoadException {

        final String path = locator.getPath().trim();
        final boolean isHtsGet = locator.getAttribute("htsget") != null && (Boolean) locator.getAttribute("htsget");

        // Check if the AWS credentials are still valid. If not, re-login and renew pre-signed urls
        if (AmazonUtils.isAwsS3Path(path)) {
            AmazonUtils.checkLogin();
        }

        log.info("Loading resource, path " + path);
        try {
            String typeString = locator.getTypeString();

            if (typeString.endsWith(".tbi")) {
                MessageUtils.showMessage("<html><b>Error:</b>File type '.tbi' is not recognized.  If this is a 'tabix' index <br>" +
                        " load the associated gzipped file, which should have an extension of '.gz'");
            }

            //This list will hold all new tracks created for this locator
            List<Track> newTracks = new ArrayList<Track>();

            if (typeString.endsWith(".gmt")) {
                loadGMT(locator);
            } else if (typeString.endsWith(".vcf.list")) {
                loadVCFListFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".trio")) {
                loadTrioData(locator);
            } else if (typeString.endsWith(".gct") || typeString.endsWith("res") || typeString.endsWith("tab")) {
                loadGctFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".gbk") || typeString.endsWith(".gb")) {
                loadGbkFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".cn") || typeString.endsWith(".xcn") || typeString.endsWith(".snp") ||
                    typeString.endsWith(".igv") || typeString.endsWith(".loh")) {
                loadIGVFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".cbs") || typeString.endsWith(".seg") ||
                    typeString.endsWith("glad") || typeString.endsWith("birdseye_canary_calls")
                    || typeString.endsWith(".seg.zip")) {
                loadSegFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".gistic")) {
                loadGisticFile(locator, newTracks);
            } else if (typeString.contains(".tabblastn") || typeString.endsWith(".orthologs")) {
                loadBlastMapping(locator, newTracks);
            } else if (isAlignmentTrack(typeString) || (path.startsWith("http") && path.contains("/query.cgi?")) || isHtsGet) {
                loadAlignmentsTrack(locator, newTracks, genome);
            } else if (typeString.endsWith(".shape") || typeString.endsWith(".map")) {
                convertLoadShapeFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".wig") || typeString.endsWith(".bedgraph") || typeString.endsWith(".bdg") ||
                    typeString.endsWith("cpg.txt") || typeString.endsWith(".expr")) {
                loadWigFile(locator, newTracks, genome);
            } else if (typeString.endsWith("fpkm_tracking") || typeString.endsWith("gene_exp.diff") ||
                    typeString.endsWith("cds_exp.diff")) {
                loadCufflinksFile(locator, newTracks, genome);
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
            } else if (WiggleParser.isWiggle(locator)) {
                loadWigFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".maf.dict")) {
                loadMultipleAlignmentTrack(locator, newTracks, genome);
            } else if (typeString.endsWith("mage-tab") || ExpressionFileParser.parsableMAGE_TAB(locator)) {
                locator.setDescription("MAGE_TAB");
                loadGctFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".db") || typeString.endsWith(".dbn")) {
                convertLoadStructureFile(locator, newTracks, genome, "dotBracket");
            } else if (typeString.endsWith(".ct")) {
                convertLoadStructureFile(locator, newTracks, genome, "connectTable");
            } else if (typeString.endsWith(".dp")) {
                convertLoadStructureFile(locator, newTracks, genome, "pairingProb");
            } else if (typeString.endsWith(".bp")) {
                loadBasePairFile(locator, newTracks, genome);
            } else if (GWASParser.isGWASFile(typeString)) {
                loadGWASFile(locator, newTracks, genome);
            } else if (GobyAlignmentQueryReader.supportsFileType(path)) {
                loadAlignmentsTrack(locator, newTracks, genome);
            } else if (typeString.endsWith(".list")) {
                // This should be deprecated
                loadListFile(locator, newTracks, genome);
            } else if (typeString.endsWith(".smap")) {
                loadSMAPFile(locator, newTracks, genome);
            } else if (typeString.endsWith("dsi")) {
                loadDSIFile(locator, newTracks, genome);
            } else if (typeString.endsWith("bedpe") || typeString.endsWith("_clusters")) {
                loadBedPEFile(locator, newTracks, genome);
            } else if (typeString.endsWith("clusters")) {
                loadClusterFile(locator, newTracks, genome);
            } else if (CodecFactory.hasCodec(locator, genome) && !forceNotTribble(typeString)) {
                loadTribbleFile(locator, newTracks, genome);
            } else if (MutationTrackLoader.isMutationAnnotationFile(locator)) {
                loadMutFile(locator, newTracks, genome); // Must be tried before ".maf" test below
            } else if (typeString.endsWith(".maf")) {
                loadMultipleAlignmentTrack(locator, newTracks, genome);
            } else if (AttributeManager.isSampleInfoFile(locator)) {
                // This might be a sample information file.
                AttributeManager.getInstance().loadSampleInfo(locator);
            } else {
                MessageUtils.showMessage("<html>Unknown file type: " + path + "<br>Check file extension");
            }


            // Track line
            if (newTracks.size() > 0) {
                TrackProperties tp = null;
                String trackLine = locator.getTrackLine();
                if (trackLine != null) {
                    tp = new TrackProperties();
                    ParsingUtils.parseTrackLine(trackLine, tp);
                }

                for (Track track : newTracks) {
                    if (locator.getFeatureInfoURL() != null) {
                        track.setUrl(locator.getFeatureInfoURL());
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
            }

            return newTracks;
        } catch (Exception e) {
            if (!NOLogExceptions.contains(e.getClass())) {
                log.error(e.getMessage(), e);
            }
            throw new DataLoadException(e.getMessage());
        }

    }

    public static boolean isAlignmentTrack(String typeString) {
        return typeString.endsWith(".sam") || typeString.endsWith(".bam") || typeString.endsWith(".cram") ||
                typeString.endsWith(".sam.list") || typeString.endsWith(".bam.list") ||
                typeString.endsWith(".aligned") || typeString.endsWith(".sai") ||
                typeString.endsWith(".bai") || typeString.endsWith(".csi") || typeString.equals("alist") ||
                typeString.equals(Ga4ghAPIHelper.RESOURCE_TYPE);
    }

    private void loadSMAPFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        List<Feature> features = SMAPParser.parseFeatures(locator, genome);
        FeatureCollectionSource src = new FeatureCollectionSource(features, genome);
        FeatureTrack track = new FeatureTrack(locator, locator.getName(), src);
        track.setRendererClass(SMAPRenderer.class);
        track.setDisplayMode(Track.DisplayMode.EXPANDED);
        newTracks.add(track);
    }

    private boolean forceNotTribble(String typeString) {
        List<String> nonTribble = Arrays.asList("fpkm_tracking", "exp_diff", "_exp.diff");
        for (String s : nonTribble) {
            if (typeString.endsWith(s)) {
                return true;
            }
        }
        return false;
    }


    private void loadGMT(ResourceLocator locator) throws IOException {
        List<GeneList> lists = GeneListManager.getInstance().loadGMTFile(locator.getPath());
        if (lists.size() == 1) {
            GeneList gl = lists.get(0);
            IGV.getInstance().setGeneList(gl, true);
        } else {
            MessageUtils.showMessage("Loaded " + lists.size() + " gene lists.");
        }
    }

    private void loadVCF(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException, TribbleIndexNotFoundException {


        TribbleFeatureSource src = TribbleFeatureSource.getFeatureSource(locator, genome);


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

    private void loadVCFListFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException, TribbleIndexNotFoundException {

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


    private void loadBlastMapping(ResourceLocator locator, List<Track> newTracks) {

        List<BlastMapping> mappings = (new BlastParser()).parse(locator.getPath());
        List<htsjdk.tribble.Feature> features = new ArrayList<htsjdk.tribble.Feature>(mappings.size());
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


    private void loadBedPEFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {
        BedPEParser.Dataset features = BedPEParser.parse(locator, genome);
        newTracks.add(new InteractionTrack(locator, features, genome));
    }

    private void loadClusterFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {
        ClusterParser.ClusterSet features = ClusterParser.parse(locator.getPath());
        newTracks.add(new ClusterTrack(locator, features, genome));
    }


    /**
     * Load the input file as a feature, mutation, or maf (multiple alignment) file.
     *
     * @param locator
     * @param newTracks
     */
    private void loadTribbleFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException, TribbleIndexNotFoundException {

        String typeString = locator.getTypeString();

        // Mutation (mut, maf, vcf) files are handled special.  Check here, rather than depend on order in giant case statement.
        if (MutationTrackLoader.isMutationAnnotationFile(locator)) {
            loadMutFile(locator, newTracks, genome); // Must be tried before generic "loadIndexed" below
        } else if (VariantTrack.isVCF(typeString)) {
            loadVCF(locator, newTracks, genome);
        } else {

            TribbleFeatureSource tribbleFeatureSource = TribbleFeatureSource.getFeatureSource(locator, genome);

            FeatureSource src;
            if(GFFFeatureSource.isGFF(locator.getPath())) {
                GFFCodec codec =  (GFFCodec) CodecFactory.getCodec(locator, genome);
                 src = new GFFFeatureSource(tribbleFeatureSource, codec.getVersion()) ;
            } else {
                src = tribbleFeatureSource;
            }

            // Create feature source and track
            FeatureTrack t = new FeatureTrack(locator, src);
            t.setName(locator.getTrackName());
            //t.setRendererClass(BasicTribbleRenderer.class);

            // Set track properties from header
            Object header = tribbleFeatureSource.getHeader();
            if (header != null && header instanceof FeatureFileHeader) {
                FeatureFileHeader ffh = (FeatureFileHeader) header;
                if (ffh.getTrackType() != null) {
                    t.setTrackType(ffh.getTrackType());
                }
                if (ffh.getTrackProperties() != null) {
                    TrackProperties tp = ffh.getTrackProperties();
                    t.setProperties(tp);
                    t.setTrackLine(tp.getTrackLine());
                }
                if (ffh.getTrackType() == TrackType.REPMASK) {
                    t.setHeight(15);
                }
            }
            if (locator.getPath().contains(".narrowPeak") ||
                    locator.getPath().contains(".broadPeak") ||
                    locator.getPath().contains(".gappedPeak")||
                    locator.getPath().contains(".regionPeak") ) {
                t.setUseScore(true);
            }
            newTracks.add(t);
        }
    }

    private void loadDSIFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException, TribbleIndexNotFoundException {

        TribbleFeatureSource tribbleFeatureSource = TribbleFeatureSource.getFeatureSource(locator, genome);

        // Create feature source and track
        DSITrack t = new DSITrack(locator, tribbleFeatureSource);
        t.setName(locator.getTrackName());
        //t.setRendererClass(BasicTribbleRenderer.class);

        // Set track properties from header
        Object header = tribbleFeatureSource.getHeader();
        if (header != null && header instanceof TrackProperties) {
            TrackProperties tp = (TrackProperties) header;
            t.setProperties(tp);
            t.setTrackLine(tp.getTrackLine());
        }

        t.setRendererClass(DSIRenderer.class);

        newTracks.add(t);
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
        Map<String, List<GWASFeature>> gwasData = gwasParser.parse();

        GWASTrack gwasTrack = new GWASTrack(locator, locator.getPath(), locator.getFileName(), gwasData, gwasParser.getColumnHeaders(), genome);
        newTracks.add(gwasTrack);

    }


    private void loadGctFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        if (locator.isLocal()) {
            if (!checkSize(locator)) {
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

    /**
     * Load features from a genbank (.gbk)file.  This method ignores the fasta section.  To define a genome from
     * a genbank file use GenomeManager.
     *
     * @param newTracks
     * @param genome
     * @throws IOException
     */
    private void loadGbkFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        GenbankParser genbankParser = new GenbankParser(locator.getPath());
        genbankParser.readFeatures(false);
        FeatureCollectionSource src = new FeatureCollectionSource(genbankParser.getFeatures(), genome);
        FeatureTrack track = new FeatureTrack(locator, src);
        newTracks.add(track);
    }

    private void loadIGVFile(ResourceLocator locator, List<Track> newTracks, Genome genome) {

        if (locator.isLocal()) {
            if (!checkSize(locator)) {
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


    private void loadCufflinksFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        final String path = locator.getPath();
        final String s = path.toLowerCase();
        List<DataTrack> cuffTracks = new ArrayList<DataTrack>();
        if (s.endsWith("fpkm_tracking")) {
            FPKMTrackingCodec codec = new FPKMTrackingCodec(path);
            List<FPKMValue> values = CufflinksParser.parse(codec, path);
            for (int sampleIndex = 0; sampleIndex < codec.getNumSamples(); sampleIndex++) {
                CufflinksDataSource ds = new CufflinksDataSource(sampleIndex, values, genome);
                String supId = String.format("q%02d", sampleIndex);
                DataTrack track = new DataSourceTrack(locator, locator.getPath() + " " + supId, locator.getTrackName() + " " + supId, ds);
                cuffTracks.add(track);
            }
        } else if (s.endsWith("gene_exp.diff") || s.endsWith("cds_exp.diff")) {
            AsciiFeatureCodec<ExpDiffValue> codec = new ExpDiffCodec(path);
            List<ExpDiffValue> values = CufflinksParser.parse(codec, path);
            CufflinksDataSource ds = new CufflinksDataSource(values, genome);
            DataTrack track = new DataSourceTrack(locator, locator.getPath(), locator.getTrackName(), ds);
            cuffTracks.add(track);
        } else {
            throw new RuntimeException("Unsupported file type: " + path);
        }

        for (DataTrack track : cuffTracks) {
            track.setTrackType(TrackType.FPKM);
            CufflinksTrack.setCufflinksScale(track);
            newTracks.add(track);
        }
    }


    private static boolean checkSize(ResourceLocator locator) {

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SIZE_WARNING)) {
            return true;
        }

        final String path = locator.getPath();
        long size = FileUtils.getLength(path);
        int maxSize = 200000000;  // 200 mb
        if (path.endsWith(".gz") || path.endsWith(".bgz")) {
            maxSize /= 4;
        }

        if (size > maxSize) {

            String message = "The file " + path + " is large (" + (size / 1000000) + " mb).  It is recommended " +
                    "that large files be converted to the binary <i>.tdf</i> format using the IGVTools " +
                    "<b>toTDF</b> command. Loading  unconverted ascii fies of this size can lead to poor " +
                    "performance or unresponsiveness (freezing).  " +
                    "<br><br>IGVTools can be launched from the <b>Tools</b> menu or separately as a " +
                    "command line program. See the user guide for more details.<br><br>Click <b>Continue</b> " +
                    "to continue loading, or <b>Cancel</b> to skip this file.";

            return ConfirmDialog.optionallyShowConfirmDialog(message, SHOW_SIZE_WARNING, true);


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
            if (!checkSize(locator)) {
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

        log.debug("Loading TDF file " + locator.getPath());
        TDFReader reader = TDFReader.getReader(locator);
        TrackType type = reader.getTrackType();

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

    private void loadMultipleAlignmentTrack(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {
        MultipleAlignmentTrack t = new MultipleAlignmentTrack(locator, genome);
        t.setName("Multiple Alignments");
        newTracks.add(t);
    }

    private void loadAlignmentsTrack(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException {

        try {
            String dsName = locator.getTrackName();

            // If the user tried to load the index,  look for the file (this is a common mistake)
            if (locator.getTypeString().endsWith(".sai") ||
                    locator.getTypeString().endsWith(".bai") ||
                    locator.getTypeString().endsWith(".csi")) {
                MessageUtils.showMessage("<html><b>ERROR:</b> Loading SAM/BAM index files are not supported:  " + locator.getPath() +
                        "<br>Load the SAM or BAM file directly. ");
                return;
            }

            AlignmentDataManager dataManager = new AlignmentDataManager(locator, genome);

            // Check that alignments we loaded actually match some data.  Many BAM files will contain some sequences
            // not represented in the genome, buf if there are no matches warn the user.
            List<String> seqNames = dataManager.getSequenceNames();
            if (seqNames != null && seqNames.size() > 0) {
                if (!dataManager.hasMatchingSequences()) {
                    showMismatchSequenceNameMessage(locator.getPath(), genome, seqNames);
                }
            }

            if (locator.getTypeString().endsWith("bam") || locator.getTypeString().endsWith("cram")) {
                if (!dataManager.hasIndex()) {
                    MessageUtils.showMessage("<html>Could not load index file for: " +
                            locator.getPath() + "<br>  An index file is required for SAM & BAM files.");
                    return;
                }
            }

            AlignmentTrack alignmentTrack = new AlignmentTrack(locator, dataManager, genome);    // parser.loadTrack(locator, dsName)
            alignmentTrack.setName(dsName);
            alignmentTrack.setVisible(PreferencesManager.getPreferences().getAsBoolean(SAM_SHOW_ALIGNMENT_TRACK));

            // Create coverage track
            CoverageTrack covTrack = new CoverageTrack(locator, dsName + " Coverage", alignmentTrack, genome);
            covTrack.setVisible(PreferencesManager.getPreferences().getAsBoolean(SAM_SHOW_COV_TRACK));
            newTracks.add(covTrack);
            covTrack.setDataManager(dataManager);
            dataManager.setCoverageTrack(covTrack);

            alignmentTrack.setCoverageTrack(covTrack);

            // Search for precalculated coverage data
            // Skip for GA4GH & SU2C resources
            if (!(Ga4ghAPIHelper.RESOURCE_TYPE.equals(locator.getType()) ||
                    locator.getPath().contains("dataformat=.bam") ||
                    GoogleUtils.isGoogleCloud(locator.getPath()))) {

                String covPath = locator.getCoverage();
                if (covPath == null) {
                    boolean bypassFileAutoDiscovery = PreferencesManager.getPreferences().getAsBoolean(BYPASS_FILE_AUTO_DISCOVERY);
                    String path = locator.getPath();
                    if (!bypassFileAutoDiscovery && !path.contains("/query.cgi?")) {
                        covPath = path + ".tdf";
                    }

                }
                if (covPath != null && !covPath.equals(".")) {
                    if (FileUtils.resourceExists(covPath)) {
                        log.debug("Loading TDF for coverage: " + covPath);
                        try {
                            TDFReader reader = TDFReader.getReader(covPath);
                            TDFDataSource ds = new TDFDataSource(reader, 0, dsName + " coverage", genome);
                            covTrack.setDataSource(ds);
                        } catch (Exception e) {
                            log.error("Error loading coverage TDF file", e);
                        }
                    }

                }
            }

            boolean showSpliceJunctionTrack = PreferencesManager.getPreferences().getAsBoolean(SAM_SHOW_JUNCTION_TRACK);

            SpliceJunctionTrack spliceJunctionTrack = new SpliceJunctionTrack(locator,
                    dsName + " Junctions", dataManager, alignmentTrack, SpliceJunctionTrack.StrandOption.BOTH);
            spliceJunctionTrack.setHeight(60);
            spliceJunctionTrack.setVisible(showSpliceJunctionTrack);
            newTracks.add(spliceJunctionTrack);

            alignmentTrack.setSpliceJunctionTrack(spliceJunctionTrack);

            newTracks.add(alignmentTrack);

            log.debug("Alignment track loaded");

        } catch (IndexNotFoundException e) {
            MessageUtils.showMessage("<html>Could not find the index file for  <br><br>&nbsp;&nbsp;" + e.getSamFile() +
                    "<br><br>Note: The index file can be created using igvtools and must be in the same directory as the .sam file.");
        }
    }


    private void showMismatchSequenceNameMessage(String filename, Genome genome, List<String> seqNames) {
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


    /**
     * Load a mutation file (".mut" or ".maf").
     *
     * @param locator
     * @param newTracks
     */
    private void loadMutFile(ResourceLocator locator, List<Track> newTracks, Genome genome) throws IOException, TribbleIndexNotFoundException {

        MutationTrackLoader loader = new MutationTrackLoader();
        List<FeatureTrack> mutationTracks = loader.loadMutationTracks(locator, genome);
        for (FeatureTrack track : mutationTracks) {
            track.setTrackType(TrackType.MUTATION);
            track.setRendererClass(MutationRenderer.class);
            newTracks.add(track);
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
                ds.getSampleNames().size() > 1) {
            FreqData fd = new FreqData(ds, genome);
            String freqTrackId = path;
            String freqTrackName = "CNV Summary";
            CNFreqTrack freqTrack = new CNFreqTrack(locator, freqTrackId, freqTrackName, fd);

            if (props != null) {
                freqTrack.setProperties(props);
            }

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

    private void loadTrioData(ResourceLocator locator) throws IOException {
        PedigreeUtils.parseTrioFile(locator.getPath());
    }

    /**
     * Convert an RNA chemical reactivity file (.shape, .map) into a .wig file
     * and load.
     */
    private void convertLoadShapeFile(ResourceLocator locator,
                                      List<Track> newTracks,
                                      Genome genome) throws IOException {
        String inPath = locator.getPath();
        String fileName = locator.getFileName();
        String outPath = inPath + ".wig";
        String message = "The chemical reactivity file <br> &nbsp;&nbsp;" + fileName + "<br> needs to be converted to IGV chromosome <br>" +
                "coordinates and .wig format before loading. <br><br>Click <b>Continue</b> " +
                "to save converted file to <br> &nbsp;&nbsp;" + fileName + ".wig" +
                "<br>and load with the selected options, or <b>Cancel</b> to skip this<br>file.";

        ConvertOptions opts = ConvertFileDialog.showConvertFileDialog(message);
        if (opts.doConvert) {
            ShapeFileUtils.shapeToWigFile(inPath, outPath, opts.chrom, opts.strand, opts.start);
            loadWigFile(new ResourceLocator(outPath), newTracks, genome);
        }
    }

    /**
     * Convert various RNA structure formats to a more easily parseable format
     * in genomic coordinates, then load converted file.
     */
    private void convertLoadStructureFile(ResourceLocator locator,
                                          List<Track> newTracks,
                                          Genome genome,
                                          String fileType) throws IOException {
        String inPath = locator.getPath();
        String fileName = locator.getFileName();
        String outPath = inPath + ".bp";

        String message = "The RNA structure file <br> &nbsp;&nbsp;" + fileName + "<br> needs to be converted to IGV chromosome <br>" +
                "coordinates and .bp format before loading. <br><br>Click <b>Continue</b> " +
                "to save converted file to <br> &nbsp;&nbsp;" + fileName + ".bp" +
                "<br>and load with the selected options, or <b>Cancel</b> to skip this<br>file.";

        ConvertOptions opts = ConvertFileDialog.showConvertFileDialog(message);

        if (opts.doConvert) {
            if (fileType == "connectTable") {
                BasePairFileUtils.connectTableToBasePairFile(inPath, outPath, opts.chrom, opts.strand, opts.start);
            } else if (fileType == "pairingProb") {
                BasePairFileUtils.pairingProbToBasePairFile(inPath, outPath, opts.chrom, opts.strand, opts.start);
            } else if (fileType == "dotBracket") {
                BasePairFileUtils.dotBracketToBasePairFile(inPath, outPath, opts.chrom, opts.strand, opts.start);
            }
            loadBasePairFile(new ResourceLocator(outPath), newTracks, genome);
        }
    }

    private void loadBasePairFile(ResourceLocator locator,
                                  List<Track> newTracks,
                                  Genome genome) throws IOException {
        String name = locator.getTrackName();
        String path = locator.getPath();
        String id = path + "_" + name;
        newTracks.add(new BasePairTrack(locator, id, name, genome));
    }

    public static boolean isIndexed(ResourceLocator locator, Genome genome) {

        // Checking for the index is expensive over HTTP.  First see if this is an indexable format by fetching the codec
        String fullPath = locator.getPath();
        String pathNoQuery = locator.getURLPath();
        if (!CodecFactory.hasCodec(locator, genome)) {
            return false;
        }

        String indexExtension = pathNoQuery.endsWith("gz") ? ".tbi" : ".idx";

        String indexPath = fullPath + indexExtension;
        if (HttpUtils.isRemoteURL(fullPath)) {
            //Handle query string, if it exists
            String[] toks = fullPath.split("\\?", 2);
            if (toks.length == 2) {
                indexPath = String.format("%s%s?%s", toks[0], indexExtension, toks[1]);
            }
        }
        return FileUtils.resourceExists(indexPath);

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
