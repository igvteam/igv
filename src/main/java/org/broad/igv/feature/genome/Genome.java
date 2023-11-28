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


/*
 * Genome.java
 *
 * Created on November 9, 2007, 9:05 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.NamedFeature;
import org.apache.commons.math3.stat.StatUtils;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Cytoband;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndex;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.feature.genome.load.ChromSizesParser;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.ucsc.twobit.TwoBitSequence;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.liftover.Liftover;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Simple model of a genome.  Keeps an ordered list of Chromosomes, an alias table, and genome position offsets
 * for each chromosome to support a whole-genome view.
 */
public class Genome {

    private static final int MAX_WHOLE_GENOME_LONG = 100;
    private static Logger log = LogManager.getLogger(Genome.class);
    private String id;
    private String displayName;
    private List<String> chromosomeNames;
    private List<String> longChromosomeNames;
    private LinkedHashMap<String, Chromosome> chromosomeMap;
    private long totalLength = -1;
    private long nominalLength = -1;
    private Map<String, Long> cumulativeOffsets = new HashMap();
    private Map<String, String> chrAliasCache = new HashMap<>();
    private ChromAliasSource chromAliasSource;
    private Sequence sequence;
    private FeatureTrack geneTrack;
    private String species;
    private String ucscID;
    private String blatDB;
    private List<ResourceLocator> annotationResources;
    private Map<ResourceLocator, List<Track>> annotationTracks;
    private boolean showWholeGenomeView = true;
    private Map<String, Liftover> liftoverMap;
    private CytobandSource cytobandSource;
    private String homeChromosome;
    private String defaultPos;

    private String nameSet;

    public Genome(GenomeConfig config) throws IOException {

        id = config.id;
        displayName = config.name;
        nameSet = config.nameSet;
        blatDB = config.blatDB;
        if (config.ucsdID == null) {
            ucscID = ucsdIDMap.containsKey(id) ? ucsdIDMap.get(id) : id;
        } else {
            ucscID = config.ucsdID;
        }
        blatDB = (config.blatDB != null) ? config.blatDB : ucscID;
        defaultPos = config.defaultPos;

        // Load the sequence object.  Some configurations will specify both 2bit and fasta references.  The 2 bit
        // has preference
        if (config.sequence != null) {
            sequence = config.sequence;
        } else if (config.twoBitURL != null) {
            sequence = (config.twoBitBptURL != null) ?
                    new TwoBitSequence(config.twoBitURL, config.twoBitBptURL) :
                    new TwoBitSequence(config.twoBitURL);
        } else if (config.fastaURL != null) {
            String fastaPath = config.fastaURL;
            String indexPath = config.indexURL;
            String gziIndexPath = config.gziIndexURL;
            sequence = fastaPath.endsWith(".gz") ?
                    new FastaBlockCompressedSequence(fastaPath, gziIndexPath, indexPath) :
                    new FastaIndexedSequence(fastaPath, indexPath);
        } else {
            throw new RuntimeException("Genomes require either a .2bit or fasta reference ");
        }

        // Search for chromosomes.  Chromosome names are required to support the chromosome pulldown, names and
        // lengths are required to support whole genome view.  Both can be obtained from fasta index files, but
        // for .2bit sequences a 'chromSizes" file is required.  If not supplied the chr pulldown and wg view are disabled.
        List<Chromosome> chromosomeList = null;
        if (config.chromSizesURL != null) {
            chromosomeList = ChromSizesParser.parse(config.chromSizesURL);
        } else if (sequence.hasChromosomes()) {
            chromosomeList = sequence.getChromosomes();
        } else if (config.indexURL != null) {
            FastaIndex index = new FastaIndex(config.indexURL);
            chromosomeList = index.getChromosomes();
        }

        // Whole genome view is initially enabled by default if we have the chromosome information.  This might
        // still get disabled if we can't determine a sufficiently small number of wg chromosomes.
        showWholeGenomeView = (config.wholeGenomeView == null || config.wholeGenomeView) && chromosomeList.size() > 1;

        this.chromosomeMap = new LinkedHashMap<>();
        this.chromosomeNames = new ArrayList<>();
        if (chromosomeList != null) {

            chromosomeList = sortChromosomeList(chromosomeList, config.chromosomeOrder);

            for (Chromosome c : chromosomeList) {
                this.chromosomeMap.put(c.getName(), c);
                this.chromosomeNames.add(c.getName());
            }

            if (config.chromosomeOrder != null) {
                this.longChromosomeNames = Arrays.asList(config.chromosomeOrder);
            } else {
                this.longChromosomeNames = computeLongChromosomeNames();
                if (longChromosomeNames.size() > MAX_WHOLE_GENOME_LONG) {
                    showWholeGenomeView = false;
                }
            }
        }


        if (showWholeGenomeView) {
            homeChromosome = Globals.CHR_ALL;
        } else if (config.defaultPos != null) {
            int idx = config.defaultPos.indexOf(":");
            homeChromosome = idx > 0 ? config.defaultPos.substring(0, idx) : config.defaultPos;
        } else if (this.chromosomeNames != null && this.chromosomeNames.size() > 0) {
            homeChromosome = this.chromosomeNames.get(0);
        } else {
            // TODO -- no place to go
        }


        // Cytobands
        if (config.cytobands != null) {
            cytobandSource = new CytobandMap(config.cytobands);    // Directly supplied, from .genome file
        } else if (config.cytobandBbURL != null) {
            cytobandSource = new CytobandSourceBB(config.cytobandBbURL, this);
        } else if (config.cytobandURL != null) {
            cytobandSource = new CytobandMap(config.cytobandURL);
        }

        // Chromosome aliases
        if (config.aliasURL != null) {
            chromAliasSource = (new ChromAliasFile(config.aliasURL, chromosomeNames));
        } else if (config.chromAliasBbURL != null) {
            chromAliasSource = (new ChromAliasBB(config.chromAliasBbURL, this));
        } else {
            chromAliasSource = (new ChromAliasDefaults(id, chromosomeNames));
        }
        if (config.chromAliases != null) {
            addChrAliases(config.chromAliases);
        }


        addTracks(config);

    }

    private static List<Chromosome> sortChromosomeList(List<Chromosome> chromosomeList, String[] orderedNames) {
        if (orderedNames == null) {
            Collections.sort(chromosomeList, new ChromosomeComparator());
            return chromosomeList;
        } else {
            // Order chromosomes in orderedNames, leave others in original order
            Map<String, Chromosome> chrMap = new HashMap<>();
            for (Chromosome c : chromosomeList) {
                chrMap.put(c.getName(), c);
            }
            List<Chromosome> orderedChromosomes = new ArrayList<>(chromosomeList.size());
            Set<String> orderedNameSet = new HashSet<>(Arrays.asList(orderedNames));
            for (String nm : orderedNames) {
                orderedChromosomes.add(chrMap.get(nm));
            }
            for (Chromosome c : chromosomeList) {
                if (!orderedNameSet.contains(c.getName())) {
                    orderedChromosomes.add(c);
                }
            }
            return orderedChromosomes;
        }
    }

    private void addTracks(GenomeConfig config) {
        // Tracks and hidden tracks
        ArrayList<ResourceLocator> tracks = new ArrayList<>();
        ArrayList<ResourceLocator> hiddenTracks = new ArrayList<>();

        List<TrackConfig> trackConfigs = config.tracks;
        if (trackConfigs == null) {
            trackConfigs = config.annotations;
        }

        if (trackConfigs != null) {

            trackConfigs.forEach((TrackConfig trackConfig) -> {

                String trackPath = trackConfig.url;
                ResourceLocator res = new ResourceLocator(trackPath);
                res.setName(trackConfig.name);
                res.setIndexPath(trackConfig.indexURL);
                res.setFormat(trackConfig.format);
                Integer vw = trackConfig.visibilityWindow;
                if (vw != null) {
                    res.setVisibilityWindow(vw);
                }

                if (trackConfig.searchTrix != null) {
                    res.setTrixURL(trackConfig.searchTrix);
                }

                res.setFeatureInfoURL(trackConfig.infoURL);
                Boolean indexed = trackConfig.indexed;
                if (indexed != null) {
                    res.setIndexed(indexed);
                }

                // Track properties
                TrackProperties properties = new TrackProperties();
                String color = trackConfig.color;
                if (color != null) {
                    try {
                        properties.setColor(ColorUtilities.stringToColor(color.toString()));
                    } catch (Exception e) {
                        log.error("Error parsing color string: " + color, e);
                    }
                }
                String altColor = trackConfig.altColor;
                if (altColor != null) {
                    try {
                        properties.setAltColor(ColorUtilities.stringToColor(altColor.toString()));
                    } catch (Exception e) {
                        log.error("Error parsing color string: " + altColor, e);
                    }
                }
                String displayMode = trackConfig.displayMode;
                if (displayMode != null) {
                    try {
                        Track.DisplayMode dp = Track.DisplayMode.valueOf(stripQuotes(displayMode.toString()));
                        properties.setDisplayMode(dp);
                    } catch (Exception e) {
                        log.error("Error parsing displayMode " + displayMode, e);
                    }
                }
                Integer vizwindow = trackConfig.visibilityWindow;
                if (vizwindow != null) {
                    properties.setFeatureVisibilityWindow(vizwindow);
                } else {
                    // If not explicitly set, assume whole chromosome viz window for annotations
                    properties.setFeatureVisibilityWindow(-1);
                }

                if (trackConfig.min != null) {
                    properties.setMinValue(trackConfig.min);
                }
                if (trackConfig.max != null) {
                    properties.setMaxValue(trackConfig.max);
                }
                res.setTrackProperties(properties);


                Boolean hidden = trackConfig.hidden;

                if (hidden != null && hidden) {
                    hiddenTracks.add(res);
                } else {
                    tracks.add(res);
                }

            });
        }

        setAnnotationResources(tracks);

        if (hiddenTracks.size() > 0) {
            addToFeatureDB(hiddenTracks, this);
        }
    }


    /**
     * Alternate constructor for defining a minimal genome, usually from parsing a chrom.sizes file.  Used to
     * create mock genomes for testing.
     *
     * @param id
     * @param chromosomes
     */
    public Genome(String id, List<Chromosome> chromosomes) {
        this.id = id;
        this.displayName = id;
        this.chrAliasCache = new HashMap<>();
        this.sequence = null;

        chromosomeNames = new ArrayList<>(chromosomes.size());
        chromosomeMap = new LinkedHashMap<>(chromosomes.size());
        for (Chromosome chromosome : chromosomes) {
            chromosomeNames.add(chromosome.getName());
            chromosomeMap.put(chromosome.getName(), chromosome);
        }
    }

    /**
     * Return the canonical chromosome name for the (possibly) alias
     *
     * @param str chromosome or alias name
     * @return the canonical chromsoome name -if the chromosome exists.
     */
    public String getCanonicalChrName(String str) {
        if (str == null) {
            return str;
        } else if (chrAliasCache.containsKey(str)) {
            return chrAliasCache.get(str);
        } else if (chromAliasSource != null) {
            try {
                ChromAlias aliasRecord = chromAliasSource.search(str);
                if (aliasRecord != null) {
                    String chr = aliasRecord.getChr();
                    chrAliasCache.put(str, chr);
                    return chr;
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return str;
    }

    public String getChromosomeDisplayName(String chr) {
        if (this.nameSet != null && this.chromAliasSource != null) {
            String nm = this.chromAliasSource.getChromosomeAlias(chr, this.nameSet);
            return nm != null ? nm : chr;
        } else {
            return chr;
        }
    }

    @Deprecated
    public boolean isKnownChr(String str) {
        return chrAliasCache.containsKey(str);
    }

    /**
     * Add user-defined chromosome aliases directly to the cache.  The input is a collection of chromosome synonym lists.
     * The directionality is determined by the "true" chromosome names, or if chromosome names are not
     * defined by the first entry in the line.  This is a bit complex since, apparently, we allowed synonym lists
     * that did not include the canonical genome.
     *
     * @param synonymsList
     */
    public void addChrAliases(List<List<String>> synonymsList) {

        if (synonymsList == null) return;

        // Convert names to a set for fast "contains" testing.
        Set<String> chrNameSet = new HashSet<>(chromosomeNames);

        for (List<String> synonyms : synonymsList) {

            // See if there is an existing chrom alias record
            ChromAlias chromAlias = null;
            if (chromAliasSource != null) {
                for (String syn : synonyms) {
                    try {
                        chromAlias = chromAliasSource.search(syn);
                    } catch (IOException e) {
                        log.error("Error searching chromosome alias", e);
                    }
                    if (chromAlias != null) break;
                }
            }

            String chr = null;
            if (chromAlias == null) {
                // Find the chromosome name as used in this genome
                for (String syn : synonyms) {
                    if (chrNameSet.contains(syn)) {
                        chr = syn;
                        break;
                    }
                }
                if (chr == null) {
                    chr = synonyms.get(0);
                }
                chromAlias = new ChromAlias(chr);
            } else {
                chr = chromAlias.getChr();
            }

            // Append the new synonyms to the chromAlias
            int idx = chromAlias.values().size() + 1;
            for (String syn : synonyms) {
                chrAliasCache.put(syn, chr);
                chromAlias.put(String.valueOf(idx++), syn);
            }
            if (chromAliasSource != null) {
                chromAliasSource.add(chromAlias);
            }
        }
    }


    /**
     * Extract the user friendly name from an NCBI accession
     * example: gi|125745044|ref|NC_002229.3|  =>  NC_002229.3
     */
    public static String getNCBIName(String name) {

        String[] tokens = name.split("\\|");
        return tokens[tokens.length - 1];
    }


    /**
     * Return the chromosome name associated with the "home" button,  usually the whole genome chromosome.
     *
     * @return
     */
    public String getHomeChromosome() {
        if (homeChromosome == null) {
            if (showWholeGenomeView == false || chromosomeNames.size() == 1 || getLongChromosomeNames().size() > MAX_WHOLE_GENOME_LONG) {
                homeChromosome = chromosomeNames.get(0);
            } else {
                homeChromosome = Globals.CHR_ALL;
            }

        }
        return homeChromosome;
    }

    public String getDefaultPos() {
        return defaultPos == null ? homeChromosome : defaultPos;
    }


    public Chromosome getChromosome(String name) {
        String chrName = getCanonicalChrName(name);
        if (chromosomeMap.containsKey(chrName)) {
            return chromosomeMap.get(chrName);
        } else {
            int length = this.sequence.getChromosomeLength(chrName);
            if (length > 0) {
                int idx = this.chromosomeMap.size();
                Chromosome chromosome = new Chromosome(idx, chrName, length);
                chromosomeMap.put(chrName, chromosome);
                return chromosome;
            } else {
                return null;
            }
        }
    }


    public List<String> getAllChromosomeNames() {
        return chromosomeNames;
    }


    public Collection<Chromosome> getChromosomes() {
        return chromosomeMap.values();
    }


    public long getTotalLength() {
        if (totalLength < 0) {
            totalLength = 0;
            for (Chromosome chr : chromosomeMap.values()) {
                totalLength += chr.getLength();
            }
        }
        return totalLength;
    }


    public long getCumulativeOffset(String chr) {

        Long cumOffset = cumulativeOffsets.get(chr);
        if (cumOffset == null) {
            long offset = 0;
            for (String c : getLongChromosomeNames()) {
                if (chr.equals(c)) {
                    break;
                }
                offset += getChromosome(c).getLength();
            }
            cumOffset = offset;
            cumulativeOffsets.put(chr, cumOffset);
        }
        return cumOffset;
    }

    /**
     * Covert the chromosome coordinate in basepairs to genome coordinates in kilo-basepairs
     *
     * @param chr
     * @param locationBP
     * @return The overall genome coordinate, in kilo-bp
     */
    public int getGenomeCoordinate(String chr, int locationBP) {
        return (int) ((getCumulativeOffset(chr) + locationBP) / 1000);
    }

    /**
     * Translate a genome coordinate, in kilo-basepairs, to a chromosome & position in basepairs.
     *
     * @param genomeKBP The "genome coordinate" in kilo-basepairs.  This is the distance in kbp from the start of the
     *                  first chromosome.
     * @return the position on the corresponding chromosome
     */
    public ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP) {

        long cumOffset = 0;
        List<String> wgChrNames = getLongChromosomeNames();
        for (String c : wgChrNames) {
            int chrLen = getChromosome(c).getLength();
            if ((cumOffset + chrLen) / 1000 > genomeKBP) {
                int bp = (int) (genomeKBP * 1000 - cumOffset);
                return new ChromosomeCoordinate(c, bp);
            }
            cumOffset += chrLen;
        }


        String c = wgChrNames.get(wgChrNames.size() - 1);
        int bp = (int) (genomeKBP - cumOffset) * 1000;
        return new ChromosomeCoordinate(c, bp);
    }


    public String getId() {
        return id;
    }

    public String getBlatDB() {
        return blatDB;
    }


    public void setUcscID(String ucscID) {
        this.ucscID = ucscID;
    }

    public String getUCSCId() {
        return ucscID == null ? id : ucscID;
    }

    public String getSpecies() {
        if (species == null) {
            species = Genome.getSpeciesForID(getUCSCId());
        }
        return species;
    }

    public String getNextChrName(String chr) {
        List<String> chrList = getAllChromosomeNames();
        for (int i = 0; i < chrList.size() - 1; i++) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i + 1);
            }
        }
        return null;
    }

    public String getPrevChrName(String chr) {
        List<String> chrList = getAllChromosomeNames();
        for (int i = chrList.size() - 1; i > 0; i--) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i - 1);
            }
        }
        return null;
    }

    /**
     * Return the nucleotide sequence on the + strand for the genomic interval.  This method can return null
     * if sequence is not available.
     *
     * @param chr
     * @param start start position in "zero-based" coordinates
     * @param end   end position
     * @return sequence, or null if not available
     * @api
     */

    public byte[] getSequence(String chr, int start, int end) {
        return getSequence(getCanonicalChrName(chr), start, end, true);
    }

    public byte[] getSequence(String chr, int start, int end, boolean useCache) {


        if (sequence == null) {
            return null;
        }

        Chromosome c = getChromosome(chr);
        if (c == null) {
            return null;
        }
        end = Math.min(end, c.getLength());
        if (end <= start) {
            return null;
        }
        return sequence.getSequence(chr, start, end);
    }

    public String getDisplayName() {
        return displayName;
    }

    /**
     * Return the reference base at the given position.  Can return null if reference sequence is unknown
     *
     * @param chr
     * @param pos
     * @return the reference base, or null if unknown
     */
    public byte getReference(String chr, int pos) {
        return sequence == null ? null : sequence.getBase(chr, pos);
    }


    public void setCytobands(LinkedHashMap<String, List<Cytoband>> chrCytoMap) {
        this.cytobandSource = new CytobandMap(chrCytoMap);
    }

    public void setGeneTrack(FeatureTrack geneFeatureTrack) {
        this.geneTrack = geneFeatureTrack;
    }

    /**
     * Return the annotation track associated with this genome.  This is a legacy ".genome" artifact, can return null.
     *
     * @return a FeatureTrack, or null
     */
    public FeatureTrack getGeneTrack() {
        return geneTrack;
    }

    /**
     * Return names of "long" chromosomes relative to a fraction of the median length.  The intent here is to
     * remove small contigs in an otherwise well assembled genome to support a "whole genome" view.  We first sort
     * chromosomes in length order, then look for the first large break in size.
     *
     * @return
     */
    public List<String> getLongChromosomeNames() {
        return longChromosomeNames;
    }

    public long getWGLength() {

        if (nominalLength < 0) {
            nominalLength = 0;
            for (String chrName : getLongChromosomeNames()) {
                Chromosome chr = getChromosome(chrName);
                nominalLength += chr.getLength();
            }
        }
        return nominalLength;
    }


    // TODO A hack (obviously),  we need to record a species in the genome definitions to support old style
    // blat servers.
    private static Map<String, String> ucscSpeciesMap;

    private static synchronized String getSpeciesForID(String id) {
        if (ucscSpeciesMap == null) {
            ucscSpeciesMap = new HashMap<>();

            InputStream is = null;

            try {
                is = Genome.class.getResourceAsStream("speciesMapping.txt");
                BufferedReader br = new BufferedReader(new InputStreamReader(is));

                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    if (nextLine.startsWith("#")) continue;
                    String[] tokens = Globals.tabPattern.split(nextLine);
                    if (tokens.length == 2) {
                        ucscSpeciesMap.put(tokens[0].trim(), tokens[1].trim());
                    } else {
                        log.error("Unexpected number of tokens in species mapping file for line: " + nextLine);
                    }
                }
            } catch (IOException e) {
                log.error("Error reading species mapping table", e);
            } finally {
                if (is != null) try {
                    is.close();
                } catch (IOException e) {
                    log.error("", e);
                }
            }

        }

        for (Map.Entry<String, String> entry : ucscSpeciesMap.entrySet()) {
            if (id.startsWith(entry.getKey())) {
                return entry.getValue();
            }
        }
        return null;
    }

    // Map some common IGV genome IDs to UCSC equivalents.  Primarily for BLAT usage
    private static Map<String, String> ucsdIDMap;

    static {
        ucsdIDMap = new HashMap<>();
        ucsdIDMap.put("1kg_ref", "hg18");
        ucsdIDMap.put("1kg_v37", "hg19");
        ucsdIDMap.put("b37", "hg19");
    }

    public void setAnnotationResources(List<ResourceLocator> annotationResources) {
        this.annotationResources = annotationResources;
    }

    public List<ResourceLocator> getAnnotationResources() {
        return annotationResources;
    }

    public Map<ResourceLocator, List<Track>> getAnnotationTracks() {
        return annotationTracks;
    }

    public void setAnnotationTracks(Map<ResourceLocator, List<Track>> annotationTracks) {
        this.annotationTracks = annotationTracks;
    }

    /**
     * Mock genome for unit tests
     */

    private Genome(String id) {
        this.id = id;
    }

    ;

    public static Genome mockGenome = new Genome("hg19");

    public void setShowWholeGenomeView(boolean showWholeGenomeView) {
        this.showWholeGenomeView = showWholeGenomeView;
    }

    public boolean getShowWholeGenomeView() {
        return showWholeGenomeView;
    }

    public void setLongChromosomeNames(List<String> chrNames) {
        this.longChromosomeNames = chrNames;
    }

    public Map<String, Liftover> getLiftoverMap() {
        return liftoverMap;
    }

    public void setLiftoverMap(Map<String, Liftover> liftoverMap) {
        this.liftoverMap = liftoverMap;
    }

    public void setChromAliasSource(ChromAliasSource chromAliasSource) {
        this.chromAliasSource = chromAliasSource;
    }

    public ChromAlias getAliasRecord(String chr) throws IOException {
        return chromAliasSource == null ? null : chromAliasSource.search(chr);
    }

    public List<Cytoband> getCytobands(String chrName) {
        if (cytobandSource != null) {
            try {
                return cytobandSource.getCytobands(chrName);
            } catch (IOException e) {
                log.error("Error fetching cytobands for chr: " + chrName, e);
                return null;
            }
        } else {
            return null;
        }
    }

    private static String stripQuotes(String str) {
        if (str.startsWith("\"")) {
            return str.substring(1, str.length() - 1);  // Assume also ends with
        } else {
            return str;
        }
    }

    private static void addToFeatureDB(List<ResourceLocator> locators, Genome genome) {
        for (ResourceLocator locator : locators) {
            try {
                FeatureReader featureReader = TribbleFeatureSource.getBasicReader(locator, genome);
                CloseableTribbleIterator<Feature> iter = featureReader.iterator();
                while (iter.hasNext()) {
                    Feature f = iter.next();
                    if (f instanceof NamedFeature) {
                        FeatureDB.addFeature((NamedFeature) f, genome);
                    }
                }
            } catch (IOException e) {
                log.error("Error loading " + locator.getPath());
            }
        }
    }

    public List<String> computeLongChromosomeNames() {

        List<String> longChromosomeNames = new ArrayList<>();
        if (chromosomeMap.size() < 100) {
            // Keep all chromosomes > 10% of the mean in length
            double[] lengths = new double[chromosomeMap.size()];
            int idx = 0;
            for (Chromosome c : chromosomeMap.values()) {
                lengths[idx++] = c.getLength();
            }
            double mean = StatUtils.mean(lengths);
            double min = 0.1 * mean;
            for (String chr : getAllChromosomeNames()) {
                if (chromosomeMap.get(chr).getLength() > min) {
                    longChromosomeNames.add(chr);
                }
            }
        } else {
            // Long list, likely many small contigs.  Search for a break between long (presumably assembled)
            // chromosomes and small contigs.
            List<Chromosome> allChromosomes = new ArrayList<>(chromosomeMap.values());
            allChromosomes.sort((c1, c2) -> c2.getLength() - c1.getLength());

            Chromosome lastChromosome = null;
            Set<String> tmp = new HashSet<>();
            for (Chromosome c : allChromosomes) {
                if (lastChromosome != null) {
                    double delta = lastChromosome.getLength() - c.getLength();
                    if (delta / lastChromosome.getLength() > 0.7) {
                        break;
                    }
                }
                tmp.add(c.getName());
                lastChromosome = c;
            }


            for (String chr : getAllChromosomeNames()) {
                if (tmp.contains(chr)) {
                    longChromosomeNames.add(chr);
                }
            }
        }

        return longChromosomeNames;

    }
}
