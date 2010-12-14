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

package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.renderer.GeneTrackRenderer;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.*;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.net.URLDecoder;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Oct 21, 2009
 * Time: 10:14:24 PM
 * To change this template use File | Settings | File Templates.
 */


public class GFFParser implements FeatureParser {

    static Logger log = Logger.getLogger(GFFParser.class);

    static HashSet exonTerms = new HashSet();
    static HashSet utrTerms = new HashSet();
    static HashSet<String> geneParts = new HashSet();
    static HashSet<String> ignoredTypes = new HashSet();

    static String[] nameFields = {"gene", "Name", "name", "primary_name", "Locus", "locus", "alias", "systematic_id", "ID"};


    static String[] tokens = new String[200];

    static {
        utrTerms.add("five_prime_UTR");
        utrTerms.add("three_prime_UTR");
        utrTerms.add("5'-utr");
        utrTerms.add("3'-utr");
        utrTerms.add("3'-UTR");
        utrTerms.add("5'-UTR");
        utrTerms.add("5utr");
        utrTerms.add("3utr");
    }

    static {
        exonTerms.addAll(utrTerms);
        exonTerms.add("exon");
        exonTerms.add("coding_exon");
        exonTerms.add("CDS");
        exonTerms.add("cds");

    }


    static {
        geneParts.addAll(exonTerms);
        geneParts.add("transcript");
        geneParts.add("processed_transcript");
        geneParts.add("mrna");
        geneParts.add("mRNA");
        geneParts.add("promoter");
        geneParts.add("intron");
        geneParts.add("CDS_parts");

    }

    static {
        ignoredTypes.add("start_codon");
        ignoredTypes.add("stop_codon");
        ignoredTypes.add("Contig");
        ignoredTypes.add("RealContig");
        ignoredTypes.add("intron");
    }

    /*
    static {
        ignoreAttributes.add("ID");
        ignoreAttributes.add("Parent");
        ignoreAttributes.add("wormprep");
        ignoreAttributes.add("Index");
        ignoreAttributes.add("color");
        ignoreAttributes.add("colour");
        ignoreAttributes.add("exonNumber");

    }
    */


    Helper helper;

    Map<String, GFF3Transcript> transcriptCache = new HashMap(50000);
    Map<String, BasicFeature> geneCache = new HashMap(50000);
    private TrackProperties trackProperties = null;

    static StringBuffer buf = new StringBuffer(1000);


    public GFFParser(String path) {
        // Assume V2 until proven otherwise
        if (path.toLowerCase().endsWith("gff3") || path.toLowerCase().endsWith("gff3.gz")) {
            helper = new GFF3Helper();
        } else {
            helper = new GFF2Helper();
        }
    }

    /**
     * By definition this is a feature file
     */
    public boolean isFeatureFile(ResourceLocator locator) {
        return true;
    }


    public List<FeatureTrack> loadTracks(ResourceLocator locator) {

        if (locator.getPath().toLowerCase().endsWith("gff3") || locator.getPath().toLowerCase().endsWith("gff3.gz")) {
            helper = new GFF3Helper();
        }

        AsciiLineReader reader = null;
        try {
            reader = ParsingUtils.openAsciiReader(locator);

            List<org.broad.tribble.Feature> features = loadFeatures(reader);

            FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features));
            track.setName(locator.getTrackName());
            track.setRendererClass(IGVFeatureRenderer.class);
            track.setMinimumHeight(35);
            track.setHeight(45);
            track.setRendererClass(GeneTrackRenderer.class);

            if (trackProperties != null) {
                track.setTrackProperties(trackProperties);
            }

            List<FeatureTrack> tracks = new ArrayList();
            tracks.add(track);
            return tracks;

        } catch (Exception ex) {
            log.error(ex);
            throw new RuntimeException(ex);

        } finally {
            if (reader != null) {
                reader.close();

            }
        }
    }

    /**
     * Method description
     *
     * @param reader
     * @return
     */


    public List<org.broad.tribble.Feature> loadFeatures(AsciiLineReader reader) {
        List<org.broad.tribble.Feature> features = new ArrayList();
        URLDecoder d = new URLDecoder();
        String line = null;
        try {

            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            Set<String> featuresToHide = new HashSet();
            while ((line = reader.readLine()) != null) {

                if (line.startsWith("##gff-version") && line.endsWith("3")) {
                    helper = new GFF3Helper();
                }


                if (line.startsWith("#")) {
                    if (line.startsWith("#track")) {
                        trackProperties = new TrackProperties();
                        ParsingUtils.parseTrackLine(line, trackProperties);
                    } else if (line.startsWith("#nodecode")) {
                        helper.setUrlDecoding(false);
                    } else if (line.startsWith("#hide")) {
                        String[] kv = line.split("=");
                        if (kv.length > 1) {
                            featuresToHide.addAll(Arrays.asList(kv[1].split(",")));
                        }
                    }
                    continue;
                }


                int nTokens = ParsingUtils.split(line, tokens, '\t');

                // GFF files have 9 tokens
                if (nTokens < 9) {
                    continue;
                }

                // The type
                String featureType = new String(tokens[2].trim());

                if (ignoredTypes.contains(featureType)) {
                    continue;
                }


                String chromosome = genome.getChromosomeAlias(tokens[0]);

                // GFF coordinates are 1-based inclusive (length = end - start + 1)
                // IGV (UCSC) coordinates are 0-based exclusive.  Adjust start and end accordingly
                int start;
                int end;
                try {
                    start = Integer.parseInt(tokens[3]) - 1;
                }
                catch (NumberFormatException ne) {
                    throw new ParserException("Column 4 must contain a numeric value", reader.getCurrentLineNumber(), line);
                }

                try {
                    end = Integer.parseInt(tokens[4]);
                }
                catch (NumberFormatException ne) {
                    throw new ParserException("Column 5 must contain a numeric value", reader.getCurrentLineNumber(), line);
                }

                Strand strand = convertStrand(tokens[6]);

                String attributeString = tokens[8];

                LinkedHashMap<String, String> attributes = new LinkedHashMap();
                //attributes.put("Type", featureType);
                helper.parseAttributes(attributeString, attributes);

                String description = getDescription(attributes, featureType);

                String id = helper.getID(attributes);

                String[] parentIds = helper.getParentIds(attributes, attributeString);

                if (featureType.equals("CDS_parts")) {
                    for (String pid : parentIds) {
                        getGFF3Transcript(pid).addCDSParts(chromosome, start, end, description);
                    }

                } else if (featureType.equals("intron")) {

                    for (String pid : parentIds) {
                        getGFF3Transcript(pid).addCDSParts(chromosome, start, end, description);
                    }
                } else if (exonTerms.contains(featureType) && parentIds != null && parentIds.length > 0) {

                    String name = getName(attributes);
                    int phase = -1;
                    String phaseString = tokens[7].trim();
                    if (!phaseString.equals(".")) {
                        try {
                            phase = Integer.parseInt(phaseString);

                        } catch (NumberFormatException numberFormatException) {

                            // Just skip setting the phase
                            log.error("GFF3 error: non numeric phase: " + phaseString);
                        }
                    }

                    // Make a copy of the exon record for each parent
                    for (String pid : parentIds) {

                        Exon exon = new Exon(chromosome, start, end, strand);
                        exon.setDescription(description);
                        exon.setUTR(utrTerms.contains(featureType));

                        if (phase >= 0) {
                            exon.setPhase(phase);

                        }
                        exon.setName(name);

                        if (featureType.equals("exon")) {
                            getGFF3Transcript(pid).addExon(exon);
                        } else if (featureType.equals("CDS")) {
                            getGFF3Transcript(pid).addCDS(exon);
                        } else if (featureType.equals("five_prime_UTR") || featureType.equals("5'-UTR")) {
                            getGFF3Transcript(pid).setFivePrimeUTR(exon);
                        } else if (featureType.equals("three_prime_UTR") || featureType.equals("3'-UTR")) {
                            getGFF3Transcript(pid).setThreePrimeUTR(exon);
                        }
                    }


                } else {
                    BasicFeature f = new BasicFeature(chromosome, start, end, strand);
                    f.setName(getName(attributes));
                    f.setDescription(description);

                    if (attributes.containsKey("color")) {
                        f.setColor(ColorUtilities.getColorFromString(attributes.get("color")));
                    }
                    if (attributes.containsKey("Color")) {
                        f.setColor(ColorUtilities.getColorFromString(attributes.get("Color")));
                    }


                    if (id == null) {
                        features.add(f);
                    } else {
                        f.setIdentifier(id);

                        if (featureType.equals("gene")) {
                            geneCache.put(id, f);
                            if (featuresToHide.contains(featureType)) {
                                FeatureDB.addFeature(f);
                            } else {
                                features.add(f);
                            }
                        } else if (featureType.equals("mRNA") || featureType.equals("transcript")) {
                            String pid = null;
                            if (parentIds != null && parentIds.length > 0) {
                                pid = parentIds[0];
                            }
                            getGFF3Transcript(id).transcript(f, pid);
                        } else if (featuresToHide.contains(featureType)) {
                            FeatureDB.addFeature(f);
                        } else {
                            features.add(f);
                        }
                    }
                }
            }

            // Create and add IGV genes
            for (GFF3Transcript transcript : transcriptCache.values()) {
                BasicFeature igvGene = transcript.createTranscript();
                if (igvGene != null) {
                    features.add(igvGene);
                }
            }

            FeatureDB.addFeatures(features);

        }
        catch (ParserException
                e) {
            throw e;
        }
        catch (Exception
                ex) {
            log.error("Error parsing GFF file", ex);
            if (line != null && reader.getCurrentLineNumber() != 0) {

                throw new ParserException(ex.getMessage(), ex, reader.getCurrentLineNumber(), line);
            } else {
                throw new RuntimeException(ex);
            }
        }

        return features;
    }


    static String getDescription(Map<String, String> attributes, String type) {
        buf.setLength(0);
        buf.append(type);
        buf.append("<br>");
        for (Map.Entry<String, String> att : attributes.entrySet()) {

            //if (!ignoreAttributes.contains(att.getKey())) {
            String attValue = att.getValue().replaceAll(";", "<br>");
            buf.append(att.getKey());
            buf.append(" = ");
            buf.append(attValue);
            buf.append("<br>");
            //}
        }
        String description = buf.toString();

        return description;
    }

    private Strand convertStrand(String strandString) {
        Strand strand = Strand.NONE;
        if (strandString.equals("-")) {
            strand = Strand.NEGATIVE;
        } else if (strandString.equals("+")) {
            strand = Strand.POSITIVE;
        }

        return strand;
    }

    private GFF3Transcript getGFF3Transcript(String id) {
        GFF3Transcript transcript = transcriptCache.get(id);
        if (transcript == null) {
            transcript = new GFF3Transcript(id);
            transcriptCache.put(id, transcript);
        }
        return transcript;
    }


    String getName(Map<String, String> attributes) {

        if (attributes.size() == 0) {
            return null;
        }

        for (String nf : nameFields) {
            if (attributes.containsKey(nf)) {
                return attributes.get(nf);
            }
        }

        return attributes.values().iterator().next();

    }


    public static void splitFileByType(String gffFile, String outputDirectory) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(gffFile));
        String nextLine;

        Map<String, PrintWriter> writers = new HashMap();
        while ((nextLine = br.readLine()) != null) {
            nextLine = nextLine.trim();
            if (!nextLine.startsWith("#")) {
                int nTokens = ParsingUtils.split(nextLine.trim().replaceAll("\"", ""), tokens, '\t');

                // GFF files have 9 columns
                String type = tokens[2];
                if (geneParts.contains(type)) {
                    type = "gene";
                }
                if (!writers.containsKey(type)) {
                    writers.put(type,
                            new PrintWriter(new FileWriter(new File(outputDirectory,
                                    type + ".gff"))));
                }
            }
        }
        br.close();

        br = new BufferedReader(new FileReader(gffFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            nextLine = nextLine.trim();
            if (nextLine.startsWith("#")) {
                ParsingUtils.split(nextLine.trim().replaceAll("\"", ""), tokens, '\t');

                // GFF files have 9 columns
                for (PrintWriter pw : writers.values()) {
                    pw.println(nextLine);
                }
            } else {
                int nTokens = ParsingUtils.split(nextLine.trim().replaceAll("\"", ""), tokens, '\t');
                String type = tokens[2];
                if (geneParts.contains(type)) {
                    type = "gene";
                }
                currentWriter = writers.get(type);

                if (currentWriter != null) {
                    currentWriter.println(nextLine);
                } else {
                    System.out.println("No writer for: " + type);
                }
            }

        }

        br.close();
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    class GFF3Transcript {

        private String id;
        private Set<Exon> exons = new HashSet();
        private List<Exon> cdss = new ArrayList();
        private List<Exon> cdsParts = null;
        private Exon fivePrimeUTR;
        private Exon threePrimeUTR;
        private BasicFeature transcript;
        private String parentId;
        String chr = null;
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        String description;

        GFF3Transcript(String id) {
            this.id = id;
        }

        void transcript(BasicFeature mRNA, String parent) {
            this.transcript = mRNA;
            this.parentId = parent;
            if (mRNA.getName() == null) {
                mRNA.setName(mRNA.getIdentifier());
            }
            int prefixIndex = mRNA.getName().indexOf(":");
            if (prefixIndex > 0) {
                mRNA.setName(mRNA.getName().substring(prefixIndex + 1));
            }

        }

        void setFivePrimeUTR(Exon exon) {
            fivePrimeUTR = exon;
            this.start = Math.min(exon.getStart(), start);
            this.end = Math.max(exon.getEnd(), end);
        }

        void setThreePrimeUTR(Exon exon) {
            threePrimeUTR = exon;
            this.start = Math.min(exon.getStart(), start);
            this.end = Math.max(exon.getEnd(), end);
        }

        void addExon(Exon exon) {
            exons.add(exon);
            this.start = Math.min(exon.getStart(), start);
            this.end = Math.max(exon.getEnd(), end);
        }

        void addCDS(Exon cds) {
            cdss.add(cds);
        }

        void addCDSParts(String chr, int start, int end, String desc) {
            this.chr = chr;
            this.start = Math.min(this.start, start);
            this.end = Math.max(this.end, end);
            appendDescription(desc);
        }

        void appendDescription(String desc) {
            if (description == null) {
                description = "<html>";
            } else {
                description += "<br>---------<br>";
            }
            description += desc;
        }

        /**
         * Create a transcript from its constituitive parts. "
         *
         * @return
         */
        BasicFeature createTranscript() {

            Strand strand = Strand.NONE;
            String name = null;

            // Combine exon & cds  if feature contains both CDS and CDS_parts use CDS_parts
            List<Exon> cdsList = cdsParts;
            if (cdsList == null) {
                cdsList = cdss;
            }
            if (exons.size() > 0) {
                while (!cdsList.isEmpty()) {
                    Exon cds = cdsList.get(0);
                    Exon exon = findMatchingExon(cds);
                    if (exon == null) {
                        exons.add(cds);
                    } else {
                        exon.setCodingStart(cds.getStart());
                        exon.setCodingEnd(cds.getEnd());
                        exon.setReadingFrame(cds.getReadingShift());
                    }
                    cdsList.remove(0);
                }
                for (Exon exon : exons) {
                    chr = exon.getChr();
                    strand = exon.getStrand();
                    start = Math.min(exon.getStart(), start);
                    end = Math.max(exon.getEnd(), end);
                    name = exon.getName();
                }
            }

            if (transcript == null) {
                // transcript is implied
                transcript = new BasicFeature(chr, start, end, strand);
                transcript.setIdentifier(id);
                transcript.setName(name == null ? id : name);
                transcript.setDescription(description);
            }

            if ((parentId != null) && geneCache.containsKey(parentId)) {
                BasicFeature gene = geneCache.get(parentId);
                geneCache.remove(parentId);
                if (transcript.getName() == null && gene.getName() != null) {
                    transcript.setName(gene.getName());
                }
                transcript.setDescription("Transcript<br>" + transcript.getDescription() + "<br>--------<br>Gene<br>" + gene.getDescription());

                // mRNA.setName(gene.getName());
            }

            for (Exon exon : exons) {
                transcript.addExon(exon);
            }


            transcript.sortExons();

            // If 5'UTR is represented by an exon, adjust its start, else add an exon to represnet 5'utr
            if (fivePrimeUTR != null) {
                fivePrimeUTR.setUTR(true);
                transcript.addExon(fivePrimeUTR);
                Exon exon = findMatchingExon(fivePrimeUTR);
                if (exon != null) {
                    if (exon.getStrand() == Strand.POSITIVE) {
                        exon.setStart(fivePrimeUTR.getEnd());
                    } else {
                        exon.setEnd(fivePrimeUTR.getStart());
                    }
                }
            }

            if (threePrimeUTR != null) {
                threePrimeUTR.setUTR(true);
                transcript.addExon(threePrimeUTR);
                Exon exon = findMatchingExon(threePrimeUTR);
                if (exon != null) {
                    if (exon.getStrand() == Strand.POSITIVE) {
                        exon.setEnd(threePrimeUTR.getStart());
                    } else {
                        exon.setStart(threePrimeUTR.getEnd());
                    }
                }
            }

            //transcript.setDescription(getDescription(transcript.getAttributes()));

            return transcript;
        }

        Exon findMatchingExon(IGVFeature cds) {
            for (Exon exon : exons) {
                if (exon.contains(cds)) {
                    return exon;
                }
            }
            return null;
        }
    }

    protected interface Helper {

        String[] getParentIds(Map<String, String> attributes, String attributeString);

        void parseAttributes(String attributeString, Map<String, String> map);

        String getID(Map<String, String> attributes);

        void setUrlDecoding(boolean b);

    }

    class GFF2Helper implements Helper {

        String[] idFields = {"systematic_id", "ID", "transcript_id", "Name", "name", "primary_name", "gene", "Locus", "locus", "alias"};


        public void setUrlDecoding(boolean b) {
            // Ignored,  GFF files are never url DECODED
        }


        public void parseAttributes(String description, Map<String, String> kvalues) {

            List<String> kvPairs = StringUtils.breakQuotedString(description.trim(), ';');

            for (String kv : kvPairs) {
                List<String> tokens = StringUtils.breakQuotedString(kv, ' ');
                if (tokens.size() >= 2) {
                    String key = tokens.get(0).trim().replaceAll("\"", "");
                    String value = tokens.get(1).trim().replaceAll("\"", "");
                    kvalues.put(key, value);
                }
            }
        }


        public String[] getParentIds(Map<String, String> attributes, String attributeString) {

            String[] parentIds = new String[1];
            if (attributes.isEmpty()) {
                parentIds[0] = attributeString;
            } else {
                parentIds[0] = attributes.get("id");
                if (parentIds[0] == null) {
                    parentIds[0] = attributes.get("mRNA");
                }
                if (parentIds[0] == null) {
                    parentIds[0] = attributes.get("systematic_id");
                }
                if (parentIds[0] == null) {
                    parentIds[0] = attributes.get("transcript_id");
                }
                if (parentIds[0] == null) {
                    parentIds[0] = attributes.get("gene");
                }
                if (parentIds[0] == null) {
                    parentIds[0] = attributes.get("transcriptId");
                }
                if (parentIds[0] == null) {
                    parentIds[0] = attributes.get("proteinId");
                }
            }
            return parentIds;

        }


        public String getID(Map<String, String> attributes) {
            for (String nf : idFields) {
                if (attributes.containsKey(nf)) {
                    return attributes.get(nf);
                }
            }
            return getName(attributes);
        }
    }

    class GFF3Helper implements Helper {


        private boolean useUrlDecoding = true;


        public String[] getParentIds(Map<String, String> attributes, String ignored) {
            String parentIdString = attributes.get("Parent");
            if (parentIdString != null) {
                return attributes.get("Parent").split(",");
            } else {
                return null;
            }
        }


        public void parseAttributes(String description, Map<String, String> kvalues) {

            List<String> kvPairs = StringUtils.breakQuotedString(description.trim(), ';');
            for (String kv : kvPairs) {
                //int nValues = ParsingUtils.split(kv, tmp, '=');
                List<String> tmp = StringUtils.breakQuotedString(kv, '=');
                int nValues = tmp.size();
                if (nValues > 0) {
                    String key = tmp.get(0).trim();
                    String value = ((nValues == 1) ? "" : tmp.get(1).trim());

                    if (useUrlDecoding) {
                        key = URLDecoder.decode(key);
                        value = URLDecoder.decode(value);
                    }
                    kvalues.put(key, value);
                } else {
                    log.info("No attributes: " + description);
                }
            }
        }

        public void setUrlDecoding(boolean useUrlDecoding) {
            this.useUrlDecoding = useUrlDecoding;
        }


        public String getID(Map<String, String> attributes) {
            String id = attributes.get("ID");
            return id;
        }
    }


    public static void main(String[] args) throws IOException {
        if (args.length < 2) {
            System.out.println("SpitFilesByType <gffFile> <outputDirectory>");
            return;
        }
        splitFileByType(args[0], args[1]);
    }
}
