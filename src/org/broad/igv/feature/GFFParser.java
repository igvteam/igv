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

package org.broad.igv.feature;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.GeneTrackRenderer;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.CI;
import org.broad.tribble.Feature;

import java.io.*;
import java.util.*;

/**
 *
 * User: jrobinso
 */


public class GFFParser implements FeatureParser {

    static Logger log = Logger.getLogger(GFFParser.class);

    static HashSet exonTerms = new HashSet();
    static HashSet utrTerms = new HashSet();
    static HashSet<String> geneParts = new HashSet();
    static HashSet<String> ignoredTypes = new HashSet();

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


    Helper helper;

    Map<String, GFF3Transcript> transcriptCache = new HashMap(50000);
    Map<String, BasicFeature> geneCache = new HashMap(50000);
    private TrackProperties trackProperties = null;
    Set<String> featuresToHide = new HashSet();

    public static boolean isGFF(String path) {
        String lowpath = path.toLowerCase();
        lowpath = lowpath.replace(".gz", "");
        return lowpath.endsWith("gff3") ||
                lowpath.endsWith("gvf") ||
                lowpath.endsWith("gff") ||
                lowpath.endsWith("gtf");
    }

    private static boolean isGFF3(String path){
    String lowpath = path.toLowerCase().replace(".gz", "");
    return lowpath.endsWith("gff3") || lowpath.endsWith("gvf");
    }


    public GFFParser(String path) {
        // Assume V2 until proven otherwise

        if (isGFF3(path)) {
            helper = new GFF3Helper();
        } else {
            helper = new GFF2Helper();
        }
    }

    public List<FeatureTrack> loadTracks(ResourceLocator locator, Genome genome) {

        String path = locator.getPath().toLowerCase();
        if (isGFF3(path)) {
            helper = new GFF3Helper();
        }

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(locator);

            List<org.broad.tribble.Feature> features = loadFeatures(reader, genome);

            FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features, genome));
            track.setName(locator.getTrackName());
            track.setRendererClass(IGVFeatureRenderer.class);
            track.setMinimumHeight(35);
            track.setHeight(45);
            track.setRendererClass(GeneTrackRenderer.class);

            if (trackProperties != null) {
                track.setProperties(trackProperties);
            }

            List<FeatureTrack> tracks = new ArrayList();
            tracks.add(track);
            return tracks;

        } catch (Exception ex) {
            log.error(ex);
            throw new RuntimeException(ex);

        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }

            }
        }
    }

    private void readHeaderLine(String line){
        if (line.startsWith("#track") || line.startsWith("##track")) {
            trackProperties = new TrackProperties();
            ParsingUtils.parseTrackLine(line, trackProperties);
        } else if (line.startsWith("#nodecode") || line.startsWith("##nodecode")) {
            helper.setUrlDecoding(false);
        } else if (line.startsWith("#hide") || line.startsWith("##hide")) {
            String[] kv = line.split("=");
            if (kv.length > 1) {
                featuresToHide.addAll(Arrays.asList(kv[1].split(",")));
            }
        } else if (line.startsWith("#displayName") || line.startsWith("##displayName")) {
            String[] nameTokens = line.split("=");
            if (nameTokens.length < 2) {
                helper.setNameFields(null);
            } else {
                String[] fields = nameTokens[1].split(",");
                helper.setNameFields(fields);
            }
        }
    }

    private void processExon(String featureType, String[] parentIds, String chromosome,
                             int start, int end, Strand strand, Map<String, String> attributes, String phaseString){
        String name = getName(attributes);
        int phase = -1;
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
            exon.setAttributes(attributes);
            exon.setUTR(utrTerms.contains(featureType));

            if (phase >= 0) {
                exon.setPhase(phase);

            }
            exon.setName(name);

            if (featureType.equalsIgnoreCase("exon")) {
                getGFF3Transcript(pid).addExon(exon);
            } else if (featureType.equals("CDS")) {
                getGFF3Transcript(pid).addCDS(exon);
            } else if (featureType.equals("five_prime_UTR") || featureType.equals("5'-UTR")) {
                getGFF3Transcript(pid).setFivePrimeUTR(exon);
            } else if (featureType.equals("three_prime_UTR") || featureType.equals("3'-UTR")) {
                getGFF3Transcript(pid).setThreePrimeUTR(exon);
            }
        }
    }

    private BasicFeature generateFeature(String featureType, String[] parentIds, String chromosome,
                                         int start, int end, Strand strand, Map<String, String> attributes){


        BasicFeature f = new BasicFeature(chromosome, start, end, strand);
        String name = getName(attributes);
        f.setName(name);
        f.setAttributes(attributes);

        if (attributes.containsKey("color")) {
            f.setColor(ColorUtilities.stringToColor(attributes.get("color")));
        }

        String id = helper.getID(attributes);
        if (id == null) {
            return f;
        } else {
            f.setIdentifier(id);
        }

        if (featuresToHide.contains(featureType)) {
            if (IGV.hasInstance()) FeatureDB.addFeature(f);
            return null;
        }

        if (featureType.equalsIgnoreCase("gene")) {
            geneCache.put(id, f);
        } else if (featureType.equalsIgnoreCase("mRNA") || featureType.equalsIgnoreCase("transcript")) {
            String pid = null;
            if (parentIds != null && parentIds.length > 0) {
                pid = parentIds[0];
            }
            getGFF3Transcript(id).transcript(f, pid);
        }
        return f;
    }


    public List<org.broad.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {
        List<org.broad.tribble.Feature> features = new ArrayList();
        String line = null;
        int lineNumber = 0;
        try {


            while ((line = reader.readLine()) != null) {
                lineNumber++;

                if (line.startsWith("##gff-version") && line.endsWith("3")) {
                    helper = new GFF3Helper();
                }


                if (line.startsWith("#")) {
                    readHeaderLine(line);
                    continue;
                }

                String[] tokens = Globals.tabPattern.split(line, -1);
                int nTokens = tokens.length;

                // GFF files have 9 tokens
                if (nTokens < 9) {
                    // Maybe its a track line?
                    if (line.startsWith("track")) {
                        trackProperties = new TrackProperties();
                        ParsingUtils.parseTrackLine(line, trackProperties);
                        continue;
                    } else {
                        String msg = String.format("GFF line expected to have 9 tokens, but has %d", nTokens);
                        throw new ParserException(msg, lineNumber, line);
                    }
                }

                String featureType = new String(tokens[2].trim());

                if (ignoredTypes.contains(featureType)) {
                    continue;
                }

                String chromosome = genome == null ? tokens[0] : genome.getChromosomeAlias(tokens[0]);

                // GFF coordinates are 1-based inclusive (length = end - start + 1)
                // IGV (UCSC) coordinates are 0-based exclusive.  Adjust start and end accordingly
                int start;
                int end;
                int col = 3;
                try {
                    start = Integer.parseInt(tokens[col]) - 1;
                    col++;
                    end = Integer.parseInt(tokens[col]);
                } catch (NumberFormatException ne) {
                    String msg = String.format("Column %d must contain a numeric value. %s", col + 1, ne.getMessage());
                    throw new ParserException(msg, lineNumber, line);
                }

                Strand strand = convertStrand(tokens[6]);
                String attributeString = tokens[8];

                CI.CILinkedHashMap<String> attributes = new CI.CILinkedHashMap();

                helper.parseAttributes(attributeString, attributes);
                String[] parentIds = helper.getParentIds(attributes, attributeString);


                if (featureType.equals("CDS_parts") || featureType.equals("intron")) {
                    for (String pid : parentIds) {
                        getGFF3Transcript(pid).addCDSParts(chromosome, start, end);
                    }

                } else if (exonTerms.contains(featureType) && parentIds != null && parentIds.length > 0 &&
                        parentIds[0] != null && parentIds[0].length() > 0 && !parentIds[0].equals(".")) {

                    String phaseString = tokens[7].trim();
                    processExon(featureType, parentIds, chromosome, start, end, strand, attributes, phaseString);

                } else {

                    BasicFeature f = generateFeature(featureType, parentIds, chromosome, start, end, strand, attributes);
                    if(f != null){
                        features.add(f);
                    }

                }
            }

            // Create and add IGV genes
            for (GFF3Transcript transcript : transcriptCache.values()) {
                Feature igvTranscript = transcript.createTranscript();
                if (igvTranscript != null) {
                    features.add(igvTranscript);
                }
            }


        } catch (ParserException e) {
            throw e;
        } catch (IOException ex) {
            log.error("Error reading GFF file", ex);
            if (line != null && lineNumber != 0) {
                throw new ParserException(ex.getMessage(), ex, lineNumber, line);
            } else {
                throw new RuntimeException(ex);
            }
        }

        return features;
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

        return helper.getName(attributes);

    }


    public static void splitFileByType(String gffFile, String outputDirectory) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(gffFile));
        String nextLine;
        String ext = "." + gffFile.substring(gffFile.length() - 4);

        Map<String, PrintWriter> writers = new HashMap();

        while ((nextLine = br.readLine()) != null) {
            nextLine = nextLine.trim();
            if (!nextLine.startsWith("#")) {
                String[] tokens = Globals.tabPattern.split(nextLine.trim().replaceAll("\"", ""), -1);

                // GFF files have 9 columns
                String type = tokens[2];
                if (geneParts.contains(type)) {
                    type = "gene";
                }
                if (!writers.containsKey(type)) {
                    writers.put(type,
                            new PrintWriter(new FileWriter(new File(outputDirectory, type + ext))));
                }
            }
        }
        br.close();

        br = new BufferedReader(new FileReader(gffFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            nextLine = nextLine.trim();
            if (nextLine.startsWith("#")) {
                // GFF files have 9 columns
                for (PrintWriter pw : writers.values()) {
                    pw.println(nextLine);
                }
            } else {
                String[] tokens = Globals.tabPattern.split(nextLine.trim().replaceAll("\"", ""), -1);
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
        private Exon fivePrimeUTR;
        private Exon threePrimeUTR;
        private BasicFeature transcript;
        private String parentId;
        String chr = null;
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        String description;
        Map<String, String> attributes;

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
            this.start = Math.min(cds.getStart(), start);
            this.end = Math.max(cds.getEnd(), end);
        }

        void addCDSParts(String chr, int start, int end) {
            this.chr = chr;
            this.start = Math.min(this.start, start);
            this.end = Math.max(this.end, end);
        }

        /**
         * Create a transcript from its constituitive parts. "
         *
         * @return
         */
        Feature createTranscript() {

            Strand strand = Strand.NONE;
            String name = null;

            // Combine CDS and exons
            while (!cdss.isEmpty()) {
                Exon cds = cdss.get(0);
                Exon exon = findMatchingExon(cds);
                if (exon == null) {
                    cds.setCodingStart(cds.getStart());
                    cds.setCodingEnd(cds.getEnd());
                    exons.add(cds);
                } else {
                    exon.setCodingStart(cds.getStart());
                    exon.setCodingEnd(cds.getEnd());
                    exon.setReadingFrame(cds.getReadingShift());
                }
                cdss.remove(0);
            }
            for (Exon exon : exons) {
                chr = exon.getChr();
                strand = exon.getStrand();
                start = Math.min(exon.getStart(), start);
                end = Math.max(exon.getEnd(), end);
                name = exon.getName();
            }


            if (transcript == null) {
                // transcript is implied

                transcript = new BasicFeature(chr, start, end, strand);
                transcript.setIdentifier(id);
                transcript.setName(name == null ? id : name);
                transcript.setDescription(description);
                transcript.setAttributes(attributes);
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

        String getName(Map<String, String> attributes);

        void setNameFields(String[] fields);
    }

    public static class GFF2Helper implements Helper {

        //TODO Almost identical
        static String[] idFields = {"systematic_id", "ID", "transcript_id", "name", "primary_name", "gene", "locus", "alias"};
        static String[] DEFAULT_NAME_FIELDS = {"gene", "name", "primary_name", "locus", "alias", "systematic_id", "ID"};

        private String[] nameFields;

        GFF2Helper() {
            this(DEFAULT_NAME_FIELDS);
        }

        GFF2Helper(String[] nameFields) {
            if (nameFields != null) {
                this.nameFields = nameFields;
            }

        }

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

        public String getName(Map<String, String> attributes) {

            if (attributes.size() > 0 && nameFields != null) {
                for (String nf : nameFields) {
                    if (attributes.containsKey(nf)) {
                        return attributes.get(nf);
                    }
                }
            }

            return null;
        }

        public void setNameFields(String[] nameFields) {
            this.nameFields = nameFields;
        }
    }

    public static class GFF3Helper implements Helper {

        static String[] DEFAULT_NAME_FIELDS = {"Name", "Alias", "ID", "Gene", "gene", "Locus", "locus"};
        private boolean useUrlDecoding = true;

        private String[] nameFields;

        public GFF3Helper() {
            this(DEFAULT_NAME_FIELDS);
        }

        GFF3Helper(String[] nameFields) {
            if (nameFields != null) {
                this.nameFields = nameFields;
            }

        }


        public String[] getParentIds(Map<String, String> attributes, String ignored) {
            String parentIdString = attributes.get("Parent");
            if (parentIdString != null) {
                return attributes.get("Parent").split(",");
            } else {
                return null;
            }
        }

        /**
         * Parse the column 9 attributes.  Attributes are separated by semicolons.
         *
         * @param description
         * @param kvalues
         */
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
                        key = StringUtils.decodeURL(key);
                        value = StringUtils.decodeURL(value);
                        // Limit values to 50 characters
                        if (value.length() > 50) {
                            value = value.substring(0, 50) + " ...";
                        }
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

        public String getName(Map<String, String> attributes) {

            if (attributes.size() > 0 && nameFields != null) {
                for (String nf : nameFields) {
                    if (attributes.containsKey(nf)) {
                        return attributes.get(nf);
                    }
                }
            }

            return null;
        }

        public String getID(Map<String, String> attributes) {
            String id = attributes.get("ID");
            return id;
        }

        public String[] getNameFields() {
            return nameFields;
        }

        public void setNameFields(String[] nameFields) {
            this.nameFields = nameFields;
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
