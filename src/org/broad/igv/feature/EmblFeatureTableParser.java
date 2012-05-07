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

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.renderer.GeneTrackRenderer;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.util.*;

/**
 * Parses exon information from the feature table section of an embl file.
 * <p/>
 * Note -- this is no a general parser,  it is built specifically to handle the S Pombe annotation
 * from Sanger.  Example entry follows
 * <p/>
 * FT   CDS             join(180423..180528,180586..180721,180890..180911)
 * FT                   /product="TIM22 inner membrane protein import complex
 * FT                   subunit Tim8 (predicted)"
 * FT                   /gene="tim8"
 * FT                   /gene="SPAC13G6.04"
 */
public class EmblFeatureTableParser implements FeatureParser {

    /**
     * Method description
     *
     * @param gene
     */
    public static void computeReadingShifts(IGVFeature gene) {
        List<org.broad.igv.feature.Exon> exons = gene.getExons();
        if (exons.size() == 0) {
            return;
        }

        int startIndex = (gene.getStrand() == Strand.POSITIVE) ? 0 : exons.size() - 1;
        int endIndex = (gene.getStrand() == Strand.POSITIVE) ? exons.size() : -1;
        int increment = (gene.getStrand() == Strand.POSITIVE) ? 1 : -1;
        int cds = 0;
        int exonNumber = 1;
        for (int i = startIndex; i != endIndex; i += increment) {
            org.broad.igv.feature.Exon exon = exons.get(i);
            exon.setNumber(exonNumber);
            int modCds = cds % 3;
            int phase = (modCds == 0) ? 0 : 3 - modCds;
            exon.setPhase(phase);
            cds += exon.getCodingLength();
            exonNumber++;
        }
    }

    public TrackProperties getTrackProperties() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Method description
     *
     * @param locator
     * @return
     */
    public List<FeatureTrack> loadTracks(ResourceLocator locator, Genome genome) {

        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(locator);

            List<org.broad.tribble.Feature> features = loadFeatures(reader, genome);

            if (features.isEmpty()) {
                return null;
            } else {
                FeatureTrack track = new FeatureTrack(locator, new FeatureCollectionSource(features, genome));
                track.setName(locator.getTrackName());
                track.setRendererClass(IGVFeatureRenderer.class);
                track.setMinimumHeight(35);
                track.setHeight(45);
                track.setRendererClass(GeneTrackRenderer.class);

                List<FeatureTrack> newTracks = new ArrayList();

                newTracks.add(track);

                FeatureDB.addFeatures(features);
                return newTracks;
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {

                }
            }
        }
    }

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    public List<org.broad.tribble.Feature> loadFeatures(BufferedReader reader, Genome genome) {

        List<BasicFeature> features = new ArrayList();

        String nextLine = null;
        String chromosome = null;
        EmblRecord currentRecord = null;
        try {
            while ((nextLine = reader.readLine()) != null) {

                if (nextLine.startsWith("ID"))    // Chromosome change
                {
                    String chr = getFirstWord(nextLine.substring(2));
                    chromosome = chr.replace("chromosome", "chr").replace("_", "");
                } else if (nextLine.startsWith("FT")) {
                    String featureKey = nextLine.substring(5, 19).trim();
                    if (featureKey.length() == 0) {
                        if (currentRecord != null) {
                            currentRecord.append(nextLine);
                        }
                    } else {

                        // New feature started.
                        if ((currentRecord != null)) {
                            features.add(createFeature(currentRecord));
                        }

                        String temp = nextLine.substring(21);
                        boolean isNegative = temp.contains("complement");
                        String lociString = parseJoinString(temp, reader).replace("<",
                                "").replace(">", "").trim();
                        currentRecord = new EmblRecord(featureKey.trim(), chromosome, lociString,
                                isNegative);
                    }


                }
            }

            // features.addAll(convertEmblToGenes(emblGenes));
            return combineGeneParts(features);
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }

    static BasicFeature createFeature(EmblRecord emblRecord) {
        BasicFeature feature = new BasicFeature(emblRecord.getChromosome(), emblRecord.getStart(),
                emblRecord.getEnd());
        feature.setType(emblRecord.getType());
        feature.setIdentifier(emblRecord.getIdentifier());
        feature.setName(emblRecord.getIdentifier());
        feature.setStrand(emblRecord.getStrand());
        feature.setDescription(emblRecord.getDescription());
        if (emblRecord.getAlias() != null) {
            feature.setName(emblRecord.getAlias());
        }

        // If this is a "gene part" add the exons
        for (Exon exon : emblRecord.getExons()) {
            feature.addExon(exon);
        }

        return feature;

    }

    static HashSet<String> geneParts = new HashSet();

    static {
        geneParts.add("CDS");
        geneParts.add("5'UTR");
        geneParts.add("3'UTR");
        geneParts.add("intron");
    }

    // Combine the constituitve gene parts into a single feature (3', 5', and cds).  This is
    // necccessary with the current IGV gene model, which is heavily influence by ucsc conventions
    static List<org.broad.tribble.Feature> combineGeneParts(List<BasicFeature> features) {

        List<org.broad.tribble.Feature> newFeatureList = new ArrayList(features.size());
        Map<String, BasicFeature> genes = new HashMap();

        // Step 1 -- find the coding regions

        for (BasicFeature feature : features) {
            if (feature.getType().equals("CDS")) {
                String geneId = feature.getIdentifier();
                BasicFeature gene = genes.get(geneId);
                if (gene == null) {
                    gene = feature.copy();
                    gene.setType("transcript");
                    newFeatureList.add(gene);
                    genes.put(geneId, gene);
                }
                Exon exon = new Exon(feature.getChr(), feature.getStart(), feature.getEnd(),
                        feature.getStrand());
                gene.addExon(exon);
            }
        }

        // Step 2 -- add the UTRs and non-gene features
        for (BasicFeature feature : features) {
            String type = feature.getType();
            if (type.equals("CDS")) {

                // already accounted for
            } else if (type.equals("3'UTR") || type.equals("5'UTR")) {
                BasicFeature gene = genes.get(feature.getIdentifier());
                if (gene != null) {
                    Exon exon = new Exon(feature.getChr(), feature.getStart(),
                            feature.getEnd(), feature.getStrand());
                    exon.setUTR(true);
                    gene.addExon(exon);
                }

            } else {
                newFeatureList.add(feature);
            }
        }

        // Compute the phases (reading shifts) for the genes
        // TODO -- should we be doing this?  It generally works but could there be reading shifts
        // between exons that would throw this off?
        for (IGVFeature gene : genes.values()) {
            computeReadingShifts(gene);
        }

        return newFeatureList;

    }

    /**
     * FT   CDS             join(complement(5000933..5001976),
     * FT                   complement(5000325..5000891),complement(5000024..5000272))
     * FT                   /product="GTPase activating protein (predicted)"
     * FT                   /gene="SPAC1952.17c"
     * FT                   /gene="SPAC890.01c"
     *
     * @param joinString
     * @param reader
     * @return
     * @throws IOException
     */
    public static String parseJoinString(String joinString, BufferedReader reader)
            throws IOException {

        if (joinString.startsWith("join") || joinString.startsWith("complement")) {
            int leftParenCount = countChar(joinString, '(');
            int rightParenCount = countChar(joinString, ')');
            while (leftParenCount != rightParenCount) {
                joinString += reader.readLine().replace("FT", "").trim();
                leftParenCount = countChar(joinString, '(');
                rightParenCount = countChar(joinString, ')');
            }

            // join and complement functions irrelevant
            joinString = joinString.replace("join", "");
            joinString = joinString.replace("complement", "");
            joinString = joinString.replace("(", "");
            joinString = joinString.replace(")", "");
            joinString = joinString.replace('<', ' ');
            return joinString;

        } else {
            return joinString;
        }

    }

    /**
     * This must exist in the jdk ?
     *
     * @param string
     * @return
     */
    static int countChar(String string, char c) {
        int cnt = 0;
        for (int i = 0; i < string.length(); i++) {
            if (c == string.charAt(i)) {
                cnt++;
            }
        }
        return cnt;

    }

    static class EmblRecord {

        private static Logger log = Logger.getLogger(EmblRecord.class);

        boolean isNegative;
        private String type;
        private String chromosome;
        private String identifier;
        private String alias;
        private String description;
        private int start = Integer.MAX_VALUE;
        private int end;
        List<Exon> exons;


        EmblRecord(String type, String chromosome, String lociString, boolean isNegative) {
            this.isNegative = isNegative;
            this.type = type;
            this.chromosome = chromosome;
            createExons(lociString, isNegative);
        }

        /**
         * Method description
         *
         * @return
         */
        public int getStart() {
            return start;
        }

        /**
         * Method description
         *
         * @return
         */
        public int getEnd() {
            return end;
        }

        /**
         * Method description
         *
         * @return
         */
        public boolean isGenePart() {
            return type.equals("CDS") || type.equals("3'UTR") || type.equals("5'UTR");
        }

        /**
         * Method description
         *
         * @return
         */
        public Strand getStrand() {
            return isNegative ? Strand.NEGATIVE : Strand.POSITIVE;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getType() {
            return type;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getIdentifier() {
            return identifier;
        }

        /**
         * Method description
         *
         * @param identifier
         */
        public void setIdentifier(String identifier) {
            this.identifier = identifier;
        }

        /**
         * Method description
         *
         * @return
         */
        public String getAlias() {
            return alias;
        }

        /**
         * Method description
         *
         * @param alias
         */
        public void setAlias(String alias) {
            this.alias = alias;
        }

        /**
         * Method description
         *
         * @return
         */
        public List<Exon> getExons() {
            return exons;
        }

        /**
         * Method description
         *
         * @param nextLine
         */
        public void append(String nextLine) {
            String attrString = nextLine.substring(21);
            if (attrString.startsWith("/gene=")) {
                String[] kv = attrString.split("=");
                String geneName = kv[1].replace("\"", "");
                if (geneName.startsWith("SP")) {

                    // Some genes have multiple identifiers.  Only use the first one
                    if (getIdentifier() == null) {
                        setIdentifier(geneName);
                    }
                } else {
                    setAlias(geneName);
                }
            } else if (attrString.startsWith("/systematic_id=")) {
                String[] kv = attrString.split("=");
                String id = kv[1].replace("\"", "");
                setIdentifier(id);
                setAlias(id);

            } else {
                appendToDescription(nextLine.substring(22).trim());
            }
        }

        /**
         * Method description
         *
         * @param note
         */
        public void appendToDescription(String note) {
            if (description == null) {
                description = note;
            } else {
                description += "<br>" + note;
            }
        }

        /**
         * Method description
         *
         * @return
         */
        public String getDescription() {
            return description;
        }

        /**
         * Create a list of Exon objects from the Embl join string.  Apparently exons in embl
         * format are represented by a single CDS record.
         *
         * @param joinString
         * @param isNegative
         */
        void createExons(String joinString, boolean isNegative) {
            String[] lociArray = joinString.split(",");
            exons = new ArrayList(lociArray.length);
            for (String loci : lociArray) {
                try {
                    String[] tmp = loci.split("\\.\\.");
                    int exonStart = Integer.parseInt(tmp[0]) - 1;    // - (isNegative ? 0 : 1);
                    int exonEnd = exonStart + 1;
                    if (tmp.length > 1) {
                        exonEnd = Integer.parseInt(tmp[1]);
                    }

                    Strand strand = isNegative ? Strand.NEGATIVE : Strand.POSITIVE;
                    Exon r = new Exon(chromosome, exonStart, exonEnd, strand);
                    start = Math.min(start, exonStart);
                    end = Math.max(end, exonEnd);
                    exons.add(r);


                } catch (NumberFormatException e) {
                    log.error("Error parsing exon number; " + joinString, e);
                }
            }
        }

        /**
         * Method description
         *
         * @return
         */
        public String getChromosome() {
            return chromosome;
        }
    }

    static void combineFeatureTables(File[] emblFiles, String outputFile) throws IOException {

        PrintWriter pw = new PrintWriter(new FileWriter(outputFile));
        for (File ifile : emblFiles) {
            BufferedReader br = new BufferedReader(new FileReader(ifile));
            String nextLine = null;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("ID") || nextLine.startsWith("FT")) {
                    pw.println(nextLine);
                }
            }
            br.close();
        }
        pw.close();
    }

    static void splitByType(String emblFile, String outputDirectory) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(emblFile));
        String nextLine;
        Set<String> codes = new HashSet();

        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("FT") && (nextLine.length() > 19)) {
                String code = nextLine.substring(5, 19).trim();
                if (code.length() > 0) {
                    codes.add(code);
                }
            }
        }
        br.close();

        Map<String, PrintWriter> writers = new HashMap();
        for (String code : codes) {
            writers.put(code,
                    new PrintWriter(new FileWriter(new File(outputDirectory, code + ".embl"))));
        }

        br = new BufferedReader(new FileReader(emblFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("ID")) {
                for (PrintWriter pw : writers.values()) {
                    pw.println(nextLine);
                }
            } else if (nextLine.startsWith("FT")) {
                String code = nextLine.substring(5, 19).trim();
                if (code.length() > 0) {
                    currentWriter = writers.get(code);
                }
                if (currentWriter != null) {
                    currentWriter.println(nextLine);
                }
            } else {
                currentWriter = null;
            }
        }

        br.close();
        for (PrintWriter pw : writers.values()) {
            pw.close();
        }
    }

    /**
     * Utility method to extract the gene records from an embl file
     *
     * @param emblFile
     * @param outputFile
     * @throws java.io.IOException
     */
    static void extractGenes(String emblFile, String outputFile) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(emblFile));
        String nextLine;

        br = new BufferedReader(new FileReader(emblFile));
        PrintWriter geneWriter = new PrintWriter(new FileWriter(outputFile));
        PrintWriter currentWriter = null;
        while ((nextLine = br.readLine()) != null) {
            if (nextLine.startsWith("ID")) {
                geneWriter.println(nextLine);
            } else if (nextLine.startsWith("FT")) {
                String code = nextLine.substring(5, 19).trim();
                if (code.equals("CDS") || code.equals("3'UTR") || code.equals("5'UTR")) {

                    currentWriter = geneWriter;
                } else if (code.length() > 0) {
                    currentWriter = null;
                }
                if (currentWriter != null) {
                    currentWriter.println(nextLine);
                }
            } else {
                currentWriter = null;
            }
        }

        br.close();
        geneWriter.close();
    }


    // TODO -- put this in a utility class
    private static String getFirstWord(String string) {
        String trimmedString = string.trim();
        char[] chars = trimmedString.toCharArray();
        int whitespaceIndex = 0;
        for (whitespaceIndex = 0; whitespaceIndex < chars.length; whitespaceIndex++) {
            if (Character.isSpaceChar(chars[whitespaceIndex])) {
                break;
            }
        }

        return trimmedString.substring(0, whitespaceIndex).trim();

    }

    /**
     * Method description
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        // String[] ifiles = { "/Users/jrobinso/IGVTestData/SPombe/chromosome1.contig.embl",
        // "/Users/jrobinso/IGVTestData/SPombe/chromosome2.contig.embl",
        // "/Users/jrobinso/IGVTestData/SPombe/chromosome3.contig.embl" };
        // splitByType("/Users/jrobinso/IGVTestData/SPombe/spombe.ft.embl",
        // "/Users/jrobinso/IGVTestData/SPombe/");

        // combineFeatureTables((new File("/Users/jrobinso/plasmodium/v2.1/Embl/")).listFiles(), "/Users/jrobinso/plasmodium/v2.1/Embl/MAL.embl");
        // extractGenes("/Users/jrobinso/plasmodium/v2.1/Embl/MAL.embl",
        // "/Users/jrobinso/plasmodium/v2.1/Embl/MAL.genes.embl");
        splitByType("/Users/jrobinso/plasmodium/v2.1/Embl/MAL.embl",
                "/Users/jrobinso/plasmodium/v2.1/Embl");

        /*
        *
        * String fnTemp =
        *   "/Users/jrobinso/IGVTestData/SPombe/chromosome$.contig.embl";
        *
        * List<IGVFeature> allFeatures = new ArrayList();
        *
        * for (int i = 1; i <= 3; i++)
        * {
        *   String chromosome = "chr" + i;
        *   String fn = fnTemp.replace("$", "" + i);
        *
        *   List<IGVFeature> features = (new EmblFeatureTableParser()).loadTracks(fn);
        *   for (IGVFeature f : features)
        *   {
        *       ((BasicFeature) f).setChromosome(chromosome);
        *       ((BasicFeature) f).sortExons();
        *       if (f.getIdentifier().equals("SPAC212.04c"))
        *       {
        *           f.getAminoAcidSequence(1);
        *       }
        *       allFeatures.add(f);
        *   }
        * }
        *
        * // SPAC212.04c
        *
        *
        * AbstractFeatureFileParser.dumpFeatures(allFeatures,
        *       "/Users/jrobinso/spombe_160708.refflat");
        *
        *
        * String test1 = "<5000933..5001976";
        *
        * System.out.println(test1.replace("<", ""));
        *
        * BufferedReader br = new BufferedReader(new StringReader(test1));
        * System.out.println(parseJoinString(test1, br));
        *
        * String test2 =
        *   "join(5000933..5001976,5000325..5000891,5000024..5000272)";
        *
        * br = new BufferedReader(new StringReader(test2));
        * System.out.println(parseJoinString(test2, br));
        */
    }
}
