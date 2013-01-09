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

package org.broad.igv.feature.tribble;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.CI;
import org.broad.igv.util.collections.MultiMap;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.*;

/**
 * Notes from GFF3 spec  http://www.sequenceontology.org/gff3.shtml
 * These tags have predefined meanings (tags are case sensitive):
 * <p/>
 * ID	   Indicates the name of the feature (unique).
 * Name   Display name for the feature.
 * Alias  A secondary name for the feature.
 * Parent Indicates the parent of the feature.
 * <p/>
 * Specs:
 * GFF3  http://www.sequenceontology.org/gff3.shtml
 * GFF2 specification: http://www.sanger.ac.uk/resources/software/gff/spec.html
 * UCSC GFF (GFF "1") http://genome.ucsc.edu/FAQ/FAQformat#format3
 * GTF  http://mblab.wustl.edu/GTF2.html
 * UCSC GTF  http://genome.ucsc.edu/FAQ/FAQformat#format4
 * Feature type definitions http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html#7.2
 */
public class GFFCodec extends AsciiFeatureCodec<Feature> {

    private static Logger log = Logger.getLogger(GFFCodec.class);

    public static Set<String> exonTerms = new HashSet();
    public static Set<String> utrTerms = new HashSet();
    public static Set<String> geneParts = new HashSet();
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


    private TrackProperties trackProperties = null;
    CI.CIHashSet featuresToHide = new CI.CIHashSet();


    FeatureFileHeader header;
    Helper helper;
    Genome genome;

    public enum Version {
        GFF2, GFF3
    }


    /**
     * List of know "Name" fields.  Some important fields from the GFF3 spec are listed below.  Note GFF3
     * is case sensitive, however GFF2, GTF, and other variants might not be.
     * <p/>
     * ID	  Indicates the ID of the feature.
     * Name   Display name for the feature.
     * Alias  A secondary name for the feature.
     */
    static String[] nameFields = {"Name", "name", "Alias", "gene", "primary_name", "locus", "alias", "systematic_id", "ID"};


    public GFFCodec(Genome genome) {
        super(Feature.class);
        // Assume GFF2 until shown otherwise
        helper = new GFF2Helper();
        this.genome = genome;
    }

    public GFFCodec(Version version, Genome genome) {
        super(Feature.class);
        this.genome = genome;
        if (version == Version.GFF2) {
            helper = new GFF2Helper();
        } else {
            helper = new GFF3Helper();
        }
    }

    public void readHeaderLine(String line) {
        if(header == null) {
            header = new FeatureFileHeader();
        }
        if (line.startsWith("#track") || line.startsWith("##track")) {
            trackProperties = new TrackProperties();
            ParsingUtils.parseTrackLine(line, trackProperties);
            header.setTrackProperties(trackProperties);
        } else if (line.startsWith("##gff-version") && line.endsWith("3")) {
            helper = new GFF3Helper();
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

    public Object readHeader(LineReader reader) {

        if(header == null) {
            header = new FeatureFileHeader();
        }
        header = new FeatureFileHeader();
        String line;
        int nLines = 0;
        try {
            while ((line = reader.readLine()) != null) {

                if (line.startsWith("#")) {
                    nLines++;
                    readHeaderLine(line);
                } else {
                    break;
                }
            }

            header.setTrackProperties(trackProperties);
            return header;
        } catch (IOException e) {
            throw new CodecLineParsingException("Error parsing header", e);
        }
    }

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(String path) {
        final String pathLowerCase = path.toLowerCase();
        return pathLowerCase.endsWith(".gff") || pathLowerCase.endsWith(".gff3") ||
                pathLowerCase.endsWith(".gvf");
    }

    public BasicFeature decodeLoc(String line) {
        return decode(line);
    }

    public BasicFeature decode(String line) {

        if (line.startsWith("#")) {
            // This should not be possible as this line would be parsed as a header.  But just in case
            return null;
        }

        String[] tokens = Globals.tabPattern.split(line, -1);
        int nTokens = tokens.length;

        // GFF3 files have 9 tokens,
        // TODO -- the attribute column is optional for GFF 2 and earlier (8 tokens required)
        if (nTokens < 9) {
            return null;
        }

        String chrToken = tokens[0].trim();
        String featureType = StringUtils.intern(tokens[2].trim());

        if (ignoredTypes.contains(featureType)) {
            return null;
        }

        String chromosome = genome == null ? StringUtils.intern(chrToken) : genome.getChromosomeAlias(chrToken);

        // GFF coordinates are 1-based inclusive (length = end - start + 1)
        // IGV (UCSC) coordinates are 0-based exclusive.  Adjust start and end accordingly
        int start;
        int end;
        int col = 3;
        try {
            start = Integer.parseInt(tokens[col]) - 1;
            if(start < 0) throw new ParserException("Start index must be 1 or larger; GFF is 1-based", -1, line);
            col++;
            end = Integer.parseInt(tokens[col]);
        } catch (NumberFormatException ne) {
            String msg = String.format("Column %d must contain a numeric value. %s", col + 1, ne.getMessage());
            throw new ParserException(msg, -1, line);
        }

        Strand strand = convertStrand(tokens[6]);
        String attributeString = tokens[8];

        //CI.CILinkedHashMap<String> attributes = new CI.CILinkedHashMap();
        MultiMap<String, String> attributes = new MultiMap<String, String>();

        helper.parseAttributes(attributeString, attributes);

        String description = getDescription(attributes, featureType);
        String id = helper.getID(attributes);
        String[] parentIds = helper.getParentIds(attributes, attributeString);

        if (exonTerms.contains(featureType) && parentIds != null && parentIds.length > 0 &&
                parentIds[0] != null && parentIds[0].length() > 0 && !parentIds[0].equals(".")) {

            //Somewhat tacky, but we need to store the phase somewhere in the feature
            String phaseString = tokens[7].trim();
            //String old = attributes.put(GFFFeatureSource.PHASE_STRING, phaseString);
            //if(old != null){
            //    log.debug("phase string attribute was overwritten internally; old value was: " + old);
            //}
        }

        BasicFeature f = new BasicFeature(chromosome, start, end, strand);

        f.setName(getName(attributes));
        f.setType(featureType);
        f.setDescription(description);

        id = id != null ? id : "igv_" + UUID.randomUUID().toString();
        f.setIdentifier(id);

        f.setParentIds(parentIds);
        f.setAttributes(attributes);

        String[] colorNames = new String[]{"color", "Color", "colour", "Colour"};
        for(String colorName: colorNames){
            if (attributes.containsKey(colorName)) {
                f.setColor(ColorUtilities.stringToColor(attributes.get(colorName)));
                break;
            }
        }

        if (featuresToHide.contains(featureType)) {
            if (IGV.hasInstance()) FeatureDB.addFeature(f);
            return null;
        }

        return f;

    }

    public Object getHeader() {
        return header;
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


    String getName(MultiMap<String, String> attributes) {

        if (attributes == null || attributes.size() == 0) {
            return null;
        }
        for (String nf : nameFields) {
            if (attributes.containsKey(nf)) {
                return attributes.get(nf);
            }
        }
        return "";
    }

    static StringBuffer buf = new StringBuffer();

    static String getDescription(MultiMap<String, String> attributes, String type) {
        buf.setLength(0);
        buf.append(type);
        buf.append("<br>");
        attributes.printHtml(buf, 100);
        return buf.toString();
    }


    protected interface Helper {

        String[] getParentIds(MultiMap<String, String> attributes, String attributeString);

        void parseAttributes(String attributeString, MultiMap<String, String> map);

        String getID(MultiMap<String, String> attributes);

        void setUrlDecoding(boolean b);

        String getName(MultiMap<String, String> attributes);

        void setNameFields(String[] fields);

    }

    public static class GFF2Helper implements Helper {

        //TODO Almost identical
        static String[] idFields = {"systematic_id", "ID", "transcript_id", "name", "primary_name", "gene", "locus", "alias"};
        static String[] DEFAULT_NAME_FIELDS = {"gene", "name", "primary_name", "locus", "alias", "systematic_id", "ID"};
        static String[] possParentNames = new String[]{"id", "mRna", "systematic_id", "transcript_id", "gene", "transcriptId", "Parent", "proteinId"};

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
            // Ignored,  GFF2 files are never url DECODED
        }


        public void parseAttributes(String description, MultiMap<String, String> kvalues) {

            List<String> kvPairs = StringUtils.breakQuotedString(description.trim(), ';');
            for (String kv : kvPairs) {
                String[] tokens = kv.split(" ");
                if (tokens.length == 1) {
                    //Not space delimited, check =
                    tokens = kv.split("=");
                }
                if (tokens.length >= 2) {
                    String key = tokens[0].trim().replaceAll("\"", "");
                    String value = tokens[1].trim().replaceAll("\"", "");
                    kvalues.put(StringUtils.intern(key), value);
                }
            }
        }

        /**
         * parentIds[0] = attributes.get("id");
         * if (parentIds[0] == null) {
         * parentIds[0] = attributes.get("mRNA");
         * }
         * if (parentIds[0] == null) {
         * parentIds[0] = attributes.get("systematic_id");
         * }
         * if (parentIds[0] == null) {
         * parentIds[0] = attributes.get("transcript_id");
         * }
         * if (parentIds[0] == null) {
         * parentIds[0] = attributes.get("gene");
         * }
         * if (parentIds[0] == null) {
         * parentIds[0] = attributes.get("transcriptId");
         * }
         * if (parentIds[0] == null) {
         * parentIds[0] = attributes.get("proteinId");
         * }
         *
         * @param attributes
         * @param attributeString
         * @return
         */

        public String[] getParentIds(MultiMap<String, String> attributes, String attributeString) {

            String[] parentIds = new String[1];
            if (attributes.size() == 0) {
                parentIds[0] = attributeString;
            } else {
                for (String possName : possParentNames) {
                    if (attributes.containsKey(possName)) {
                        parentIds[0] = attributes.get(possName);
                        break;
                    }
                }
            }
            return parentIds;
        }


        public String getID(MultiMap<String, String> attributes) {
            for (String nf : idFields) {
                if (attributes.containsKey(nf)) {
                    return attributes.get(nf);
                }
            }
            return getName(attributes);
        }

        public String getName(MultiMap<String, String> attributes) {

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

        static String[] DEFAULT_NAME_FIELDS = {"Name", "Alias", "ID", "gene", "locus"};
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


        public String[] getParentIds(MultiMap<String, String> attributes, String ignored) {
            String parentIdString = attributes.get("Parent");
            if (parentIdString != null) {
                return parentIdString.split(",");
            } else {
                return null;
            }
        }

        /**
         * Parse the column 9 attributes.  Attributes are separated by semicolons.
         * <p/>
         * TODO -- quotes (column 9) are explicitly forbidden in GFF3 -- should breakQuotedString be used?
         *
         * @param description
         * @param kvalues
         */
        public void parseAttributes(String description, MultiMap<String, String> kvalues) {

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
                    kvalues.put(StringUtils.intern(key), value);
                } else {
                    log.info("No attributes: " + description);
                }
            }
        }

        public void setUrlDecoding(boolean useUrlDecoding) {
            this.useUrlDecoding = useUrlDecoding;
        }

        public String getName(MultiMap<String, String> attributes) {

            if (attributes.size() > 0 && nameFields != null) {
                for (String nf : nameFields) {
                    if (attributes.containsKey(nf)) {
                        return attributes.get(nf);
                    }
                }
            }

            return null;
        }

        public String getID(MultiMap<String, String> attributes) {
            return attributes.get("ID");
        }

        public void setNameFields(String[] nameFields) {
            this.nameFields = nameFields;
        }
    }


    /**
     * Helper for GTF files
     * <p/>
     * mandatory attributes
     * gene_id value;     A globally unique identifier for the genomic source of the transcript
     * transcript_id value;     A globally unique identifier for the predicted transcript.
     * <p/>
     * Attributes must end in a semicolon which must then be separated from the start of any subsequent
     * attribute by exactly one space character (NOT a tab character).
     * <p/>
     * Textual attributes should be surrounded by doublequotes.
     */
    public static class GTFHelper implements Helper {

        @Override
        public void parseAttributes(String description, MultiMap<String, String> kvalues) {

            List<String> kvPairs = StringUtils.breakQuotedString(description.trim(), ';');
            for (String kv : kvPairs) {
                List<String> tokens = StringUtils.breakQuotedString(kv, ' ');
                if (tokens.size() >= 2) {
                    String key = tokens.get(0).trim().replaceAll("\"", "");
                    String value = tokens.get(1).trim().replaceAll("\"", "");
                    kvalues.put(StringUtils.intern(key), value);
                }
            }
        }

        @Override
        public String[] getParentIds(MultiMap<String, String> attributes, String attributeString) {
            return new String[0];  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public String getID(MultiMap<String, String> attributes) {
            return attributes.get("transcript_id");  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public void setUrlDecoding(boolean b) {
            //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public String getName(MultiMap<String, String> attributes) {
            return attributes.get("transcript_id");  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public void setNameFields(String[] fields) {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }

}
