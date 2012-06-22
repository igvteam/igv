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
import org.broad.igv.track.GFFFeatureSource;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.CI;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.*;

/**
 * Notes from GFF3 spec
 * These tags have predefined meanings:
 * <p/>
 * ID	   Indicates the name of the feature.  IDs must be unique
 * within the scope of the GFF file.
 * <p/>
 * Name   Display name for the feature.  This is the name to be
 * displayed to the user.  Unlike IDs, there is no requirement
 * that the Name be unique within the file.
 * <p/>
 * Alias  A secondary name for the feature.  It is suggested that
 * this tag be used whenever a secondary identifier for the
 * feature is needed, such as locus names and
 * accession numbers.  Unlike ID, there is no requirement
 * that Alias be unique within the file.
 * <p/>
 * Parent Indicates the parent of the feature.  A parent ID can be
 * used to group exons into transcripts, transcripts into
 * genes, an so forth.  A feature may have multiple parents.
 * Parent can *only* be used to indicate a partof
 * relationship.
 */
public class GFFCodec extends AsciiFeatureCodec<Feature> {

    private static Logger log = Logger.getLogger(GFFCodec.class);

    public static CI.CIHashSet exonTerms = new CI.CIHashSet();
    public static CI.CIHashSet utrTerms = new CI.CIHashSet();
    static CI.CIHashSet geneParts = new CI.CIHashSet();
    static CI.CIHashSet ignoredTypes = new CI.CIHashSet();

    static {
        utrTerms.add("five_prime_UTR");
        utrTerms.add("three_prime_UTR");
        utrTerms.add("5'-utr");
        utrTerms.add("3'-utr");
        utrTerms.add("5utr");
        utrTerms.add("3utr");
    }

    static {
        exonTerms.addAll(utrTerms);
        exonTerms.add("exon");
        exonTerms.add("coding_exon");
        exonTerms.add("CDS");
    }


    static {
        geneParts.addAll(exonTerms);
        geneParts.add("transcript");
        geneParts.add("processed_transcript");
        geneParts.add("mrna");
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

    static String[] nameFields = {"name", "gene", "primary_name", "locus",
            "alias", "systematic_id", "ID"};


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

    public void readHeaderLine(String line){
        header = new FeatureFileHeader();
        if (line.startsWith("#track") || line.startsWith("##track")) {
            trackProperties = new TrackProperties();
            ParsingUtils.parseTrackLine(line, trackProperties);
            header.setTrackProperties(trackProperties);
        } else if (line.startsWith("##gff-version") && line.endsWith("3")) {
            helper = new GFF3Helper();
        }else if (line.startsWith("#nodecode") || line.startsWith("##nodecode")) {
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


        if (line.startsWith("##gff-version") && line.endsWith("3")) {
            helper = new GFF3Helper();
        }

        if (line.startsWith("#")) {
            return null;
        }

        String[] tokens = Globals.tabPattern.split(line, -1);
        int nTokens = tokens.length;

        // GFF files have 9 tokens
        if (nTokens < 9) {
            return null;
        }

        String chrToken = tokens[0].trim();
        String featureType = tokens[2].trim();
        String chromosome = genome == null ? chrToken : genome.getChromosomeAlias(chrToken);

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
            throw new ParserException(msg, -1, line);
        }

        Strand strand = convertStrand(tokens[6]);
        String attributeString = tokens[8];

        CI.CILinkedHashMap<String> attributes = new CI.CILinkedHashMap();

        helper.parseAttributes(attributeString, attributes);

        String description = getDescription(attributes, featureType);
        String id = helper.getID(attributes);
        String[] parentIds = helper.getParentIds(attributes, attributeString);

        if (exonTerms.contains(featureType) && parentIds != null && parentIds.length > 0 &&
                parentIds[0] != null && parentIds[0].length() > 0 && !parentIds[0].equals(".")) {

            //Somewhat tacky, but we need to store the phase somewhere in the feature
            String phaseString = tokens[7].trim();
            String old = attributes.put(GFFFeatureSource.PHASE_STRING, phaseString);
            if(old != null){
                log.debug("phase string attribute was overwritten internally; old value was: " + old);
            }
        }

        BasicFeature f = new BasicFeature(chromosome, start, end, strand);

        f.setName(getName(attributes));
        f.setType(featureType);
        f.setDescription(description);
        f.setIdentifier(id);
        f.setParentIds(parentIds);
        f.setAttributes(attributes);

        if (attributes.containsKey("color")) {
            f.setColor(ColorUtilities.stringToColor(attributes.get("color")));
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


    String getName(Map<String, String> attributes) {

        if (attributes == null || attributes.size() == 0) {
            return null;
        }
        for (String nf : nameFields) {
            if (attributes.containsKey(nf)) {
                return attributes.get(nf);
            }
        }

        // If still nothing return the first attribute value
        return attributes.values().iterator().next();
    }

    static StringBuffer buf = new StringBuffer();

    static String getDescription(Map<String, String> attributes, String type) {
        buf.setLength(0);
        buf.append(type);
        buf.append("<br>");
        for (Map.Entry<String, String> att : attributes.entrySet()) {
            String attValue = att.getValue().replaceAll(";", "<br>");
            buf.append(att.getKey());
            buf.append(" = ");
            buf.append(attValue);
            buf.append("<br>");
        }

        String description = buf.toString();

        return description;
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
            // Ignored,  GFF2 files are never url DECODED
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
                String[] possNames = new String[]{"id", "mrna", "systematic_id", "transcript_id", "gene", "transcriptid", "proteinid"};
                for(String possName: possNames){
                    if(attributes.containsKey(possName)){
                        parentIds[0] = attributes.get(possName);
                        break;
                    }
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


        public String[] getParentIds(Map<String, String> attributes, String ignored) {
            String parentIdString = attributes.get("Parent");
            if (parentIdString != null) {
                return parentIdString.split(",");
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
            return attributes.get("ID");
        }

        public void setNameFields(String[] nameFields) {
            this.nameFields = nameFields;
        }
    }

}
