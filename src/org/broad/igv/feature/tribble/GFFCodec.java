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

package org.broad.igv.feature.tribble;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.SequenceOntology;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.CI;
import org.broad.igv.util.collections.MultiMap;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.exception.CodecLineParsingException;
import htsjdk.tribble.readers.LineIterator;

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


    static HashSet<String> ignoredTypes = new HashSet();

    static {
        ignoredTypes.add("start_codon");
        ignoredTypes.add("stop_codon");
        ignoredTypes.add("Contig");
        ignoredTypes.add("RealContig");
        ignoredTypes.add("CDS_parts");
    }


    private TrackProperties trackProperties = null;
    private CI.CIHashSet featuresToHide = new CI.CIHashSet();
    private FeatureFileHeader header;
    private Helper helper;
    private Genome genome;
    private boolean fastaSection = false;

    public enum Version {
        GFF2, GFF3
    }


    /**
     * List of known "Name" fields.  Some important fields from the GFF3 spec are listed below.  Note GFF3
     * is case sensitive, however GFF2, GTF, and other variants might not be.
     * <p/>
     * ID	  Indicates the ID of the feature.
     * Name   Display name for the feature.
     * Alias  A secondary name for the feature.
     */
    static String[] nameFields = {"Name", "name", "Alias", "gene", "primary_name", "locus", "alias", "systematic_id", "ID", "transcript_id"};


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
        if (header == null) {
            header = new FeatureFileHeader();
        }
        if (line.startsWith("#track") || line.startsWith("##track")) {
            trackProperties = new TrackProperties();
            ParsingUtils.parseTrackLine(line, trackProperties);
            header.setTrackProperties(trackProperties);
        } else if (line.startsWith("##gff-version") && line.contains("3")) {
            String[] tokens = Globals.whitespacePattern.split(line);
            if (tokens.length > 1 && tokens[1].startsWith("3")) {
                helper = new GFF3Helper();
            }
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

    public Object readActualHeader(LineIterator reader) {

        if (header == null) {
            header = new FeatureFileHeader();
        }
        String line;
        int nLines = 0;
        try {
            while (reader.hasNext()) {
                line = reader.peek();
                if (line.startsWith("#")) {
                    nLines++;
                    readHeaderLine(line);
                    reader.next();
                } else {
                    break;
                }
            }

            header.setTrackProperties(trackProperties);
            return header;
        } catch (Exception e) {
            throw new CodecLineParsingException("Error parsing header: " + e.getMessage(), e);
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
        String pathLowerCase = path.toLowerCase();
        if (pathLowerCase.endsWith(".gz")) pathLowerCase = pathLowerCase.substring(0, pathLowerCase.length() - 3);
        return pathLowerCase.endsWith(".gff") || pathLowerCase.endsWith(".gff3") ||
                pathLowerCase.endsWith(".gvf") || pathLowerCase.endsWith(".gtf");
    }

    public BasicFeature decodeLoc(String line) {
        return decode(line);
    }

    public BasicFeature decode(String line) {

        if (fastaSection) {
            return null;
        }
        if (line.startsWith("#")) {
            if (line.toUpperCase().startsWith("##FASTA")) {
                fastaSection = true;
            }
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

        String chromosome = genome == null ? StringUtils.intern(chrToken) : genome.getCanonicalChrName(chrToken);

        // GFF coordinates are 1-based inclusive (length = end - start + 1)
        // IGV (UCSC) coordinates are 0-based exclusive.  Adjust start and end accordingly
        int start;
        int end;
        int col = 3;
        try {
            start = Integer.parseInt(tokens[col]) - 1;
            if (start < 0) throw new ParserException("Start index must be 1 or larger; GFF is 1-based", -1, line);
            col++;
            end = Integer.parseInt(tokens[col]);
        } catch (NumberFormatException ne) {
            String msg = String.format("Column %d must contain a numeric value. %s", col + 1, ne.getMessage());
            throw new ParserException(msg, -1, line);
        }
        Strand strand = convertStrand(tokens[6]);

        String attributeString = tokens[8];

        MultiMap<String, String> attributes = new MultiMap<String, String>();

        helper.parseAttributes(attributeString, attributes);

        String id = helper.getID(attributes, featureType);
        String[] parentIds = helper.getParentIds(attributes, attributeString);

        BasicFeature f = new BasicFeature(chromosome, start, end, strand);


        // Set "thick start/end" => corresponds to coding start & end, for UTRs
        if (SequenceOntology.utrTypes.contains(featureType)) {
            boolean plus = (SequenceOntology.fivePrimeUTRTypes.contains(featureType) && strand == Strand.POSITIVE) ||
                    (SequenceOntology.threePrimeUTRTypes.contains(featureType) && strand == Strand.NEGATIVE);
            if (plus) {
                f.setThickStart(end);
            } else {
                f.setThickEnd(end);
            }
        }

        String phaseString = tokens[7].trim();
        if (!phaseString.equals(".")) {
            int phaseNum = Integer.parseInt(phaseString);
            f.setReadingFrame(phaseNum);
        }

        f.setName(helper.getName(attributes));
        f.setType(featureType);

        id = id != null ? id : "igv_" + UUID.randomUUID().toString();
        f.setIdentifier(id);

        f.setParentIds(parentIds);
        f.setAttributes(attributes);

        String[] colorNames = new String[]{"color", "Color", "colour", "Colour"};
        for (String colorName : colorNames) {
            if (attributes.containsKey(colorName)) {
                f.setColor(ColorUtilities.stringToColor(attributes.get(colorName)));
                break;
            }
        }

        if (featuresToHide.contains(featureType)) {
            if (IGV.hasInstance()) FeatureDB.addFeature(f, genome);
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

        String getID(MultiMap<String, String> attributes, String type);

        void setUrlDecoding(boolean b);

        String getName(MultiMap<String, String> attributes);

        void setNameFields(String[] fields);

    }

    public static class GFF2Helper implements Helper {

        //TODO Almost identical
        static String[] DEFAULT_NAME_FIELDS = {"alias", "gene", "ID", "Locus", "locus", "Name", "name", "gene_name", "primary_name", "systematic_id", "transcript_id"};
        static List<String> idFields = new ArrayList<String>(Arrays.asList(DEFAULT_NAME_FIELDS));

        static {
            idFields.add("transcript_id");
        }

        static String[] possParentNames = new String[]{"transcript_id", "id", "mRNA", "systematic_id", "gene", "transcriptId", "Parent", "proteinId"};

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
         * @param attributes
         * @param attributeString
         * @return
         */

        public String[] getParentIds(MultiMap<String, String> attributes, String attributeString) {

            if (attributes.size() > 0) {
                for (String possName : possParentNames) {
                    if (attributes.containsKey(possName)) {
                        String parent = attributes.get(possName).trim();
                        if (parent.length() > 0) {
                            return new String[]{parent};
                        }
                    }
                }
            }
            return null;
        }


        public String getID(MultiMap<String, String> attributes, String type) {

            //Search for an attribute == type,  take this as ID
            String id = attributes.get(type);
            if (id != null && id.length() > 0) {
                return id;
            }

            for (String nf : idFields) {
                if (attributes.containsKey(nf)) {
                    String tmp = attributes.get(nf).trim();
                    if (tmp.length() > 0) return tmp;
                }
            }

            String tmp = getName(attributes);
            if (tmp != null && tmp.trim().length() > 0) {
                return tmp.trim();
            }

            return null;
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

        static String[] DEFAULT_NAME_FIELDS = {"Name", "Alias", "ID", "gene", "locus", "gene_name"};
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

        public String getID(MultiMap<String, String> attributes, String ignore) {
            return attributes.get("ID");
        }

        public void setNameFields(String[] nameFields) {
            this.nameFields = nameFields;
        }
    }


}
