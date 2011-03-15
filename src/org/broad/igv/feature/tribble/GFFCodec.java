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

package org.broad.igv.feature.tribble;

import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.StringUtils;
import org.broad.tribble.*;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;


import java.io.IOException;
import java.net.URLDecoder;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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
public class GFFCodec implements org.broad.tribble.FeatureCodec {

    private static Logger log = Logger.getLogger(GFFCodec.class);

    public enum Version {
        GFF2, GFF3
    }

    static String[] nameFields = {"Name", "name", "gene", "primary_name", "Locus", "locus",
            "alias", "systematic_id", "ID"};


    FeatureFileHeader header;

    Helper helper;

    String[] tokens = new String[10];

    public GFFCodec() {
        // Assume GFF2 until shown otherwise
        helper = new GFF2Helper();
    }

    public GFFCodec(Version version) {
        if (version == Version.GFF2) {
            helper = new GFF2Helper();
        } else {
            helper = new GFF3Helper();
        }
    }


    public Object readHeader(LineReader reader) {


        header = new FeatureFileHeader();
        String line;
        TrackProperties trackProperties;
        int nLines = 0;

        try {
            while ((line = reader.readLine()) != null) {

                if (line.startsWith("#")) {
                    nLines++;
                    if (line.startsWith("##gff-version") && line.endsWith("3")) {
                        helper = new GFF3Helper();
                    } else if (line.startsWith("#track") || line.startsWith("##track")) {
                        trackProperties = new TrackProperties();
                        ParsingUtils.parseTrackLine(line, trackProperties);
                    } else if (line.startsWith("#nodecode") || line.startsWith("##nodecode")) {
                        helper.setUrlDecoding(false);
                    } else if (line.startsWith("#hide") || line.startsWith("##hide")) {
                        String[] kv = line.split("=");
                        if (kv.length > 1) {
                            header.setFeaturesToHide(Arrays.asList(kv[1].split(",")));
                        }
                    }
                    continue;
                } else {
                    break;
                }
            }

            return header;
        }
        catch (IOException e) {
            throw new CodecLineParsingException("Error parsing header", e);
        }
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


        int nTokens = ParsingUtils.split(line, tokens, '\t');

        // GFF files have 9 tokens
        if (nTokens < 9) {
            return null;
        }

        // The type
        String featureType = new String(tokens[2].trim());


        String chromosome = tokens[0];  //genome.getChromosomeAlias(tokens[0]);

        // GFF coordinates are 1-based inclusive (length = end - start + 1)
        // IGV (UCSC) coordinates are 0-based exclusive.  Adjust start and end accordingly
        int start;
        int end;
        try {
            start = Integer.parseInt(tokens[3]) - 1;
        }
        catch (NumberFormatException ne) {
            throw new DataLoadException("Column 4 must contain a numeric value", line);
        }

        try {
            end = Integer.parseInt(tokens[4]);
        }
        catch (NumberFormatException ne) {
            throw new DataLoadException("Column 5 must contain a numeric value", line);
        }

        Strand strand = convertStrand(tokens[6]);

        String attributeString = tokens[8];

        LinkedHashMap<String, String> attributes = new LinkedHashMap();
        //attributes.put("Type", featureType);
        helper.parseAttributes(attributeString, attributes);

        String description = getDescription(attributes, featureType);

        String id = helper.getID(attributes);

        String[] parentIds = helper.getParentIds(attributes, attributeString);

        BasicFeature f = new BasicFeature(chromosome, start, end, strand);
        f.setName(getName(attributes));
        f.setType(featureType);
        f.setDescription(description);
        f.setIdentifier(id);
        f.setParentIds(parentIds);

        if (attributes.containsKey("color")) {
            f.setColor(ColorUtilities.getColorFromString(attributes.get("color")));
        }
        if (attributes.containsKey("Color")) {
            f.setColor(ColorUtilities.getColorFromString(attributes.get("Color")));
        }
        return f;

    }

    public Class getFeatureType() {
        return Feature.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public Object getHeader()  {
        return header;  //To change body of implemented methods use File | Settings | File Templates.
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

        if(attributes == null || attributes.size() == 0) {
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

    protected interface Helper {

        String[] getParentIds(Map<String, String> attributes, String attributeString);

        void parseAttributes(String attributeString, Map<String, String> map);

        String getID(Map<String, String> attributes);

        void setUrlDecoding(boolean b);

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

    class GFF2Helper implements Helper {

        String[] idFields = {"systematic_id", "ID", "Name", "name", "primary_name", "gene", "Locus", "locus", "alias"};


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
}
