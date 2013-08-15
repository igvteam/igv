/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.gwas;

import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.FeatureFileHeader;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.StringUtils;
import org.broad.tribble.AsciiFeatureCodec;
import org.broad.tribble.Feature;
import org.broad.tribble.readers.LineIterator;

import java.util.HashMap;
import java.util.Map;

/**
 * CODEC for GTex project eQTL files
 * <p/>
 * SNP	SNP_Chr	SNP_Pos	Gen_ID	Gene_Name	Gene_Pos	T_Stat	P_Val	Q_Val
 * rs1569471	1	169564130	ENSG00000000460.11	C1orf112	169631245	-4.187361378	7.79E-05	0.04794564
 */
public class EQTLCodec extends AsciiFeatureCodec<EQTLFeature> {

    String[] columnNames;
    Genome genome;
    private FeatureFileHeader header;

    public EQTLCodec(Genome genome) {
        super(EQTLFeature.class);
        this.genome = genome;
    }

    //@Override
    public Feature decodeLoc(String line) {
        String[] tokens = Globals.tabPattern.split(line);
        if (tokens[0].equals("SNP")) return null;
        String chr = tokens[1];
        int position = Integer.parseInt(tokens[2]) - 1;
        return new BasicFeature(chr, position, position + 1);
    }


    @Override
    public EQTLFeature decode(String s) {

        String[] tokens = Globals.tabPattern.split(s);
        if (tokens[0].equals("SNP")) {
            // This is the header
            columnNames = tokens;
            return null;
        }

        String snp = tokens[0];
        String chr = genome == null ? StringUtils.intern(tokens[1]) : genome.getChromosomeAlias(tokens[1]);

        int position = Integer.parseInt(tokens[2]) - 1;
        String geneId = tokens[3];
        String geneName = tokens[4];

        Map<String, String> attributes = null;
        if (columnNames != null) {
            attributes = new HashMap<String, String>();
            for (int i = 5; i < tokens.length; i++) {
                if (columnNames.length < i) {
                    attributes.put(columnNames[i], tokens[i]);
                }
            }
        }


        return new EQTLFeature(snp, chr, position, geneId, geneName, attributes);
    }


    /**
     * EQTL files do not have interesting headers, we are using (abusing) this method to set some initial
     * properties.  The fact we have to do this => the design needs more work.
     *
     * @param reader
     * @return
     */
    @Override
    public Object readActualHeader(LineIterator reader) {

        if (header == null) {
            header = new FeatureFileHeader();
            TrackProperties tp = new TrackProperties();
            tp.setDisplayMode(Track.DisplayMode.EXPANDED);
            header.setTrackProperties(tp);
        }
        return header;

    }
}
