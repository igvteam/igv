/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.methyl;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;

import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Data source backe by Michael Ziller's custom BB format
 *
 * @author Jim Robinson
 * @date 4/19/12
 */
public class BBMethylDataSource implements MethylDataSource {

    static Pattern percentPattern = Pattern.compile("%");

    enum Type {ZILLER, USC};

    BBFileReader reader;
    Type type;

    /**
     * Map of chr name in genome definition -> chr name in BB file.  Used for reverse-lookup during queries.
     */
    Map<String, String> chrNameMap;

    public BBMethylDataSource(BBFileReader reader, Type type, Genome genome) throws IOException {
        this.reader = reader;
        this.type = type;
        init(genome);
    }

    public Iterator<MethylScore> query(String chr, int start, int end){
        String tmp = chrNameMap.get(chr);
        String querySeq = tmp == null ? chr : tmp;
        BigBedIterator bedIterator = reader.getBigBedIterator(querySeq, start, chr, end, false);
        return new WrappedIterator(bedIterator, type);
    }

    private void init(Genome genome) {
        chrNameMap = new HashMap<String, String>();
        if (genome != null) {
            Collection<String> seqNames = reader.getChromosomeNames();
            if (seqNames != null)
                for (String seqName : seqNames) {
                    String igvChr = genome.getChromosomeAlias(seqName);
                    if (igvChr != null && !igvChr.equals(seqName)) {
                        chrNameMap.put(igvChr, seqName);
                    }
                }
        }
    }


    public static class WrappedIterator implements Iterator<MethylScore> {

        BigBedIterator bedIterator;
        Type type;

        public WrappedIterator(BigBedIterator bedIterator, Type type) {
            this.bedIterator = bedIterator;
            this.type = type;
        }

        public boolean hasNext() {
            return bedIterator.hasNext();  //To change body of implemented methods use File | Settings | File Templates.
        }

        public MethylScore next() {

            BedFeature feat = null;

            while (feat == null && bedIterator.hasNext()) {
                feat = bedIterator.next();
                String[] restOfFields = feat.getRestOfFields();
                MethylScore score = type == Type.ZILLER ?
                        createZillerScore(feat, restOfFields) :
                        createUSCScore(feat, restOfFields);


                return score;
            }

            return null;

        }

        private MethylScore createZillerScore(BedFeature feat, String[] restOfFields) {
            String name = restOfFields[0];
            //'92%[51]'
            String[] tokens = percentPattern.split(name.replace("'", "").replace("[", "").replace("]", ""));
            float percent = Float.parseFloat(tokens[0]);
            int count = Integer.parseInt(tokens[1]);
            return new MethylScore(feat.getChromosome(), feat.getStartBase(), feat.getEndBase(),
                    Strand.NONE, percent, count);
        }

        private MethylScore createUSCScore(BedFeature feat, String[] restOfFields) {
            String name = restOfFields[0];
            int score = Integer.parseInt(restOfFields[1]);
            char strandChar = restOfFields[2].charAt(0);
            Strand strand = strandChar == '+'  ? Strand.POSITIVE :
                    (strandChar == '-' ? Strand.NEGATIVE : Strand.NONE);
            float percentMethyl = Float.parseFloat(restOfFields[3]);
            int count = Integer.parseInt(restOfFields[4]);
            return new MethylScore(feat.getChromosome(), feat.getStartBase(), feat.getEndBase(), strand, percentMethyl, count);
         }

        public void remove() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
