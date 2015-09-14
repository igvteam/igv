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

package org.broad.igv.methyl;

import org.broad.igv.Globals;
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

    enum Type {ZILLER, USC}

    ;

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

    public Iterator<MethylScore> query(String chr, int start, int end) {
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
                        createZillerScore(feat, restOfFields) : createUSCScore(feat, restOfFields);
                return score;
            }
            return null;

        }

        private MethylScore createZillerScore(BedFeature feat, String[] restOfFields) {

            float percent;
            float count;
            String name = restOfFields[0];
            if (name.contains("%")) {
                //'92%[51]'
                String[] tokens = percentPattern.split(name.replace("'", "").replace("[", "").replace("]", ""));
                percent = Float.parseFloat(tokens[0]);
                count = Float.parseFloat(tokens[1]);
            } else {
                //  methylatedReads/totalreads
                String[] tokens = Globals.forwardSlashPattern.split(name.replace("'", ""));
                float methylatedReads = Float.parseFloat(tokens[0]);
                count = Float.parseFloat(tokens[1]);
                percent = (methylatedReads / count) * 100;
            }
            return new MethylScore(feat.getChromosome(), feat.getStartBase(), feat.getEndBase(), Strand.NONE, percent, (int) count);
        }

        private MethylScore createUSCScore(BedFeature feat, String[] restOfFields) {
            //String name = restOfFields[0];
            //int score = Integer.parseInt(restOfFields[1]);
            char strandChar = restOfFields[2].charAt(0);
            Strand strand = strandChar == '+' ? Strand.POSITIVE : (strandChar == '-' ? Strand.NEGATIVE : Strand.NONE);
            float percentMethyl = Float.parseFloat(restOfFields[3]);
            int count = Integer.parseInt(restOfFields[4]);
            return new MethylScore(feat.getChromosome(), feat.getStartBase(), feat.getEndBase(), strand, percentMethyl, count);
        }

        public void remove() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
