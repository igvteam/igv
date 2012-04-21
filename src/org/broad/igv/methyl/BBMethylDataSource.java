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
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.tribble.Feature;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.awt.geom.FlatteningPathIterator;
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

    BBFileReader reader;

    /**
     * Map of chr name in genome definition -> chr name in BB file.  Used for reverse-lookup during queries.
     */
    Map<String, String> chrNameMap;

    public BBMethylDataSource(String path, Genome genome) throws IOException {
        reader = new BBFileReader(path);
        init(genome);
    }

    public Iterator<MethylScore> query(String chr, int start, int end) {
        String tmp = chrNameMap.get(chr);
        String querySeq = tmp == null ? chr : tmp;
        BigBedIterator bedIterator = reader.getBigBedIterator(querySeq, start, chr, end, false);
        return new WrappedIterator(bedIterator);
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

        public WrappedIterator(BigBedIterator bedIterator) {
            this.bedIterator = bedIterator;
        }

        public boolean hasNext() {
            return bedIterator.hasNext();  //To change body of implemented methods use File | Settings | File Templates.
        }

        public MethylScore next() {

            BedFeature feat = bedIterator.next();
            BasicFeature feature = new BasicFeature(feat.getChromosome(), feat.getStartBase(), feat.getEndBase());
            String[] restOfFields = feat.getRestOfFields();
            String name = restOfFields[0];
            //‘92%[51]’

            String[] tokens = percentPattern.split(name.replace("'", "").replace("[", "").replace("]", ""));
            if (tokens.length != 2) {
                // What to do, throw exception?
            }
            short percent = Short.parseShort(tokens[0]);
            short count;
            try {
                count = Short.parseShort(tokens[1]);
            } catch (NumberFormatException e) {
                count = Short.MAX_VALUE;    // Protect against the occassional massive coverage depth
            }
            return new MethylScore(feat.getChromosome(), feat.getStartBase(), feat.getEndBase(), percent, count);

        }

        public void remove() {
            //To change body of implemented methods use File | Settings | File Templates.
        }
    }
}
