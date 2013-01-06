package org.broad.igv.dev.affective;

import org.broad.igv.Globals;
import org.broad.igv.data.BasicScore;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.*;

import java.util.*;
import java.util.regex.Pattern;

/**
 * @author Jim Robinson
 * @date 1/23/12
 */
public class AffectiveDataSource extends TDFDataSource {

    static Pattern colonPattern = Pattern.compile(":");
    Map<String, Integer> startTimeMap = new HashMap<String, Integer>();
    Map<String, Integer> samplingRateMap = new HashMap<String, Integer>();
    AffectiveGenome affectiveGenome;

    // TODO -- invalidate when new data is loaded (genome changes)
    List<LocusScore> wholeGenomeScores = null;

    public AffectiveDataSource(TDFReader reader, int trackNumber, String trackName, Genome genome) {

        super(reader, trackNumber, trackName, genome);

        this.affectiveGenome = (AffectiveGenome) genome;

        // "dynamic" genome that is defined by the data files
        if (genome instanceof AffectiveGenome) {
            for (String chrName : reader.getChromosomeNames()) {
                if (!chrName.equals(Globals.CHR_ALL)) {
                    ((AffectiveGenome) genome).createChromosome(chrName);
                }
            }
        }

        // Get the start times for each date
        TDFGroup rootGroup = reader.getGroup("/");
        String prefix = "ATTR:";
        for (String attrKey : rootGroup.getAttributeNames()) {
            if (attrKey.startsWith(prefix) && attrKey.endsWith("startTime")) {
                String[] tokens = colonPattern.split(attrKey);
                if (tokens.length == 3) {
                    String chr = tokens[1];
                    String value = rootGroup.getAttribute(attrKey);
                    Integer startTime = new Integer(value);     // Start time in seconds from midnight
                    startTimeMap.put(chr, startTime);
                }
            } else if (attrKey.startsWith(prefix) && attrKey.endsWith("samplingRate")) {
                String[] tokens = colonPattern.split(attrKey);
                if (tokens.length == 3) {
                    String chr = tokens[1];
                    String value = rootGroup.getAttribute(attrKey);
                    Integer samplingRate = new Integer(value);     // Start time in seconds from midnight
                    samplingRateMap.put(chr, samplingRate);
                }
            }
        }
    }

    // This is a copy of the super method, wih an offset
    @Override
    protected List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {

        if (chr.equals(Globals.CHR_ALL)) {
            // todo
            return getWholeGenomeScores();
        } else {
            return super.getSummaryScores(chr, startLocation, endLocation, zoom);
        }
    }

    public List<LocusScore> getWholeGenomeScores() {

        if (wholeGenomeScores == null && affectiveGenome.getLongChromosomeNames().size() > 1) {
            wholeGenomeScores = new ArrayList<LocusScore>(1000);
            for (String chr : affectiveGenome.getLongChromosomeNames()) {
                wholeGenomeScores.addAll(getWholeGenomeScoresForChromosome(chr));
            }

        }
        return wholeGenomeScores;
    }

    List<LocusScore> getWholeGenomeScoresForChromosome(String chr) {

        long offset = affectiveGenome.getCumulativeOffset(chr);
        ArrayList<LocusScore> scores = new ArrayList<LocusScore>();
        List<LocusScore> tmp = getSummaryScores(chr, 0, Integer.MAX_VALUE, 0);
        if (tmp != null) {
            float value = 0;
            int lastWGStart = (int) ((tmp.get(0).getStart() + offset) / 1000);
            int lastWGEnd = (int) ((tmp.get(0).getEnd() + offset) / 1000);
            int numPoints = 0;
            for (LocusScore s : tmp) {
                int wgStart = (int) ((s.getStart() + offset) / 1000);
                int wgEnd = (int) ((s.getEnd() + offset) / 1000);
                if (Float.isNaN(s.getScore())) {

                }
                if (wgEnd > lastWGEnd) {
                    scores.add(new BasicScore(lastWGStart, lastWGEnd, value / numPoints));
                    lastWGStart = wgStart;
                    lastWGEnd = wgEnd;
                    value = 0;
                    numPoints = 0;
                }
                value += s.getScore();
                numPoints++;
            }
        }
        return scores;
    }



}
