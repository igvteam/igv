package org.broad.igv.dev.affective;

import org.broad.igv.Globals;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFGroup;
import org.broad.igv.tdf.TDFReader;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
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

        AffectiveChromosome chromosome = (AffectiveChromosome) affectiveGenome.getChromosome(chr);
        int chrStartTime = chromosome.startTime;
        int chrSamplingRate = chromosome.samplingRate;

        int dataStartTime = startTimeMap.get(chr);
        int dataSamplingRate = samplingRateMap.get(chr);

        // TODO -- if chrSamplingRate != dataSampingRate some adjustments need to be made.

        int offset = (dataStartTime - chrStartTime) * dataSamplingRate;

        List<LocusScore> scores = super.getSummaryScores(chr, startLocation, endLocation, zoom);
        for(LocusScore score : scores) {
            score.setStart(score.getStart() + offset);
            score.setEnd(score.getEnd()  + offset);
        }
        return scores;
    }
}
