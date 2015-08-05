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

package org.broad.igv.data;

import com.google.common.collect.Iterators;
import org.apache.log4j.Logger;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.DataTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.util.*;

/**
 * Data source which combines two other DataSources
 *
 * TODO Multiple DataSources. There is no 2
 * @author jrobinso, jacob
 */
@XmlAccessorType(XmlAccessType.NONE)
public class CombinedDataSource implements DataSource {

    private static Logger log = Logger.getLogger(CombinedDataSource.class);

    public enum Operation{
        ADD("+"),
        SUBTRACT("-"),
        MULTIPLY("*"),
        DIVIDE("/");

        private String stringRep;

        private Operation(String stringrep){
            this.stringRep = stringrep;
        }

    }

    @XmlJavaTypeAdapter(SessionXmlAdapters.DataTrackIDAdapter.class)
    @XmlAttribute(name = "source0")
    DataTrack source0;

    @XmlJavaTypeAdapter(SessionXmlAdapters.DataTrackIDAdapter.class)
    @XmlAttribute(name = "source1")
    DataTrack source1;

    @XmlAttribute
    Operation operation = Operation.ADD;

    @SubtlyImportant
    private CombinedDataSource(){}

    public CombinedDataSource(DataTrack source0, DataTrack source1, Operation operation){
        this.source0 = source0;
        this.source1 = source1;
        this.operation = operation;
    }

    public void updateTrackReferences(List<Track> allTracks) {
        //We filled in sources with placeholder tracks if not found, now find the real ones
        source0 = updateTrackReference(source0, allTracks);
        source1 = updateTrackReference(source1, allTracks);
    }

    private DataTrack updateTrackReference(DataTrack memberTrack, List<Track> allTracks){

        if(memberTrack.getName() == null && memberTrack.getResourceLocator() == null){
            DataTrack matchingTrack = (DataTrack) IGVSessionReader.getMatchingTrack(memberTrack.getId(), allTracks);
            if(matchingTrack == null) throw new IllegalStateException("Could not find track with ID " + memberTrack.getId());
            return matchingTrack;
        }else{
            return memberTrack;
        }

    }

    public List<LocusScore> getSummaryScoresForRange(String chr, int startLocation, int endLocation, int zoom){

        List<LocusScore> outerScores = this.source0.getSummaryScores(chr, startLocation, endLocation, zoom);
        List<LocusScore> innerScores = this.source1.getSummaryScores(chr, startLocation, endLocation, zoom);

        int initialSize = outerScores.size() + innerScores.size();
        List<LocusScore> combinedScoresList = new ArrayList<LocusScore>(initialSize);

        if(initialSize == 0) return combinedScoresList;

        //TODO We assume that having no data from one source is the identity operation, that may not be true
        if(innerScores.size() == 0) return outerScores;
        if(outerScores.size() == 0) return innerScores;


        /**
         * We first generate the chunks which will need to be calculated separately
         * This is the set of all start/end positions of outerScores and innerScores
         * We could be a bit smarter, but this is simpler and there's no problem with
         * skipping over intervals which don't have data later.
         *
         * Following that, for each interval generated, we search outerScores and innerScores
         * for the unique LocusScore which contains the generated interval.
         */

        //Generate the boundaries for the new combined regions
        Set<Integer> boundariesSet = new LinkedHashSet<Integer>(2*initialSize);
        Iterator<LocusScore> dualIter = Iterators.mergeSorted(Arrays.asList(innerScores.iterator(), outerScores.iterator()),
                new Comparator<LocusScore>() {

                    @Override
                    public int compare(LocusScore o1, LocusScore o2) {
                        return o1.getStart() - o2.getStart();
                    }
                });
        while(dualIter.hasNext()){
            LocusScore score = dualIter.next();
            boundariesSet.add(score.getStart());
            boundariesSet.add(score.getEnd());
        }
        Integer[] boundariesArray = boundariesSet.toArray(new Integer[0]);
        Arrays.sort(boundariesArray);

        int outerScoreInd = 0;
        int innerScoreInd = 0;
        //Calculate value for each interval
        for(int bb=0; bb < boundariesArray.length-1; bb++){
            int start = boundariesArray[bb];
            int end = boundariesArray[bb+1];
            //It shouldn't be possible for more than one LocusScore of either
            //tracks to overlap each interval, since the start/ends
            //were based on all start/ends of the inputs
            outerScoreInd = findContains(start, end, outerScores, Math.max(outerScoreInd, 0));
            innerScoreInd = findContains(start, end, innerScores, Math.max(innerScoreInd, 0));
            LocusScore outerScore = getContains(outerScores, outerScoreInd);
            LocusScore innerScore = getContains(innerScores, innerScoreInd);

            if(outerScore == null && innerScore == null) continue;
            float score = combineScores(outerScore,  innerScore);
            BasicScore newScore = new BasicScore(start, end, score);
            combinedScoresList.add(newScore);
        }
        return combinedScoresList;
    }

    /**
     * Search {@code scoresList} (must be sorted by start position) for a score which contains the interval specified
     * by start/end. The first one which satisfies this requirement is returned.
     *
     * @param start
     * @param end
     * @param scoresList
     * @param startIndex Optimization, where to start searching in {@code scoresList}
     **/
    private int findContains(int start, int end, List<LocusScore> scoresList, int startIndex) {
        for(int ii=startIndex; ii < scoresList.size(); ii++){
            LocusScore score = scoresList.get(ii);
            if(score.getStart() <= start && score.getEnd() >= end){
                return ii;
            }else if(score.getStart() >= end){
                return -1;
            }
        }
        return -1;
    }

    private LocusScore getContains(List<LocusScore> scores, int index){
        if(index >= 0 && index < scores.size()){
            return scores.get(index);
        }else{
            return null;
        }
    }

    /**
     * Combine the scores using this sources operation. Either can be null,
     * in which case it is given the coordinates of the other and a score of 0.
     * Both inputs cannot be null
     * @param score0
     * @param score1
     * @return
     */
    private float combineScores(LocusScore score0, LocusScore score1) {
        if(score0 == null && score1 == null) throw new IllegalArgumentException("Both inputs cannot be null");
        if(score0 == null){
            score0 = new BasicScore(score1.getStart(), score1.getEnd(), 0.0f);
        }else if(score1 == null){
            score1 = new BasicScore(score0.getStart(), score0.getEnd(), 0.0f);
        }

        switch(operation){
            case ADD:
                return score0.getScore() + score1.getScore();
            case SUBTRACT:
                return score0.getScore() - score1.getScore();
            case MULTIPLY:
                return score0.getScore() * score1.getScore();
            case DIVIDE:
                if(score1.getScore() == 0.0f){
                    return 0.0f;
                }
                return score0.getScore() / score1.getScore();
            default:
                throw new IllegalStateException("Operation not recognized: " + operation);
        }
    }

    public double getDataMax() {
        return 0;
    }

    public double getDataMin() {
        return 0;
    }

    public TrackType getTrackType() {
        return TrackType.PLUGIN;
    }

    public void setWindowFunction(WindowFunction statType) {
        //TODO
    }

    public boolean isLogNormalized() {
        return false;
    }

    public WindowFunction getWindowFunction() {
        return WindowFunction.none;
    }

    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return new ArrayList<WindowFunction>();
    }

    @Override
    public void dispose() {
        if(source0 != null) source0.dispose();
        if(source1 != null) source1.dispose();
    }

}
