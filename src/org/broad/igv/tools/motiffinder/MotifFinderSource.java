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

package org.broad.igv.tools.motiffinder;

import org.broad.igv.dev.api.IGVPlugin;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.CachingFeatureSource;
import org.broad.igv.feature.IGVFeature;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.tribble.Feature;

import javax.swing.*;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class for searching for pattern in a shortSeq.
 *
 * We recognize single letter codes, including ambiguity codes,
 * only.
 * See http://en.wikipedia.org/wiki/Nucleic_acid_notation
 * or http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
 * User: jacob
 * Date: 2013-Jan-22
 */
@XmlAccessorType(XmlAccessType.NONE)
public class MotifFinderSource implements FeatureSource<Feature> {

    @XmlAttribute private String pattern;

    @XmlJavaTypeAdapter(SessionXmlAdapters.Genome.class)
    @XmlAttribute private Genome genome;

    @XmlAttribute private int featureWindowSize = (int) 100e3;

    @SubtlyImportant
    private MotifFinderSource(){}

    /**
     *
     * @param pattern The regex pattern to search
     * @param genome Genome from which to get sequence data
     */
    public MotifFinderSource(String pattern, Genome genome){
        this.pattern = pattern;
        this.genome = genome;
    }

    /**
     * Search the provided sequence for the provided pattern
     * {@code chr} and {@code seqStart} are used in constructing the resulting
     * {@code Feature}s
     * @param chr
     * @param pattern
     * @param posStart The 0-based offset from the beginning of the genome that the {@code sequence} is based
     * @param sequence The nucleotide sequence
     * @return
     */
    public static Iterator<Feature> search(String chr, String pattern, int posStart, byte[] sequence){
        Pattern regex = Pattern.compile(pattern, Pattern.CASE_INSENSITIVE);
        Matcher matcher = regex.matcher(new String(sequence));
        return new MatchFeatureIterator(chr, posStart, matcher);
    }

    @Override
    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {
        byte[] seq = genome.getSequence(chr, start, end);
        if(seq == null) Collections.emptyList().iterator();
        return search(chr, this.pattern, start, seq);
    }

    @Override
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        //TODO Precalculate and/or store?
        return null;
    }

    @Override
    public int getFeatureWindowSize() {
        return this.featureWindowSize;
    }

    @Override
    public void setFeatureWindowSize(int size) {
        this.featureWindowSize = size;
    }

    public static class MotifFinderPlugin implements IGVPlugin {

        /**
         * Add menu entry for activating SequenceMatchDialog
         */
        @Override
        public void init() {
            JMenuItem menuItem = new JMenuItem("Find Motif...");
            menuItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    MotifFinderDialog dialog = new MotifFinderDialog(IGV.getMainFrame());
                    dialog.setVisible(true);

                    String trackName = dialog.getTrackName();
                    String pattern = dialog.getInputPattern();
                    if (pattern != null) {
                        MotifFinderSource source = new MotifFinderSource(pattern, GenomeManager.getInstance().getCurrentGenome());
                        CachingFeatureSource cachingFeatureSource = new CachingFeatureSource(source);
                        FeatureTrack track = new FeatureTrack(trackName, trackName, cachingFeatureSource);
                        IGV.getInstance().addTracks(Arrays.<Track>asList(track), PanelName.FEATURE_PANEL);
                    }
                }
            });

            IGV.getInstance().addOtherToolMenu(menuItem);
        }
    }

    /**
     * Iterator which turns regex Matcher results into Features
     *
     */
    private static class MatchFeatureIterator implements Iterator<Feature>{

        private String chr;
        private int posOffset;

        /**
         * We want to find overlapping matches. By default, matcher.find()
         * starts from the end of the previous match, which would preclude overlapping matches.
         * By storing the last start position we reset each time
         */
        private int lastMatchStart = -1;
        
        private Matcher matcher;

        private IGVFeature nextFeat;

        /**
         * 
         * @param chr The chromosome we are searching
         * @param posOffset The position within the chromosome that we started searching
         * @param matcher Matcher over sequence.
         */
        private MatchFeatureIterator(String chr, int posOffset, Matcher matcher){
            this.chr = chr;
            this.posOffset = posOffset;
            this.matcher = matcher;
            findNext();
        }

        private void findNext(){
            if(matcher.find(lastMatchStart + 1)){
                lastMatchStart = matcher.start();
                int start = posOffset + lastMatchStart;
                int end = posOffset + matcher.end();
                nextFeat = new BasicFeature(chr, start, end);
            }else{
                nextFeat = null;
            }
        }

        @Override
        public boolean hasNext() {
            return nextFeat != null;
        }

        @Override
        public IGVFeature next() {
            IGVFeature nF = nextFeat;
            findNext();
            return nF;
        }

        @Override
        public void remove() {
            throw new RuntimeException("Cannot remove from this iterator");
        }
    }
}
