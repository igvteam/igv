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

package org.broad.igv.util;

import org.broad.igv.feature.Locus;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureSource;
import org.broad.tribble.Feature;

import java.io.BufferedReader;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class for searching for motif in a shortSeq.
 *
 * We recognize single letter codes, including ambiguity codes,
 * only.
 * See http://en.wikipedia.org/wiki/Nucleic_acid_notation
 * or http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
 * User: jacob
 * Date: 2013-Jan-22
 */
public class SequenceMatchSource implements FeatureSource<Feature>{

    private static Map<String, String> letterToRegex;

    static{
        initLetterToRegex();
    }

    private static final String codeFilePath = "resources/iupac_regex_table.txt";
    private static void initLetterToRegex() {
        URL url = SequenceMatchSource.class.getResource(codeFilePath);
        letterToRegex = loadMap(StringUtils.decodeURL(url.getPath()));
    }

    private String motif;
    private Genome genome;

    private int featureWindowSize = (int) 100e3;

    /**
     *
     * @param motif The string to search for, which can include IUPAC ambiguity codes
     * @param genome Genome from which to get sequence data
     */
    public SequenceMatchSource(String motif, Genome genome){
        this.motif = motif;
        this.genome = genome;
    }
    /**
     * Replace the ambiguity codes in the motif
     * with regular expression equivalents
     * @param motif
     * @return
     */
    static String convertMotifToRegex(String motif){
        String output = motif;
        int outloc = 0;
        for(int inloc=0; inloc < motif.length(); inloc++){

            String inchar = motif.substring(inloc, inloc + 1);
            String rep = letterToRegex.get(inchar);

            output = output.substring(0, outloc) + rep + motif.substring(inloc + 1);
            outloc += rep.length();
        }
        return output;
    }

    /**
     * Search the provided sequence for the provided motif
     * {@code chr} and {@code seqStart} are used in constructing the resulting
     * {@code Feature}s
     * @param chr
     * @param motif
     * @param posStart The 0-based offset from the beginning of the genome that the {@code sequence} is based
     * @param sequence
     * @return
     */
    public static Iterator<Feature> search(String chr, String motif, int posStart, byte[] sequence){
        motif = convertMotifToRegex(motif);
        Pattern regex = Pattern.compile(motif, Pattern.CASE_INSENSITIVE);
        Matcher matcher = regex.matcher(new String(sequence));
        return new MatchFeatureIterator(chr, posStart, matcher);
    }

    /**
     * TODO Move this to someplace more general, use it wherever we store lots of this kind of data
     * @param path
     * @return
     */
    public static Map<String, String> loadMap(String path){
        BufferedReader reader = null;
        Map<String, String> map = new HashMap<String, String>();
        try {
            reader = ParsingUtils.openBufferedReader(path);
            String nextLine = null;
            while ((nextLine = reader.readLine()) != null) {
                if(nextLine.startsWith("#")) continue;

                String[] tokens = nextLine.split("=");
                if (tokens.length == 2) {
                    map.put(tokens[0], tokens[1]);
                }else{
                    throw new IllegalArgumentException("Incorrect number of tokens at line: " + nextLine);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

        return map;
    }

    @Override
    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {
        return search(chr, this.motif, start, genome.getSequence(chr, start, end));
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

    /**
     * Iterator which turns regex Matcher results into Features
     * Implementation note: We don't want to have the constructor run the query,
     *
     */
    private static class MatchFeatureIterator implements Iterator<Feature>{

        private String chr;
        private int posOffset;
        
        private Matcher matcher;

        private Locus nextFeat;

        /**
         * 
         * @param chr The chromosome we are searching
         * @param posOffset The position within the chromosome that we started searching
         * @param matcher Matcher over sequence. Will be reset
         */
        private MatchFeatureIterator(String chr, int posOffset, Matcher matcher){
            this.chr = chr;
            this.posOffset = posOffset;
            this.matcher = matcher;
            this.matcher.reset();
            findNext();
        }

        private void findNext(){
            if(matcher.find()){
                int start = posOffset + matcher.start();
                int end = posOffset + matcher.end();
                nextFeat = new Locus(chr, start, end);
            }else{
                nextFeat = null;
            }
        }

        @Override
        public boolean hasNext() {
            return nextFeat != null;
        }

        @Override
        public Locus next() {
            Locus nF = nextFeat;
            findNext();
            return nF;
        }

        @Override
        public void remove() {
            throw new RuntimeException("Cannot remove from this iterator");
        }
    }
}
