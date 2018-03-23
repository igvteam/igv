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

package org.broad.igv.sam.reader;

import org.broad.igv.Globals;
import org.broad.igv.ui.util.IndexCreatorDialog;
import org.broad.igv.ui.util.UIUtilities;
import htsjdk.tribble.readers.AsciiLineReader;

import javax.swing.*;
import java.io.*;

/**
 * @author jrobinso
 */
public abstract class AlignmentIndexer {

    static int DEFAULT_TILEWIDTH = 16000;
    static int MILLISECONDS_IN_MINUTE = 60 * 1000;
    /**
     * Optional progress bar.
     * TODO -- this should be done through event/listeners
     */
    IndexCreatorDialog.IndexWorker worker;
    JProgressBar progressBar;
    File samFile;

    public static AlignmentIndexer getInstance(File samFile,
                                               JProgressBar progressBar,
                                               IndexCreatorDialog.SamIndexWorker worker) {

        String fn = samFile.getName().toLowerCase();
        if (fn.endsWith(".aligned") || fn.endsWith(".bedz") || fn.endsWith(".bed") ||
                fn.endsWith(".aligned.txt") || fn.endsWith(".bedz.txt") || fn.endsWith(".bed.txt")) {

            return new DotAlignedIndexer(samFile, progressBar, worker);
        }
        return new SamIndexer(samFile, progressBar, worker);
    }

    AlignmentIndexer(File samFile,
                     JProgressBar progressBar,
                     IndexCreatorDialog.SamIndexWorker worker) {
        this.samFile = samFile;
        this.progressBar = progressBar;
        this.worker = worker;
    }

    public FeatureIndex createSamIndex() throws
            IOException,
            FileNotFoundException {
        File idxFile = new File(samFile.getParent(), samFile.getName() + ".sai");
        return createSamIndex(idxFile);
    }

    public FeatureIndex createSamIndex(File idxFile) throws
            IOException {
        return createSamIndex(idxFile, DEFAULT_TILEWIDTH);
    }

    public FeatureIndex createSamIndex(
            File idxFile, int tileWidth) throws
            IOException {

        FileInputStream is = new FileInputStream(samFile);
        InputStream bis = new BufferedInputStream(is);
        AsciiLineReader reader = new AsciiLineReader(bis);

        long fileLength = samFile.length();
        long progressIncrement = fileLength / 100;

        long lastFilePosition = 0;


        String lastChr = null;
        int lastAlignmentStart = 0;

        FeatureIndex featureIndex = new FeatureIndex(tileWidth);
        int recordCount = 0;
        long filePosition = 0;
        int currentTile = 0;
        int longestFeature = 0;

        long startTime = System.currentTimeMillis();
        int progressCounter = 1; // progress in %

        String nextLine = "";
        int lineNumber = 0;
        while ((nextLine = reader.readLine()) != null) {

            lineNumber++;

            if (worker != null && worker.isCancelled()) {
                return null;
            }

            //int nBytes = nextLine.length();
            nextLine = nextLine.trim();
            String[] fields = Globals.tabPattern.split(nextLine, -1);
            int nFields = fields.length;
            if (!nextLine.startsWith("@") && nFields > 3 && isMapped(fields)) {

                String chr = getChromosome(fields);
                int alignmentStart = getAlignmentStart(fields);
                int tileNumber = alignmentStart / tileWidth;

                if (lastChr == null) {
                    // First record
                    currentTile = tileNumber;
                    for (int i = 0; i < currentTile; i++) {
                        featureIndex.add(chr, lastFilePosition, 0, longestFeature);
                    }
                    lastChr = chr;

                } else if (!chr.equals(lastChr)) {   // New chromosome
                    featureIndex.add(lastChr, filePosition, recordCount, longestFeature);
                    filePosition = lastFilePosition;

                    currentTile = 0;
                    recordCount = 0;
                    lastAlignmentStart = 0;
                    longestFeature = 0;
                    lastChr = chr;
                } else {

                    longestFeature = Math.max(longestFeature, getAlignmentLength(fields));

                    if (alignmentStart < 0) {
                        System.out.println("Warning: negative start position at line: " + lineNumber + " : " + nextLine);
                        continue;
                    }

                    if (alignmentStart < lastAlignmentStart) {
                        throw new UnsortedFileException(" File must be sorted by start position. " +
                                "Sort test failed at: " + nextLine);
                    }

                    lastAlignmentStart = alignmentStart;

                    if (tileNumber > currentTile) {

                        // We have crossed a tile boundary.  Record index and counts for previous tile
                        featureIndex.add(lastChr, filePosition, recordCount, longestFeature);

                        // If tiles were skipped record zero counts for these.
                        for (int cnt = 0; cnt < (tileNumber - currentTile - 1); cnt++) {
                            featureIndex.add(lastChr, filePosition, 0, longestFeature);
                        }

                        filePosition = lastFilePosition;
                        currentTile = tileNumber;
                        recordCount = 0;
                    }
                    recordCount++;
                }
            }

            lastFilePosition = reader.getPosition();

            if (lastFilePosition > (progressCounter * progressIncrement)) {
                updateProgress(progressCounter, startTime);
                progressCounter++;

            }
        }

        // Record last partial tile
        featureIndex.add(lastChr, filePosition, recordCount, longestFeature);

        is.close();

        if (idxFile != null) {
            featureIndex.store(idxFile);
        }
        //Done now
        updateProgress(100, startTime);
        if (progressBar == null) {
            System.out.println("Done indexing " + samFile.getName());
        }

        return featureIndex;

    }

    abstract int getAlignmentStart(String[] fields) throws NumberFormatException;

    abstract int getAlignmentLength(String[] fields) throws NumberFormatException;

    abstract String getChromosome(String[] fields);

    abstract boolean isMapped(String[] fields);


    private void updateProgress(int progressCounter, long startTime) {
        final long timeToComplete = ((100 - progressCounter) *
                (System.currentTimeMillis() - startTime)) / progressCounter;
        final int p = progressCounter;
        if (progressBar != null) {
            UIUtilities.invokeOnEventThread(new Runnable() {

                public void run() {
                    progressBar.setValue(p);
                    if (worker != null) {
                        worker.setTimeRemaining(timeToComplete);
                    }

                }
            });
        } else {
            System.out.println(
                    "Progress: " + progressCounter + "%");
        }


    }

}

