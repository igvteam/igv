/*
 *  The MIT License (MIT)
 *
 * Copyright (c) 2016 University of California San Diego
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *
 */

package org.broad.igv.sam;

import oracle.jdbc.proxy.annotation.Pre;
import org.broad.igv.PreferenceManager;

import java.util.*;

/**
 * Created by jrobinso on 12/22/16.
 * <p>
 * Experimental class to test strategies for drawing insertions
 */
public class InsertionManager {

    private static InsertionManager theInstance = new InsertionManager();


    Map<Integer, Insertion> insertionMap;
    List<Integer> positions;
    private Integer selectedInsertion;

    public static synchronized InsertionManager getInstance() {
        return theInstance;
    }

    private InsertionManager() {
        this.insertionMap = Collections.synchronizedMap(new HashMap<>(100));
        this.positions = new ArrayList<>(100);
    }

    public void clear() {
        this.insertionMap.clear();
        this.positions.clear();
    }

    public synchronized List<Insertion> getInsertions(double start, double end) {

        List<Insertion> insertions = new ArrayList<>();
        for (int i = 0; i < positions.size(); i++) {
            final Integer key = positions.get(i);
            if (key > end) break;
            if (key >= start) {
                final Insertion insertion = insertionMap.get(key);
                //  if (insertion.size > 2) {
                insertions.add(insertion);
                //  }
            }
        }
        return insertions;

    }

    public void setSelected(int position) {
        this.selectedInsertion = position;
    }

    public void clearSelected() {
        this.selectedInsertion = null;
    }

    public Insertion getSelectedInsertion() {

        return selectedInsertion == null ? null : insertionMap.get(selectedInsertion);

    }

    public Insertion getInsertion(int position) {
        return insertionMap.get(position);
    }

    public void processAlignments(List<Alignment> alignments) {

        int minLength = 0;
        if(PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SAM_HIDE_SMALL_INDEL_BP)) {
            minLength = PreferenceManager.getInstance().getAsInt(PreferenceManager.SAM_SMALL_INDEL_BP_THRESHOLD);
        }

        for (Alignment a : alignments) {
            AlignmentBlock[] blocks = a.getInsertions();
            if (blocks != null) {
                for (AlignmentBlock block : blocks) {

                    if(block.getBases().length < minLength) continue;

                    Integer key = block.getStart();
                    Insertion insertion = insertionMap.get(key);
                    if (insertion == null) {
                        insertion = new Insertion(block.getStart(), block.getLength());

                        if(insertion.size > 1000) {
                            System.out.println();
                        }

                        insertionMap.put(key, insertion);
                        positions.add(block.getStart());
                    } else {
                        insertion.size = Math.max(insertion.size, block.getLength());
                    }
                }
            }
        }


        this.positions = new ArrayList<>(insertionMap.keySet());
        this.positions.sort((o1, o2) -> o1 - o2);
    }

    public static class Insertion {
        public int position;
        public int size;
        public int pixelPosition = -1;

        public Insertion(int position, int size) {
            this.position = position;
            this.size = size;
        }
    }


}
