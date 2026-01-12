package org.igv.sam;

import org.igv.feature.genome.Genome;
import org.igv.feature.genome.GenomeManager;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;

import java.util.*;
import java.util.stream.Collectors;

import static org.igv.prefs.Constants.SAM_HIDE_SMALL_INDEL;
import static org.igv.prefs.Constants.SAM_SMALL_INDEL_BP_THRESHOLD;

/**
 * Created by jrobinso on 12/22/16.
 * <p>
 * Experimental class to test strategies for drawing insertions
 */
public class InsertionManager {

    private static InsertionManager theInstance = new InsertionManager();

    private Map<String, Map<Integer, InsertionMarker>> insertionMaps;
    private Map<String, Integer> selectedInsertions;

    InsertionMarker selectedInsertion;

    public static InsertionManager getInstance() {
        return theInstance;
    }

    private InsertionManager() {
        this.insertionMaps = Collections.synchronizedMap(new HashMap<>(100));
        this.selectedInsertions = Collections.synchronizedMap(new HashMap<>(100));
    }

    public void clear() {
        this.insertionMaps.clear();
        this.selectedInsertions.clear();;
        this.selectedInsertion = null;
    }

    public List<InsertionMarker> getInsertions(String chrName, double start, double end) {

        Map<Integer, InsertionMarker> insertionMap = insertionMaps.get(chrName);
        if(insertionMap == null ) return null;

        return insertionMap.values().stream().filter(im -> im.position + im.size >= start && im.position <= end).collect(Collectors.toList());

    }

    public void setSelected(InsertionMarker insertionMarker) {
        this.selectedInsertion = insertionMarker;
    }

    public synchronized void processAlignments(String chr, List<Alignment> alignments) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        chr = genome == null ? chr : genome.getCanonicalChrName(chr);

        Map<Integer, InsertionMarker> insertionMap = insertionMaps.get(chr);
        if(insertionMap == null) {
            insertionMap =  Collections.synchronizedMap(new HashMap<>());
            insertionMaps.put(chr, insertionMap);
        }

        for (Alignment a : alignments) {
            AlignmentBlock[] blocks = a.getInsertions();
            if (blocks != null) {
                for (AlignmentBlock block : blocks) {
                    Integer key = block.getStart();
                    InsertionMarker insertionMarker = insertionMap.get(key);
                    if (insertionMarker == null) {
                        insertionMarker = new InsertionMarker(block.getStart(), block.getLength());
                        insertionMap.put(key, insertionMarker);
                    } else {
                        insertionMarker.size = Math.max(insertionMarker.size, block.getLength());
                    }
                }
            }
        }
    }

}
