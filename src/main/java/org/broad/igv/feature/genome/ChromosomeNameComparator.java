package org.broad.igv.feature.genome;

import java.util.Comparator;

/**
 * Comparator for chromosome names.
 * -- Previously complex logic is deprecated.  Now just does a simple string compare
 */
public class ChromosomeNameComparator implements Comparator<String> {

    private static ChromosomeNameComparator instance;

    private ChromosomeNameComparator() {
    }

    public static ChromosomeNameComparator get() {
        if (instance == null) {
            instance = new ChromosomeNameComparator();
        }
        return instance;
    }

    public int compare(String chr1, String chr2) {
        return chr1.compareTo(chr2);
    }
}
