package org.broad.igv.hic.tools;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * @author nchernia
 * Date: 8/26/12
 */
public class FragmentCalculation {

    private Map<String, int[]> fragmentMap = null;

    public static FragmentCalculation readFragments(InputStream is) throws IOException {
        Pattern pattern = Pattern.compile("\\s");
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        String nextLine;
        Map<String, int[]> fragmentMap = new LinkedHashMap<String, int[]>();

        while ((nextLine = reader.readLine()) != null) {
            String[] tokens = pattern.split(nextLine);
            if (tokens.length > 1) {
                String key = tokens[0];
                int[] sites = new int[tokens.length - 1];
                for (int i = 1; i < tokens.length; i++) {
                    sites[i - 1] = Integer.parseInt(tokens[i]);
                }

                fragmentMap.put(key, sites);
            } else {
                System.out.println("Skipping line: " + nextLine);
            }
        }

        return new FragmentCalculation(fragmentMap);
    }

    public static FragmentCalculation readFragments(String filename) throws IOException {
        InputStream is = null;
        try {
            File file = new File(filename);
            is = new FileInputStream(file);
            return readFragments(is);
        } finally {
            is.close();
        }

    }

    public FragmentCalculation(Map<String, int[]> fragmentMap) {
        this.fragmentMap = fragmentMap;
    }

    public int [] getSites(String chrName) {
        return fragmentMap.get(chrName);
    }

    public int getNumberFragments(String chrName) {
        int[] sites = fragmentMap.get(chrName);
        if (sites == null) // for "All"
            return 0;
        return sites.length;
    }

    public int getBin(String chrName, int position) {
        if (fragmentMap.containsKey(chrName)) {
            return binarySearch(fragmentMap.get(chrName), position);
        } else {
            return -1;
        }


    }

    /**
     * Return fragment that this position lies on.  Fragment 0 means position < sites[0].
     * Fragment 1 means position >= sites[0] and < sites[1].
     *
     * @param sites    The sorted array of fragment sites for the chromosome
     * @param position The position to search for within that array
     * @return The fragment location such that position >= sites[retVal-1] and position <  sites[retVal]
     */
    public static int binarySearch(int[] sites, int position) {
        int lo = 0;
        int hi = sites.length - 1;
        while (lo <= hi) {
            // Base case - found range
            int mid = lo + (hi - lo) / 2;

            if (position > sites[mid]) lo = mid + 1;
            else if (position < sites[mid]) hi = mid - 1;
            else return mid + 1;
        }
        return lo;
    }

}
