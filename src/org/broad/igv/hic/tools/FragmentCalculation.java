package org.broad.igv.hic.tools;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.ChromosomeImpl;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Created with IntelliJ IDEA.
 * User: nchernia
 * Date: 8/26/12
 * Time: 7:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class FragmentCalculation {

    Map<Chromosome, int[]> fragmentMap = null;

    public FragmentCalculation(String filename, List<Chromosome> chromosomes) throws IOException{
        InputStream is = null;
        try {
            File file = new File(filename);
            is = new FileInputStream(file);
            Pattern pattern = Pattern.compile("\\s");
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String nextLine;
            fragmentMap = new HashMap<Chromosome, int[]>();

            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = pattern.split(nextLine);
                Chromosome key = null;
                for (Chromosome chromosome: chromosomes) {
                    if (chromosome.getName().equals(tokens[0]))   {
                        key = chromosome;
                        break;
                    }
                }
                if (key == null) {
                    throw new IOException("Can't find chromosome name to match fragment site file.");
                }
                int[] sites = new int[tokens.length - 1];
                for (int i=1; i<tokens.length; i++) {
                    sites[i-1] = Integer.parseInt(tokens[i]);
                }
                fragmentMap.put(key, sites);
            }
        }
        finally {
            is.close();
        }

    }

    public int getNumberFragments(Chromosome chromosome) {
        int[] sites = fragmentMap.get(chromosome);
        if (sites == null) // for "All"
            return 0;
        return sites.length;
    }

    public int getBin(Integer integer, int position) {
        for (Chromosome chromosome : fragmentMap.keySet()) {
            if (chromosome.getIndex() == integer.intValue()) {
                return getBin(chromosome, position);
            }
        }
        return -1;
    }

    public int getBin(Chromosome chromosome, int position) {
        return binarySearch(fragmentMap.get(chromosome), position);
    }

    /**
     * Return fragment that this position lies on.  Fragment 0 means position < sites[0].
     * Fragment 1 means position >= sites[0] and <= sites[1].  (Equal should not happen in practice.)
     * @param sites  The sorted array of fragment sites for the chromosome
     * @param position The position to search for within that array
     * @return    The fragment location such that position > sites[retVal-1] and position < sites[retVal]
     */
    private int binarySearch(int[] sites, int position) {
        int lo = 0;
        int hi = sites.length - 1;
        while (lo <= hi) {
            // Base case - found range
            int mid = lo + (hi - lo) / 2;

            if (position > sites[mid])      lo = mid + 1;
            else if (position < sites[mid]) hi = mid - 1;
            else return mid;
        }
        return lo;
    }

}
