package org.broad.igv.hic.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 9/6/11
 */
public class Chromosome {
    
    private int index;
    private String name;
    private int size;

    
    // Map of genome -> chromosomes
    public static Map<String, List<Chromosome>> chromosomes = new HashMap<String, List<Chromosome>>();

    public Chromosome(int index, String name, int size) {
        this.index = index;
        this.name = name;
        this.size = size;
    }

    public int getIndex() {
        return index;
    }

    public String getName() {
        return name;
    }

    public int getSize() {
        return size;
    }

    public String toString() {
        return name;
    }
    
    
    
    static {
        List<Chromosome> hg18 = new ArrayList<Chromosome>();
        hg18.add(new Chromosome(1, "1", 247249719));
        hg18.add(new Chromosome(2, "2", 242951149));
        hg18.add(new Chromosome(3, "3", 199501827));
        hg18.add(new Chromosome(4, "4", 191273063));
        hg18.add(new Chromosome(5, "5", 180857866));
        hg18.add(new Chromosome(6, "6", 170899992));
        hg18.add(new Chromosome(7, "7", 158821424));
        hg18.add(new Chromosome(8, "8", 146274826));
        hg18.add(new Chromosome(9, "9", 140273252));
        hg18.add(new Chromosome(10, "10", 135374737));
        hg18.add(new Chromosome(11, "11", 134452384));
        hg18.add(new Chromosome(12, "12", 132349534));
        hg18.add(new Chromosome(13, "13", 114142980));
        hg18.add(new Chromosome(14, "14", 106368585));
        hg18.add(new Chromosome(15, "15", 100338915));
        hg18.add(new Chromosome(16, "16", 88827254));
        hg18.add(new Chromosome(17, "17", 78774742));
        hg18.add(new Chromosome(18, "18", 76117153));
        hg18.add(new Chromosome(19, "19", 63811651));
        hg18.add(new Chromosome(20, "20", 62435964));
        hg18.add(new Chromosome(21, "21", 46944323));
        hg18.add(new Chromosome(22, "22", 49691432));
        hg18.add(new Chromosome(23, "23", 154913754));
        hg18.add(new Chromosome(24, "24", 57772954));
        hg18.add(new Chromosome(0, "M", 16571));
        chromosomes.put("hg18", hg18);

        List<Chromosome> dmel = new ArrayList<Chromosome>();
        dmel.add(new Chromosome(1, "2L",  23011544));
        dmel.add(new Chromosome(2, "2R", 21146708));  // 2R
        dmel.add(new Chromosome(3, "3L", 24543557)); // 3L
        dmel.add(new Chromosome(4, "3R", 27905053)); //3R
        dmel.add(new Chromosome(5, "4", 1351857)); // 4
        dmel.add(new Chromosome(6, "X", 22422827)); // X
        dmel.add(new Chromosome(7, "U", 10049037)); // U
        chromosomes.put("dmel", dmel);
    }
}
