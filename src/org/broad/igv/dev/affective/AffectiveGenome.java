package org.broad.igv.dev.affective;

import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.ChromosomeCoordinate;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.FeatureTrack;

import java.text.ParseException;
import java.util.*;

/**
 * @author Jim Robinson
 * @date 1/21/12
 */
public class AffectiveGenome implements Genome {

    TreeMap<Date, Chromosome> chromosomeMap;
    long length;


    private Map<String, Date> chrDateMap = new HashMap<String, Date>();
    private Map<String, Long> cumulativeOffsets = new HashMap();


    public AffectiveGenome() {
        chromosomeMap = new TreeMap<Date, Chromosome>(new Comparator<Date>() {
            public int compare(Date date, Date date1) {
                if (date == null || date1 == null) {
                    return 0;
                }
                return date.compareTo(date1);
            }
        });
        length = 0;

        //addChromosome(new AffectiveChromosome("2011-04-06"));
    }


    /**
     * The "AffectiveGenome" represents the calendar, so any valid date string is a chromosome.  If its not
     * been added yet add it.
     *
     * @param chrName
     * @return
     */
    public Chromosome getChromosome(String chrName) {

        if (!chrDateMap.containsKey(chrName)) {
            // Parse the name to check legality
            try {
                AffectiveChromosome.dateFormat.parse(chrName);
            } catch (ParseException e) {
                // e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                return null;
            }
            createChromosome(chrName);
        }

        Date date = chrDateMap.get(chrName);
        return chromosomeMap.get(date);

    }

    public void createChromosome(String chrName) {

        if(chrDateMap.containsKey(chrName)) {
            return; // Already added
        }

        AffectiveChromosome chromosome = new AffectiveChromosome(chrName);
        chrDateMap.put(chromosome.getName(), chromosome.date);
        chromosomeMap.put(chromosome.date, chromosome);
        length += chromosome.getLength();

        if(!Globals.isHeadless()) {
            // Force update of command bar

        }
    }

    public String getId() {
        return "affective";
    }

    public String getHomeChromosome() {
        if (getLongChromosomeNames().size() == 1) {
            return getLongChromosomeNames().get(0);
        } else {
            return Globals.CHR_ALL;
        }
    }


    public Collection<Chromosome> getChromosomes() {
        return chromosomeMap.values();
    }

    public List<String> getAllChromosomeNames() {
        ArrayList<String> names = new ArrayList<String>(chromosomeMap.size());
        for (Chromosome chromosome : chromosomeMap.values()) {
            names.add(chromosome.getName());
        }
        return names;
    }

    public String getChromosomeAlias(String str) {
        return str;
    }

    @Override
    public long getTotalLength() {
        return length;
    }

    public long getNominalLength() {
        return length;
    }

    public long getCumulativeOffset(String chr) {

        Long cumOffset = cumulativeOffsets.get(chr);
        if (cumOffset == null) {
            long offset = 0;
            for (Chromosome chromosome : chromosomeMap.values()) {
                if (chr.equals(chromosome.getName())) {
                    break;
                }
                offset += chromosome.getLength();
            }
            cumOffset = new Long(offset);
            cumulativeOffsets.put(chr, cumOffset);
        }
        return cumOffset.longValue();
    }

    /**
     * Covert the chromosome coordinate in BP to genome coordinates in KBP
     *
     * @param chr
     * @param locationBP
     * @return
     */
    public int getGenomeCoordinate(String chr, int locationBP) {
        return (int) ((getCumulativeOffset(chr) + locationBP) / 1000);
    }

    /**
     * Convert the genome coordinates in KBP to a chromosome coordinate
     */
    public ChromosomeCoordinate getChromosomeCoordinate(int genomeKBP) {

        long cumOffset = 0;
        for (Chromosome chromosome : chromosomeMap.values()) {
            int chrLen = chromosome.getLength();
            if ((cumOffset + chrLen) / 1000 > genomeKBP) {
                int bp = (int) (genomeKBP * 1000 - cumOffset);
                return new ChromosomeCoordinate(chromosome.getName(), bp);
            }
            cumOffset += chrLen;
        }

        String c = chromosomeMap.lastEntry().getValue().getName();
        int bp = (int) (genomeKBP - cumOffset) * 1000;
        return new ChromosomeCoordinate(c, bp);
    }

    public String getNextChrName(String chr) {
        List<String> chrList = getAllChromosomeNames();
        for (int i = 0; i < chrList.size() - 1; i++) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i + 1);
            }
        }
        return null;
    }

    public String getPrevChrName(String chr) {
        List<String> chrList = getAllChromosomeNames();
        for (int i = chrList.size() - 1; i > 0; i--) {
            if (chrList.get(i).equals(chr)) {
                return chrList.get(i - 1);
            }
        }
        return null;
    }

    /**
     * Return the nucleotide sequence on the + strand for the genomic interval.  Not relevant for
     * Affective genome.
     */
    public byte[] getSequence(String chr, int start, int end) {
        return null;
    }

    public String getDisplayName() {
        return "Affective";
    }

    public byte getReference(String chr, int pos) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public FeatureTrack getGeneTrack() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public List<String> getLongChromosomeNames() {
        return getAllChromosomeNames();
    }

    public void addChrAliases(Map<String, String> aliases) {
        //To change body of implemented methods use File | Settings | File Templates.
    }
    public void loadUserDefinedAliases() {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
