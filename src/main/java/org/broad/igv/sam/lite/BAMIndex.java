package org.broad.igv.sam.lite;

import htsjdk.tribble.util.LittleEndianInputStream;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.util.*;


/**
 * Created by jrobinso on 3/9/17.
 */
public class BAMIndex {


    static int BAI_MAGIC = 21578050;
    static int TABIX_MAGIC = 21578324;
    private static final int SHIFT_AMOUNT = 16;
    private static final int OFFSET_MASK = 0xffff;
    private static final long ADDRESS_MASK = 0xFFFFFFFFFFFFL;


    long firstAlignmentBlock;
    long lastAlignmentBlock;
    Map<Integer, RefIndexPair> refIndexes;
    Map<String, Integer> sequenceIndexMap;   // For tabix use only
    boolean tabix;


    public BAMIndex(Map<Integer, RefIndexPair> refIndexes, long blockMin, long blockMax, Map<String, Integer> sequenceIndexMap) {
        this.firstAlignmentBlock = blockMin;
        this.lastAlignmentBlock = blockMax;
        this.refIndexes = refIndexes;
        this.sequenceIndexMap = sequenceIndexMap;
        this.tabix = false;
    }


    /**
     * Fetch blocks for a particular genomic range.  This method is public so it can be unit-tested.
     *
     * @param refId  the sequence dictionary index of the chromosome
     * @param min  genomic start position
     * @param max  genomic end position
     * @return an array of chunks {minv: {filePointer, offset}, {maxv: {filePointer, offset}}
     */

    public List<Chunk> chunksForRange(int refId, int min, int max) {

        RefIndexPair ba = this.refIndexes.get(refId);
        if (ba == null) {
            return Collections.EMPTY_LIST;
        } else {

            List<Integer> overlappingBins = reg2bins(min, max);        // List of bin #s that overlap min, max
            List<Chunk> chunks = new ArrayList<>();

            // Find chunks in overlapping bins.  Leaf bins (< 4681) are not pruned

            for(Integer bin : overlappingBins) {
                ArrayList<Chunk> binChunks = ba.binIndex.get(bin);
                if (binChunks != null) {
                    chunks.addAll(binChunks);
                }
            }

            // Use the linear index to find minimum file position of chunks that could contain alignments in the region
            int nintv = ba.linearIndex.size();
            VPointer lowest = null;
            int minLin = Math.min(min >> 14, nintv - 1);
            int maxLin = Math.min(max >> 14, nintv - 1);
            for (int i = minLin; i <= maxLin; ++i) {
                VPointer vp = ba.linearIndex.get(i);
                if (vp != null) {
                    // todo -- I think, but am not sure, that the values in the linear index have to be in increasing order.  So the first non-null should be minimum
                    if (lowest == null || vp.isLessThan(lowest)) {
                        lowest = vp;
                    }
                }
            }

            return optimizeChunks(chunks, lowest);
        }

    }

    ;

    List<Chunk> optimizeChunks(List<Chunk> chunks, VPointer lowest) {

        ArrayList mergedChunks = new ArrayList<>();
        Chunk lastChunk = null;

        if (chunks.size() ==  0) return chunks;

        chunks.sort((c0, c1) -> {
            long dif = c0.start.block - c1.start.block;
            if (dif > 0) return 1;
            else if (dif < 0) return -1;
            else {
                return c0.start.offset - c1.start.offset;
            }
        });

        for(Chunk chunk : chunks) {

            if (chunk.end.isGreaterThan(lowest)) {
                if (lastChunk ==  null) {
                    mergedChunks.add(chunk);
                    lastChunk = chunk;
                } else {
                    if ((chunk.start.block - lastChunk.end.block) < 65000) { // Merge chunks that are withing 65k of each other
                        if (chunk.end.isGreaterThan(lastChunk.end)) {
                            lastChunk.end = chunk.end;
                        }
                    } else {
                        mergedChunks.add(chunk);
                        lastChunk = chunk;
                    }
                }
            }
        }

        return mergedChunks;
    }

    /**
     * Calculate the list of bins that overlap with region [beg, end]
     */
    List<Integer> reg2bins(int beg, int end) {

        int k;

        List<Integer> list = new ArrayList<>();

        if (end >= 1 << 29) end = 1 << 29;
        --end;

        list.add(0);
        for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); ++k) list.add(k);
        for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); ++k) list.add(k);
        for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); ++k) list.add(k);
        for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); ++k) list.add(k);
        for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); ++k) list.add(k);
        return list;
    }


    static VPointer readVPointer(LittleEndianInputStream is) throws IOException {

        long vp = is.readLong();
        long block =  (vp >> SHIFT_AMOUNT) & ADDRESS_MASK;
        int offset =  (int) (vp & OFFSET_MASK);

        return new VPointer(block, offset);

    }

    public static class VPointer {
        public long block;
        public int offset;

        public VPointer(long block, int offset) {
            this.block = block;
            this.offset = offset;
        }

        boolean isLessThan(VPointer vp) {
            return this.block < vp.block ||
                    (this.block == vp.block && this.offset < vp.offset);
        }

        boolean isGreaterThan (VPointer vp) {
            return this.block > vp.block ||
                    (this.block == vp.block && this.offset > vp.offset);
        }

    }

    public static class Chunk {

        public VPointer start;
        public VPointer end;

        public Chunk(VPointer start, VPointer end) {
            this.end = end;
            this.start = start;
        }
    }

    public static class RefIndexPair {

        public Map<Integer, ArrayList<Chunk>> binIndex;
        public ArrayList<VPointer> linearIndex;

        public RefIndexPair(Map<Integer, ArrayList<Chunk>> binIndex, ArrayList<VPointer> linearIndex) {
            this.binIndex = binIndex;
            this.linearIndex = linearIndex;
        }
    }


    public static BAMIndex loadIndex(String indexURL,
                              Genome genome) throws IOException {


        boolean tabix = false;

        LittleEndianInputStream parser = new LittleEndianInputStream(new BufferedInputStream(ParsingUtils.openInputStream(indexURL)));

        // NOTE:  if we extend support to tabix must gunzip the file
//                if (tabix) {
//                    var inflate = new Zlib.Gunzip(new Uint8Array(arrayBuffer));
//                    arrayBuffer = inflate.decompress().buffer;
//                }

        long blockMin = Long.MAX_VALUE;
        long blockMax = 0;

        Map<Integer, RefIndexPair> refIndexes = new HashMap<>();
        Map<String, Integer> sequenceIndexMap = null;

        int magic = parser.readInt();

        if (magic == BAI_MAGIC || (tabix && magic == TABIX_MAGIC)) {

            int nref = parser.readInt();

            if (tabix) {
                // Tabix header parameters aren't used, but they must be read to advance the pointer
                int format = parser.readInt();
                int col_seq = parser.readInt();
                int col_beg = parser.readInt();
                int col_end = parser.readInt();
                int meta = parser.readInt();
                int skip = parser.readInt();
                int l_nm = parser.readInt();

                sequenceIndexMap = new HashMap<>();

                for (int i = 0; i < nref; i++) {

                    String seq_name = parser.readString();

                    // Translate to "official" chr name.
                    if (genome != null) {
                        seq_name = genome.getCanonicalChrName(seq_name);
                    }

                    sequenceIndexMap.put(seq_name, i);
                }
            }

            for (int ref = 0; ref < nref; ref++) {

                Map<Integer, ArrayList<Chunk>> binIndex = new HashMap<>();
                ArrayList<VPointer> linearIndex = new ArrayList<>();

                int nbin = parser.readInt();

                for (int b = 0; b < nbin; b++) {

                    int binNumber = parser.readInt();

                    if (binNumber == 37450) {
                        // This is a psuedo bin, not used but we have to consume the bytes
                        int nchnk = parser.readInt(); // # of chunks for this bin
                        readVPointer(parser);   // unmapped beg
                        readVPointer(parser);  // unmapped end
                        long n_maped = parser.readLong();
                        long nUnmapped = parser.readLong();

                    } else {

                        ArrayList<Chunk> chunks = new ArrayList<>();
                        binIndex.put(binNumber, chunks);
                        int nchnk = parser.readInt(); // # of chunks for this bin

                        for (int i = 0; i < nchnk; i++) {
                            VPointer cs = readVPointer(parser);      //chunk_beg
                            VPointer ce = readVPointer(parser);     //chunk_end


                            if (cs.block < blockMin) {
                                blockMin = cs.block;    // Block containing first alignment
                            }
                            if (ce.block > blockMax) {
                                blockMax = ce.block;
                            }
                            chunks.add(new Chunk(cs, ce));
                        }
                    }
                }


                int nintv = parser.readInt();
                for (int i = 0; i < nintv; i++) {
                    VPointer cs = readVPointer(parser);
                    linearIndex.add(cs);   // Might be null
                }

                if (nbin > 0) {
                    refIndexes.put(ref, new RefIndexPair(binIndex, linearIndex));
                }
            }


        } else {
            throw new RuntimeException(indexURL + " is not a " + (tabix ? "tabix" : "bai") + " file");
        }

        return new BAMIndex(refIndexes, blockMin, blockMax, sequenceIndexMap);

    }

}
