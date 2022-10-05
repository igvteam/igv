package org.broad.igv.feature.tribble.reader;

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.LittleEndianInputStream;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.*;

public class CSIIndex {

    // Represents a BAM index.
    // Code is based heavily on bam.js, part of the Dalliance Genome Explorer,  (c) Thomas Down 2006-2001.

    static final int CSI1_MAGIC = 21582659; // CSI\1
    static final int CSI2_MAGIC = 38359875; // CSI\2

    int MB = 1000000;
    String[] seqNames;
    Map<String, Integer> sequenceIndexMap;
    private boolean tabix;
    private long firstBlockPosition;
    private long lastBlockPosition;

    private Map<Integer, SeqIndex> indeces;
    private int minShift;
    private int depth;
    private int maxBinNumber;


    public void parse(byte[] arrayBuffer, Genome genome) throws IOException {

        // TODO -- handle big endian?
        LittleEndianInputStream lis;

        // Some indexs are bgzipped, specifically tabix, and csi.  Bam (bai) are not.  Tribble is usually not.
        // Check first 2 bytes of file for gzip magic number, and inflate if neccessary
        if (arrayBuffer[0] == 31 && arrayBuffer[1] == -117) {    // gzipped
            lis = new LittleEndianInputStream(new BlockCompressedInputStream(new ByteArrayInputStream(arrayBuffer)));
        } else {
            lis = new LittleEndianInputStream(new ByteArrayInputStream(arrayBuffer));
        }

        int magic = lis.readInt();

        if (magic != CSI1_MAGIC) {
            if (magic == CSI2_MAGIC) {
                throw new RuntimeException("CSI version 2 is not supported.  Please enter an issue at https://github.com/igvteam/igv.js");
            } else {
                throw new RuntimeException("Not a CSI index");
            }
        }

        //const indices = []
        long blockMin = Long.MAX_VALUE;
        long blockMax = 0l;

        Map<String, Integer> sequenceIndexMap = new HashMap<>();

        this.minShift = lis.readInt();
        this.depth = lis.readInt();
        int lAux = lis.readInt();

        this.maxBinNumber = (((1 << (this.depth + 1) * 3) - 1) / 7) + 1;
        List<String> seqNames = new ArrayList<>();
        if (lAux >= 28) {
            // Tabix header parameters aren't used, but they must be read to advance the pointer
            int format = lis.readInt();
            int col_seq = lis.readInt();
            int col_beg = lis.readInt();
            int col_end = lis.readInt();
            int meta = lis.readInt();
            int skip = lis.readInt();
            int l_nm = lis.readInt();

            int i = 0;
            int nameCharCount = 0;
            while (nameCharCount < l_nm) {
                String seq_name = lis.readString();
                nameCharCount += seq_name.length() + 1;
                // Translate to "official" chr name.
                if (genome != null) {
                    seq_name = genome.getCanonicalChrName(seq_name);
                }
                sequenceIndexMap.put(seq_name, i);
                seqNames.add(seq_name);
                i++;
            }
        }

        this.indeces = new HashMap<>();

        // Loop through sequence
        int nref = lis.readInt();
        for (int ref = 0; ref < nref; ref++) {

            int nbin = lis.readInt();
            Map<Integer, List<Chunk>> binIndex = new HashMap<>();
            Map<Integer, VPointer> linearOffset = new HashMap<>();

            for (int b = 0; b < nbin; b++) {

                int binNumber = lis.readInt();
                linearOffset.put(binNumber, new VPointer(lis.readLong()));

                if (binNumber > this.maxBinNumber) {
                    // This is a psuedo bin, not used but we have to consume the bytes
                    int nchnk = lis.readInt();
                    VPointer cs = new VPointer(lis.readLong());   // unmapped beg
                    VPointer ce = new VPointer(lis.readLong());   // unmapped end
                    long n_maped = lis.readLong();
                    long nUnmapped = lis.readLong();
                } else {

                    binIndex.put(binNumber, new ArrayList<>());
                    int nchnk = lis.readInt(); // # of chunks for this bin

                    for (int i = 0; i < nchnk; i++) {
                        VPointer cs = new VPointer(lis.readLong());    //chunk_beg
                        VPointer ce = new VPointer(lis.readLong());    //chunk_end

                        if (cs.block < blockMin) {
                            blockMin = cs.block;    // Block containing first alignment
                        }
                        if (ce.block > blockMax) {
                            blockMax = ce.block;
                        }
                        binIndex.get(binNumber).add(new Chunk(cs, ce));

                    }
                }
            }

            List<Integer> binNumberList = new ArrayList<>(binIndex.keySet());
            Collections.sort(binNumberList);

            if (nbin > 0) {
                this.indeces.put(ref, new SeqIndex(binIndex, linearOffset));
            }
        }


        this.firstBlockPosition = blockMin;
        this.lastBlockPosition = blockMax;
        this.seqNames = seqNames.toArray(new String[seqNames.size()]);
        this.sequenceIndexMap = sequenceIndexMap;

    }


    public Set<String> sequenceNames() {
        return this.sequenceIndexMap.keySet();
    }

    /**
     * Fetch blocks for a particular genomic range.  This method is public so it can be unit-tested.
     *
     * @param refId the sequence dictionary index of the chromosome
     * @param min   genomic start position
     * @param max   genomic end position
     */
    List<Chunk> chunksForRange(int refId, int min, int max) {

        SeqIndex ba = this.indeces.get(refId);

        if (ba == null) {
            return Collections.EMPTY_LIST;
        } else {
            List<int[]> overlappingBins = reg2bins(min, max);        // List of bin ranges that overlap min, max
            if (overlappingBins.size() == 0) return Collections.EMPTY_LIST;

            List<Chunk> chunks = new ArrayList<>();

            // Find chunks in overlapping bins.  Leaf bins (< 4681) are not pruned
            for (int[] binRange : overlappingBins) {
                for (int bin = binRange[0]; bin <= binRange[1]; bin++) {
                    if (ba.binIndex.get(bin) != null) {
                        List<Chunk> binChunks = ba.binIndex.get(bin);
                        chunks.addAll(binChunks);
                    }
                }
            }


            // Use the linear index to find minimum file position of chunks that could contain alignments in the region

            VPointer lowest = ba.linearOffsets.get(overlappingBins.get(0));

            return optimizeChunks(chunks, lowest);
        }
    }


    List<Chunk> optimizeChunks(List<Chunk> chunks, VPointer lowest) {

        List<Chunk> mergedChunks = new ArrayList<>();
        Chunk lastChunk = null;

        if (chunks.isEmpty()) return chunks;

        chunks.sort((c0, c1) -> {
            long dif = c0.left.block - c1.left.block;
            if (dif != 0) {
                return (int) dif;
            } else {
                return c0.left.offset - c1.left.offset;
            }
        });

        for (Chunk chunk : chunks) {

            if (lowest == null || chunk.right.isGreaterThan(lowest)) {
                if (lastChunk == null) {
                    mergedChunks.add(chunk);
                    lastChunk = chunk;
                } else {
                    if (canMerge(lastChunk, chunk)) {
                        if (chunk.right.isGreaterThan(lastChunk.right)) {
                            lastChunk.right = chunk.right;
                        }
                    } else {
                        mergedChunks.add(chunk);
                        lastChunk = chunk;
                    }
                }
            } else {
                //console.log(`skipping chunk ${chunk.minv.block} - ${chunk.maxv.block}`)
            }
        }

        return mergedChunks;
    }

    //
//    /**
//     * Merge 2 blocks if the gap between them is < 1kb and the total resulting size < 100mb
//     *
//     * @param chunk1
//     * @param chunk2
//     * @returns {boolean|boolean}
//     */
    boolean canMerge(Chunk chunk1, Chunk chunk2) {
        long gap = chunk2.left.block - chunk1.right.block;
        long total = chunk2.right.block - chunk1.left.block;
        return gap < 30000 && total < 10 * MB;
    }


    /**
     * Calculate the list of bins that overlap with region [beg, end]
     * reg2bins implementation adapted from GMOD/tabix-js  https://github.com/GMOD/tabix-js/blob/master/src/csi.ts
     *
     * @return List of int arrays, each defining a bin range
     */
    List<int[]> reg2bins(int beg, int end) {

        beg -= 1; // < convert to 1-based closed
        if (beg < 1) beg = 1;
        if (end > Math.pow(2, 34)) end = (int) Math.pow(2, 34) ;// 17 GiB ought to be enough for anybody
        end -= 1;
        int l = 0;
        int t = 0;
        int s = this.minShift + this.depth * 3;


        List<int[]> bins = new ArrayList<>();
        for (; l <= this.depth; s -= 3, t += (1 << l * 3), l += 1) {
            int b = t + (beg >> s);
            int e = t + (end >> s);
            if (e - b + bins.size() > this.maxBinNumber)
                throw new RuntimeException("" + (beg - end) + "is too large for current binning scheme (shift " +
                        this.minShift + ", depth " + this.depth + "), try a smaller query or a coarser index binning scheme");
            //for (let i = b; i <= e; i += 1) bins.push(i)
            bins.add(new int[]{b, e});
        }

        return bins;
    }

    public int getDepth() {
        return this.depth;
    }

    public int getMinShift() {
        return this.minShift;
    }

    public int getMaxBinNumber() {
        return this.maxBinNumber;
    }


    static class Block {
        Chunk minv;
        Chunk maxv;
        int bin;

        public Block(Chunk minv, Chunk maxv, int bin) {
            this.minv = minv;
            this.maxv = maxv;
            this.bin = bin;
        }
    }

    static class Chunk {

        VPointer left;
        VPointer right;

        public Chunk(VPointer left, VPointer right) {
            this.left = left;
            this.right = right;
        }
    }

    static class SeqIndex {
        Map<Integer, List<Chunk>> binIndex;
        Map<Integer, VPointer> linearOffsets;

        public SeqIndex(Map<Integer, List<Chunk>> binIndex, Map<Integer, VPointer> linearOffsets) {
            this.binIndex = binIndex;
            this.linearOffsets = linearOffsets;
        }
    }

}

