package org.broad.igv.feature.tribble.reader;

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.LittleEndianInputStream;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class BamIndex {

    // Represents a BAM index.
    // Code is based heavily on bam.js, part of the Dalliance Genome Explorer,  (c) Thomas Down 2006-2001.

    int BAI_MAGIC = 21578050;
    int TABIX_MAGIC = 21578324;

    int MB = 1000000;
    String[] seqNames;
    Map<String, Integer> sequenceIndexMap;
    private boolean tabix;
    private long firstBlockPosition;
    private long lastBlockPosition;

    private Map<Integer, SeqIndex> indeces;


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
        Map<String, Integer> sequenceIndexMap = new HashMap<>();


        if (!(magic == BAI_MAGIC || magic == TABIX_MAGIC)) {
            throw new RuntimeException("Unrecognized magic value: " + magic);
        }

        long blockMin = Long.MAX_VALUE;
        long blockMax = 0l;
        String[] seqNames;

        this.tabix = magic == TABIX_MAGIC;

        int nref = lis.readInt();
        seqNames = new String[nref];

        if (tabix) {
            // Tabix header parameters aren't used, but they must be read to advance the pointer
            int format = lis.readInt();
            int col_seq = lis.readInt();
            int col_beg = lis.readInt();
            int col_end = lis.readInt();
            int meta = lis.readInt();
            int skip = lis.readInt();
            int l_nm = lis.readInt();

            for (int i = 0; i < nref; i++) {

                String seq_name = lis.readString();
                // Translate to "official" chr name.
                if (genome != null) {
                    seq_name = genome.getCanonicalChrName(seq_name);
                }
                sequenceIndexMap.put(seq_name, i);
                seqNames[i] = seq_name;
            }
        }

        this.firstBlockPosition = blockMin;
        this.lastBlockPosition = blockMax;
        this.seqNames = seqNames;
        this.sequenceIndexMap = sequenceIndexMap;
        this.indeces = new HashMap<>();

        // Loop through sequences
        for (int ref = 0; ref < nref; ref++) {

            int nbin = lis.readInt();
            Map<Integer, List<Chunk>> binIndex = new HashMap<>();
            List<VPointer> linearIndex = new ArrayList<>();

            for (int b = 0; b < nbin; b++) {
                int binNumber = lis.readInt();
                if (binNumber == 37450) {
                    int nchnk = lis.readInt();
                    // This is a psuedo bin, not used but we have to consume the bytes
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

            int nintv = lis.readInt();
            for (int i = 0; i < nintv; i++) {
                VPointer cs = new VPointer(lis.readLong());
                linearIndex.add(cs);   // Might be null
            }

            if (nbin > 0) {
                this.indeces.put(ref, new SeqIndex(binIndex, linearIndex));
            }
        }


        this.firstBlockPosition = blockMin;
        this.lastBlockPosition = blockMax;
        this.sequenceIndexMap = sequenceIndexMap;
        this.tabix = tabix;

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
     *
     * @return List of int arrays, each defining a bin range
     */
    List<int[]> reg2bins(int beg, int end) {
        List<int[]> list = new ArrayList<>();
        if (end >= 1 << 29) end = 1 << 29;
        --end;
        list.add(new int[]{0, 0});
        list.add(new int[]{1 + (beg >> 26), 1 + (end >> 26)});
        list.add(new int[]{9 + (beg >> 23), 9 + (end >> 23)});
        list.add(new int[]{73 + (beg >> 20), 73 + (end >> 20)});
        list.add(new int[]{585 + (beg >> 17), 585 + (end >> 17)});
        list.add(new int[]{4681 + (beg >> 14), 4681 + (end >> 14)});
        return list;
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
        List<VPointer> linearIndex;

        public SeqIndex(Map<Integer, List<Chunk>> binIndex, List<VPointer> linearIndex) {
            this.binIndex = binIndex;
            this.linearIndex = linearIndex;
        }
    }

}

