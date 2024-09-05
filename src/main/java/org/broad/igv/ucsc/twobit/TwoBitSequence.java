package org.broad.igv.ucsc.twobit;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.feature.genome.SequenceNotFoundException;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ucsc.BPIndex;
import org.broad.igv.ucsc.BPTree;

import java.io.IOException;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;


/**
 * Reader for UCSC ".2bit" sequence files. Reference: https://genome.ucsc.edu/FAQ/FAQformat.html#format7
 * Note: Some portions of this code were adapated from the GMOD two-bit.js project, @Copyright (c) 2017 Robert Buels
 * https://github.com/GMOD/twobit-js/blob/master/src/twoBitFile.ts*
 */


public class TwoBitSequence implements Sequence {

    private static Logger log = LogManager.getLogger(TwoBitSequence.class);

    // the number 0x1A412743 in the architecture of the machine that created the file
    static int SIGNATURE = 0x1a412743;
    String path;
    private HashMap<String, SequenceRecord> sequenceRecordCache;
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;  // Until proven otherwise
    private int seqCount;
    BPIndex index;


    public TwoBitSequence(String path) throws IOException {
        init(path);
        index = new TwoBitIndex(path, this.byteOrder, this.seqCount);
    }

    public TwoBitSequence(String path, String indexPath) throws IOException {
        init(path);
        index = BPTree.loadBPTree(indexPath, 0);
    }

    UnsignedByteBufferImpl loadBinaryBuffer(long start, int size) throws IOException {
        return UnsignedByteBufferImpl.loadBinaryBuffer(path, byteOrder, start, size);
    }

    private void init(String path) throws IOException {

        this.path = path;
        this.sequenceRecordCache = new HashMap<>();

        long filePosition = 0;
        UnsignedByteBuffer buffer = UnsignedByteBufferImpl.loadBinaryBuffer(path, byteOrder, filePosition, 64);

        int signature = buffer.getInt();
        if (SIGNATURE != signature) {
            this.byteOrder = ByteOrder.BIG_ENDIAN;
            buffer.position(0);
            signature = buffer.getInt();
            if (SIGNATURE != signature) {
                throw new RuntimeException("Unexpected magic number");
            }
        }

        final int version = buffer.getInt();   // Should be zero
        this.seqCount = buffer.getInt();
        final int reserved = buffer.getInt();    // Should be zero

    }

    @Override
    public byte[] getSequence(String chr, int start, int end) {
        return readSequence(chr, start, end);
    }

    @Override
    public byte getBase(String chr, int position) {
        throw new RuntimeException("getBase is not implementd for TwoBitSequence");
    }

    @Override
    public List<String> getChromosomeNames() {
        return null;
    }

    @Override
    public int getChromosomeLength(String seq) {
        try {
            SequenceRecord sequenceRecord = getSequenceRecord(seq);
            return sequenceRecord.getDnaSize();
        } catch (SequenceNotFoundException e) {
            return -1;
        } catch (IOException e) {
            log.error("Error reading sequence " + seq, e);
            return -1;
        }
    }

    @Override
    public List<Chromosome> getChromosomes() {
        return null;
    }

    @Override
    public boolean hasChromosomes() {
        return false;
    }

    /**
     * Read the sequence requested.  Returns null if seqName is unknown
     *
     * @param seqName
     * @param regionStart
     * @param regionEnd
     * @return
     */
    public byte[] readSequence(String seqName, int regionStart, int regionEnd) {

        try {

            SequenceRecord record = getSequenceRecord(seqName);
            if (record == null) {
                return null;
            }

            if (regionStart < 0) {
                throw new RuntimeException("regionStart cannot be less than 0");
            }

            Queue<Block> nBlocks = _getOverlappingBlocks(regionStart, regionEnd, record.nBlocks);
            Queue<Block> maskBlocks = _getOverlappingBlocks(regionStart, regionEnd, record.maskBlocks);

            int baseBytesOffset = regionStart / 4;     // "int" division will automatically floor
            long start = record.packedPos + baseBytesOffset;
            int size = regionEnd / 4 - baseBytesOffset + 1;

            UnsignedByteBufferImpl buffer = loadBinaryBuffer(start, size);
            byte[] baseBytes = buffer.array();

            //new byte[size];
            //buffer.get(baseBytes)
            //this.is.seek(start);
            // this.is.readFully(baseBytes);

            byte[] sequenceBases = new byte[regionEnd - regionStart];
            for (int genomicPosition = regionStart; genomicPosition < regionEnd; genomicPosition++) {

                // function checks if  we are currently masked
                while (maskBlocks.size() > 0 && maskBlocks.peek().end <= genomicPosition) {
                    maskBlocks.remove();
                }
                Block mBlock = maskBlocks.peek();
                boolean baseIsMaked = mBlock != null && mBlock.start <= genomicPosition && mBlock.end > genomicPosition;

                // process the N block if we have one.  Masked "N" ("n")  is not supported
                Block firstBlock = nBlocks.peek();
                if (firstBlock != null && genomicPosition >= firstBlock.start && genomicPosition < firstBlock.end) {
                    Block currentNBlock = nBlocks.remove();
                    while (genomicPosition < currentNBlock.end && genomicPosition < regionEnd) {
                        sequenceBases[genomicPosition - regionStart] = 'N';
                        genomicPosition++;
                    }
                    genomicPosition--;
                } else {
                    int bytePosition = (genomicPosition / 4) - baseBytesOffset;
                    int subPosition = genomicPosition % 4;
                    int s = Byte.toUnsignedInt(baseBytes[bytePosition]);
                    int idx = genomicPosition - regionStart;
                    sequenceBases[idx] = baseIsMaked ? maskedByteTo4Bases[s][subPosition] : byteTo4Bases[s][subPosition];
                }
            }
            return sequenceBases;

        } catch (IOException e) {
            throw new RuntimeException(e);
        }


    }

    public SequenceRecord getSequenceRecord(String seqName) throws IOException {

        SequenceRecord record = this.sequenceRecordCache.get(seqName);

        if (record == null) {
            long[] offset_length = this.index.search(seqName);
            if (offset_length == null) {
                throw new SequenceNotFoundException("Unknown sequence: " + seqName);
            }
            long offset = offset_length[0];

            // Read size of dna data & # of "N" blocks
            int size = 8;
            UnsignedByteBuffer buffer = loadBinaryBuffer(offset, size);
            int dnaSize = buffer.getInt();
            int nBlockCount = buffer.getInt();
            offset += size;

            // Read "N" blocks and # of mask blocks
            size = nBlockCount * (4 + 4) + 4;
            buffer = loadBinaryBuffer(offset, size);

            int[] nBlockStarts = new int[nBlockCount];
            for (int i = 0; i < nBlockCount; i++) {
                nBlockStarts[i] = buffer.getInt();
            }
            int[] nBlockSizes = new int[nBlockCount];
            for (int i = 0; i < nBlockCount; i++) {
                nBlockSizes[i] = buffer.getInt();
            }
            int maskBlockCount = buffer.getInt();
            offset += size;

            size = maskBlockCount * (4 + 4) + 4;
            buffer = loadBinaryBuffer(offset, size);

            int[] maskBlockStarts = new int[maskBlockCount];
            for (int i = 0; i < maskBlockCount; i++) {
                maskBlockStarts[i] = buffer.getInt();
            }
            int[] maskBlockSizes = new int[maskBlockCount];
            for (int i = 0; i < maskBlockCount; i++) {
                maskBlockSizes[i] = buffer.getInt();
            }

            // Transform "N" and "mask" block data into something more useful
            //Transform "N" and "mask" block data into something more useful
            Block[] nBlocks = new Block[nBlockCount];
            for (int i = 0; i < nBlockCount; i++) {
                nBlocks[i] = new Block(nBlockStarts[i], nBlockSizes[i]);
            }
            Block[] maskBlocks = new Block[maskBlockCount];
            for (int i = 0; i < maskBlockCount; i++) {
                maskBlocks[i] = new Block(maskBlockStarts[i], maskBlockSizes[i]);
            }


            int reserved = buffer.getInt();
            if (reserved != 0) {
                throw new RuntimeException("Bad 2-bit file");
            }

            long packedPos = offset + size;
            record = new SequenceRecord(dnaSize, nBlocks, maskBlocks, packedPos);

            sequenceRecordCache.put(seqName, record);
        }
        return record;

    }

    /**
     * Return blocks overlapping the genome region [start, end]
     * <p>
     * TODO -- optimize this, currently it uses linear search
     * * *
     *
     * @param start
     * @param end
     * @param blocks
     * @returns {*[]}
     * @private
     */
    static Queue<Block> _getOverlappingBlocks(long start, long end, Block[] blocks) {

        Queue<Block> overlappingBlocks = new LinkedList<>();
        for (Block block : blocks) {
            if (block.start > end) {
                break;
            } else if (block.end < start) {
                continue;
            } else {
                overlappingBlocks.add(block);
            }
        }
        return overlappingBlocks;
    }

    static byte[] twoBit = {'T', 'C', 'A', 'G'};
    static byte[][] byteTo4Bases = new byte[256][4];
    static byte[][] maskedByteTo4Bases = new byte[256][4];

    static {
        for (int i = 0; i < 256; i++) {
            byteTo4Bases[i] =
                    new byte[]{twoBit[(i >> 6) & 3],
                            twoBit[(i >> 4) & 3],
                            twoBit[(i >> 2) & 3],
                            twoBit[i & 3]};
        }
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 4; j++) {
                maskedByteTo4Bases[i][j] = (byte) Character.toLowerCase(byteTo4Bases[i][j]);
            }

        }

    }

}

/*
dnaSize - number of bases of DNA in the sequence
nBlockCount - the number of blocks of Ns in the file (representing unknown sequence)
nBlockStarts - an array of length nBlockCount of 32 bit integers indicating the (0-based) starting position of a block of Ns
nBlockSizes - an array of length nBlockCount of 32 bit integers indicating the length of a block of Ns
maskBlockCount - the number of masked (lower-case) blocks
maskBlockStarts - an array of length maskBlockCount of 32 bit integers indicating the (0-based) starting position of a masked block
maskBlockSizes - an array of length maskBlockCount of 32 bit integers indicating the length of a masked block
reserved - always zero for now

 */

