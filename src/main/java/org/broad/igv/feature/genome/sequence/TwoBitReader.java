package org.broad.igv.feature.genome.sequence;

/**
 * Reader for UCSC ".2bit" sequence files. Reference: https://genome.ucsc.edu/FAQ/FAQformat.html#format7
 * Note: Some portions of this code were adapated from the GMOD two-bit.js project, @Copyright (c) 2017 Robert Buels
 * https://github.com/GMOD/twobit-js/blob/master/src/twoBitFile.ts*
 */


import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * Created by jrobinso on 6/13/17.
 */
public class TwoBitReader {

    // the number 0x1A412743 in the architecture of the machine that created the file
    static int SIGNATURE = 0x1a412743;


    private LinkedHashMap<String, Integer> sequenceDataOffsets;
    private HashMap<String, SequenceRecord> sequenceRecordMap;
    SeekableStream is;
    ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;  // Until proven otherwise

    public TwoBitReader(String path) throws IOException {
        this.is = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        init();
    }

    ByteBuffer loadBinaryBuffer(long start, int size) throws IOException {
        ByteBuffer bb = ByteBuffer.allocate(size);
        bb.order(this.byteOrder);
        byte[] bytes = bb.array();
        this.is.seek(start);
        this.is.readFully(bytes);
        return bb;
    }

    private void init() throws IOException {
        readIndex();
    }


    /**
     * signature - the number 0x1A412743 in the architecture of the machine that created the file
     * version - zero for now. Readers should abort if they see a version number higher than 0
     * sequenceCount - the number of sequences in the file
     * reserved - always zero for now
     *
     * @throws IOException
     */
    private void readIndex() throws IOException {

        long filePosition = 0;
        ByteBuffer buffer = loadBinaryBuffer(filePosition, 64);

        int signature = buffer.getInt();
        if (SIGNATURE != signature) {
            this.byteOrder = ByteOrder.BIG_ENDIAN;
            buffer = ByteBuffer.allocate(64);
            buffer = loadBinaryBuffer(0, 64);
            signature = buffer.getInt();
            if (SIGNATURE != signature) {
                throw new RuntimeException("Unexpected magic number");
            }
        }

        final int version = buffer.getInt();   // Should be zero
        final int seqCount = buffer.getInt();
        final int reserved = buffer.getInt();    // Should be zero

        // Loop through sequences loading name and file offset.  We don't know the precise size in bytes in advance
        // so we need to check for bytes available and reload as needed.
        final int estNameLength = 20;
        sequenceRecordMap = new HashMap<>();
        sequenceDataOffsets = new LinkedHashMap<>();
        for (int i = 0; i < seqCount; i++) {

            if (buffer.remaining() < 1) {
                filePosition += buffer.position();
                final int estSize = (seqCount - i) * estNameLength + 100;
                buffer = loadBinaryBuffer(filePosition, estSize);
            }

            final byte nameSize = buffer.get();

            if (buffer.remaining() < nameSize * 5) {
                filePosition += buffer.position();
                final int estSize = (seqCount - i) * estNameLength + 100;
                buffer = loadBinaryBuffer(filePosition, estSize);
            }

            byte[] seqNameBytes = new byte[nameSize];
            buffer.get(seqNameBytes);
            String seqName = new String(seqNameBytes);

            int offset = buffer.getInt();
            sequenceDataOffsets.put(seqName, offset);
        }
    }

    public List<String> getSequenceNames() {
        return new ArrayList<>(sequenceDataOffsets.keySet());
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
            if (sequenceDataOffsets == null) {
                readIndex();
            }

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


            ByteBuffer buffer = loadBinaryBuffer(start, size);
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

    SequenceRecord getSequenceRecord(String seqName) throws IOException {

        SequenceRecord record = sequenceRecordMap.get(seqName);

        if (record == null) {
            Integer offset = sequenceDataOffsets.get(seqName);
            if (offset == null) {
                throw new RuntimeException("Unknown sequence: " + seqName);
            }

            // Read size of dna data & # of "N" blocks
            int size = 8;
            ByteBuffer buffer = loadBinaryBuffer(offset, size);
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

            sequenceRecordMap.put(seqName, record);
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

    static class SequenceRecord {

        final int dnaSize;
        Block[] nBlocks;
        Block[] maskBlocks;

        final long packedPos;

        SequenceRecord(int dnaSize, Block[] nBlocks, Block[] maskBlocks, long packedPos) throws IOException {
            this.dnaSize = dnaSize;
            this.nBlocks = nBlocks;
            this.maskBlocks = maskBlocks;
            this.packedPos = packedPos;

        }
    }

    static class Block {
        final long start;
        final int size;
        final long end;

        Block(long start, int size) {
            this.start = start;
            this.size = size;
            this.end = start + size;
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

