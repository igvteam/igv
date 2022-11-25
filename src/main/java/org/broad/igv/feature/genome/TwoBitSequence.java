package org.broad.igv.feature.genome;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.LittleEndianInputStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by jrobinso on 6/13/17.
 */
public class TwoBitSequence implements Sequence {

    static int SIGNATURE_LE = 0x1a412743;
    static int SIGNATURE_BE = 0x4327411a;
    static int HEADER_BLOCK_SIZE = 12500;

    String path;
    private LinkedHashMap<String, Integer> sequenceDataOffsets;

    private HashMap<String, SequenceRecord> sequenceRecordMap;

    public TwoBitSequence(String path) throws IOException {
        this.path = path;
        init();

    }


    /**
     * signature - the number 0x1A412743 in the architecture of the machine that created the file
     * version - zero for now. Readers should abort if they see a version number higher than 0
     * sequenceCount - the number of sequences in the file
     * reserved - always zero for now
     *
     * @throws IOException
     */
    private void init() throws IOException {

        SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        SeekableStream bis = IGVSeekableStreamFactory.getInstance().getBufferedStream(is, 1000);
        LittleEndianInputStream lis = new LittleEndianInputStream(bis);

        int signature = lis.readInt();
        boolean littleEndian = SIGNATURE_LE == signature;

        int version = lis.readInt();   // Should be zero
        int seqCount = lis.readInt();
        int reserved = lis.readInt();    // Should be zero

        sequenceRecordMap = new HashMap<>();
        sequenceDataOffsets = new LinkedHashMap<>();
        for (int i = 0; i < seqCount; i++) {
            byte nameSize = lis.readByte();
            byte[] seqNameBytes = new byte[nameSize];
            lis.readFully(seqNameBytes);
            String seqName = new String(seqNameBytes);
            int offset = lis.readInt();
            sequenceDataOffsets.put(seqName, offset);
        }
    }

    @Override
    public byte[] getSequence(String chr, int start, int end, boolean useCache) {

        try {
            SequenceRecord record = getSequenceRecord(chr);

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
return null;

    }

    private SequenceRecord getSequenceRecord(String chr) throws IOException {

        SequenceRecord record = sequenceRecordMap.get(chr);
        if(record == null) {
            Integer offset = sequenceDataOffsets.get(chr);
            if(offset == null) {
                throw new RuntimeException("Unknown sequence: " + chr);
            }
            SeekableStream is = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
            is.seek(offset);
            SeekableStream bis = IGVSeekableStreamFactory.getInstance().getBufferedStream(is, 1000);
            LittleEndianInputStream lis = new LittleEndianInputStream(bis);
            record = new SequenceRecord(lis, offset);
        }
        return record;

    }

    @Override
    public byte getBase(String chr, int position) {
        return 0;
    }

    @Override
    public List<String> getChromosomeNames() {
        return null;
    }

    @Override
    public int getChromosomeLength(String chrname) {
        return 0;
    }

    @Override
    public boolean isRemote() {
        return FileUtils.isRemote(path);
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
class SequenceRecord {

    final int dnaSize;
    final int[] nBlockStarts;
    final int[] nBlockSizes;
    final int[] maskBlockStarts;
    final int[] maskBlockSizes;

    final int packedDataOffset;

    SequenceRecord(LittleEndianInputStream lis, int offset) throws IOException {

        dnaSize = lis.readInt();

        int nBlockCount = lis.readInt();
        nBlockStarts = new int[nBlockCount];
        for (int i = 0; i < nBlockCount; i++) {
            nBlockStarts[i] = lis.readInt();
        }
        nBlockSizes = new int[nBlockCount];
        for (int i = 0; i < nBlockCount; i++) {
            nBlockSizes[i] = lis.readInt();
        }

        int maskBlockCount = lis.readInt();
        maskBlockStarts = new int[maskBlockCount];
        for (int i = 0; i < maskBlockCount; i++) {
            maskBlockStarts[i] = lis.readInt();
        }
        maskBlockSizes = new int[maskBlockCount];
        for (int i = 0; i < maskBlockCount; i++) {
            maskBlockSizes[i] = lis.readInt();
        }

        packedDataOffset = offset + 12 + nBlockCount * 8 + maskBlockCount * 8;

    }

}
