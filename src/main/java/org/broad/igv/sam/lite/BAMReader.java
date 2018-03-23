package org.broad.igv.sam.lite;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.util.LittleEndianInputStream;
import org.apache.log4j.Logger;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.sam.*;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;
import java.util.zip.DataFormatException;

/**
 * Created by jrobinso on 3/9/17.
 */
public class BAMReader implements AlignmentReader<Alignment> {

    private static Logger log = Logger.getLogger(BAMReader.class);
    private final String path;
    private final String indexPath;
    int BAM_MAGIC = 21840194;
    int BAI_MAGIC = 21578050;
    char[] SECRET_DECODER = {'=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'};
    char[] CIGAR_DECODER = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'};

    int PAIRED_FLAG = 0x1;
    int READ_STRAND_FLAG = 0x10;
    int MATE_STRAND_FLAG = 0x20;

    int MATE_IS_MAPPED_FLAG = 0x8;
    int FIRST_OF_PAIR_FLAG = 0x40;
    int SECOND_OF_PAIR_FLAG = 0x80;
    int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
    int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    int DUPLICATE_READ_FLAG = 0x400;
    int SUPPLEMENTARY_FLAG = 0x800;

    int MAX_GZIP_BLOCK_SIZE = 65536;   //  APPARENTLY.  Where is this documented???
    int DEFAULT_SAMPLING_WINDOW_SIZE = 100;
    int DEFAULT_SAMPLING_DEPTH = 50;
    int MAXIMUM_SAMPLING_DEPTH = 2500;

    BAMIndex bamIndex = null;

    Genome genome;
    Map<String, Integer> chrToIndex;
    private String[] indexToChr;


    public BAMReader(String path) {

        this.path = path;
        this.indexPath = path + ".bai";

//        this.config = config;
//
//        this.filter = config.filter || new igv.BamFilter();
//
//        this.bamPath = config.url;
//        // Todo - deal with Picard convention.  WHY DOES THERE HAVE TO BE 2?
//        this.baiPath = config.indexURL || this.bamPath + ".bai"; // If there is an indexURL provided, use it!
//        this.headPath = config.headURL || this.bamPath;
//
//
//        this.samplingWindowSize = config.samplingWindowSize == = undefined ? DEFAULT_SAMPLING_WINDOW_SIZE : config.samplingWindowSize;
//        this.samplingDepth = config.samplingDepth == = undefined ? DEFAULT_SAMPLING_DEPTH : config.samplingDepth;
//        if (this.samplingDepth > MAXIMUM_SAMPLING_DEPTH) {
//            igv.log("Warning: attempt to set sampling depth > maximum value of 2500");
//            this.samplingDepth = MAXIMUM_SAMPLING_DEPTH;
//        }
//
//        if (config.viewAsPairs) {
//            this.pairsSupported = true;
//        } else {
//            this.pairsSupported = config.pairsSupported == = undefined ? true : config.pairsSupported;
//        }


        try {
            bamIndex = BAMIndex.loadIndex(this.indexPath, null);
            readHeader();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }


    @Override
    public void close() throws IOException {

    }

    @Override
    public List<String> getSequenceNames() throws IOException {
        return Arrays.asList(indexToChr);
    }

    @Override
    public SAMFileHeader getFileHeader() {
        return null;
    }

    @Override
    public Set<String> getPlatforms() {
        return null;
    }

    @Override
    public CloseableIterator<Alignment> iterator() {
        return null;
    }


    @Override
    public boolean hasIndex() {
        return bamIndex != null;
    }


    @Override
    public CloseableIterator<Alignment> query(String chr, int start, int end, boolean contained) throws IOException {
        return new CIterator(chr, start, end);
    }


    class CIterator implements CloseableIterator<Alignment> {

        String chr;
        Integer chrId;
        int start;
        int end;
        Iterator<BAMIndex.Chunk> chunks;
        Iterator<Alignment> currentChunkAlignments;
        Alignment nextAlignment;

        public CIterator(String chr, int start, int end) {
            this.chr = chr;
            this.start = start;
            this.end = end;
            init();
        }

        private void init() {

            Map<String, Integer> chrToIndex = BAMReader.this.chrToIndex;

            chrId = chrToIndex.get(chr);

            if (chrId == null) {
                // TODO
            } else {
                chunks = bamIndex.chunksForRange(chrId, start, end).iterator();
            }

            advance();

        }

        @Override
        public void close() {

        }

        @Override
        public boolean hasNext() {
            return nextAlignment != null;
        }

        @Override
        public Alignment next() {
            Alignment tmp = nextAlignment;
            advance();
            return tmp;
        }

        void advance() {
            nextAlignment = null;
            try {
                while (nextAlignment == null) {
                    if (currentChunkAlignments != null && currentChunkAlignments.hasNext()) {
                        nextAlignment = currentChunkAlignments.next();
                    } else {
                        if (chunks.hasNext()) {
                            BAMIndex.Chunk c = chunks.next();
                            currentChunkAlignments = readAlignments(c, chrId, start, end).iterator();
                        }
                        else {
                            break;
                        }
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
                nextAlignment = null;
            }
        }
    }

    public List<Alignment> readAlignments(BAMIndex.Chunk c, int chrId, int start, int end) throws IOException {

        List<Alignment> alignmentContainer = new ArrayList<>(10000);

        long fetchMin = c.start.block;
        long fetchMax = c.end.block + 65000; // Make sure we get the whole block.

        SeekableStream ss = IGVSeekableStreamFactory.getInstance().getStreamFor(this.path);

        ss.seek(fetchMin);

        byte[] buffer = new byte[(int) (fetchMax - fetchMin + 1)];

        try {
            ss.readFully(buffer);
        } catch (EOFException e) {
            // Can happen with small files
        }

        byte[] unc = BGUnzip.blockUnzip(buffer);

        decodeBamRecords(unc, c.start.offset, alignmentContainer, start, end, chrId); //, self.filter);


        return alignmentContainer;
    }


    public List<Alignment> readAlignments(String chr, int bpStart, int bpEnd) throws IOException {

        List<Alignment> alignmentContainer = new ArrayList<>(10000);

        if (chrToIndex == null) {
            readHeader();
        }

        Map<String, Integer> chrToIndex = this.chrToIndex;

        Integer chrId = chrToIndex.get(chr);

        if (chrId == null) {
            return alignmentContainer;
        } else {

            List<BAMIndex.Chunk> chunks = bamIndex.chunksForRange(chrId, bpStart, bpEnd);

            if (chunks == null || chunks.isEmpty()) {
                return alignmentContainer;
            }

            for (BAMIndex.Chunk c : chunks) {

                long fetchMin = c.start.block;
                long fetchMax = c.end.block + 65000; // Make sure we get the whole block.

                SeekableStream ss = IGVSeekableStreamFactory.getInstance().getStreamFor(this.path);

                ss.seek(fetchMin);

                byte[] buffer = new byte[(int) (fetchMax - fetchMin + 1)];

                try {
                    ss.readFully(buffer);
                } catch (EOFException e) {
                    // Can happen with small files
                }

                byte[] unc = BGUnzip.blockUnzip(buffer);

                decodeBamRecords(unc, c.start.offset, alignmentContainer, bpStart, bpEnd, chrId); //, self.filter);


            }
            return alignmentContainer;
        }
    }

    void decodeBamRecords(byte[] ba, int offset, List<Alignment> alignmentContainer, int min, int max, int chrId) {  //, filter){


        while (true) {

            if (offset >= ba.length) {
                return;
            }

            int blockSize = readInt(ba, offset);
            int blockEnd = offset + blockSize + 4;

            if (blockEnd > ba.length) {
                return;
            }

            int refID = readInt(ba, offset + 4);
            int pos = readInt(ba, offset + 8);

            if (refID < 0) {
                return;   // unmapped reads
            } else if (refID > chrId || pos > max) {
                return;    // off right edge, we're done
            } else if (refID < chrId) {
                continue;   // to left of start, not sure this is possible
            }

            int bmn = readInt(ba, offset + 12);
            int bin = (bmn & 0xffff0000) >> 16;
            int mq = (bmn & 0xff00) >> 8;
            int nl = bmn & 0xff;

            int flag_nc = readInt(ba, offset + 16);
            int flag = (flag_nc & 0xffff0000) >> 16;
            int nc = flag_nc & 0xffff;


            int lseq = readInt(ba, offset + 20);

            int mateRefID = readInt(ba, offset + 24);
            int matePos = readInt(ba, offset + 28);


            byte[] readNameBytes = Arrays.copyOfRange(ba, offset + 36, offset + 36 + nl);
            String readName = new String(readNameBytes);

            int p = offset + 36 + nl;
            byte[] cigarBytes = Arrays.copyOfRange(ba, p, p + 4 * nc);
            p += 4 * nc;

            CigarOperator[] cigarArray = new CigarOperator[nc];
            int lengthOnRef = 0;
            for (int c = 0; c < nc; ++c) {
                int cigop = readInt(cigarBytes, c);
                int opLen = (cigop >> 4);
                char opLtr = CIGAR_DECODER[cigop & 0xf];
                if (opLtr == 'M' || opLtr == 'X' || opLtr == 'D' || opLtr == 'N' || opLtr == '=')
                    lengthOnRef += opLen;
                cigarArray[c] = new CigarOperator(opLen, opLtr);
            }


            if (pos + lengthOnRef < min) {
                offset = blockEnd;
                continue;  // Record out-of-range "to the left", skip to next one
            }

            int seqSize = (lseq + 1) >> 1;
            byte[] seqBytes = Arrays.copyOfRange(ba, p, p + seqSize);
            byte[] sequence = new byte[lseq];

            for (int j = 0; j < seqSize; ++j) {
                byte sb = seqBytes[j];
                sequence[2 * j] = (byte) (SECRET_DECODER[(sb & 0xf0) >> 4]);
                if ((2 * j + 1) < lseq) sequence[2 * j + 1] += SECRET_DECODER[(sb & 0x0f)];
            }


            p += seqSize;


            byte[] qualities;
            if (lseq == 1 && sequence[0] == '*') {
                qualities = new byte[]{Byte.MAX_VALUE}; // TODO == how to represent this?
            } else {
                qualities = Arrays.copyOfRange(ba, p, p + lseq);
            }
            p += lseq;


            ReadMate mate = null;
            boolean isPaired = (flag & PAIRED_FLAG) != 0;
            if (isPaired) {
                boolean mateIsMapped = (flag & MATE_IS_MAPPED_FLAG) != 0;
                String mateChr = mateRefID > 0 ? this.indexToChr[mateRefID] : "";
                boolean mateIsNegativeStrand = ((flag & MATE_STRAND_FLAG) != 0);
                mate = new ReadMate(mateChr, matePos, mateIsNegativeStrand, mateIsMapped);
            }


            byte[] tagBytes = Arrays.copyOfRange(ba, p, blockEnd);


            if (pos + lengthOnRef >= min && pos <= max) {   // && pass filter

                BAMAlignment alignment = new BAMAlignment();
                alignment.start = pos;
                alignment.flags = flag;
                alignment.fragmentLength = readInt(ba, offset + 32);
                alignment.cigarBytes = cigarBytes;
                alignment.lengthOnRef = lengthOnRef;
                alignment.sequence = sequence;
                alignment.mq = mq;
                alignment.readName = readName;
                alignment.chr = this.indexToChr[refID];
                alignment.qualities = qualities;
                alignment.mate = mate;
                alignment.tagBytes = tagBytes;   // Decode these on demand

                makeBlocks(cigarArray, alignment);

                alignmentContainer.add(alignment);
            }

            offset = blockEnd;
        }
    }


    void readHeader() throws IOException {


        int len = (int) bamIndex.firstAlignmentBlock + MAX_GZIP_BLOCK_SIZE;   // Insure we get the complete compressed block containing the header

        //LittleEndianInputStream parser = new LittleEndianInputStream(new BufferedInputStream(ParsingUtils.openInputStream(indexURL)));
        SeekableStream ss = IGVSeekableStreamFactory.getInstance().getStreamFor(this.path);

        byte[] buffer = new byte[len];
        try {
            ss.readFully(buffer);
        } catch (EOFException e) {
            // This can happen with small files
        }

        byte[] uncba = BGUnzip.blockUnzip(buffer);

        int magic = readInt(uncba, 0);
        int samHeaderLen = readInt(uncba, 4);

        String samHeader = new String(uncba, 8, samHeaderLen);

        int nRef = readInt(uncba, samHeaderLen + 8);
        int p = samHeaderLen + 12;

        Map<String, Integer> chrToIndex = new HashMap<>();
        String[] indexToChr = new String[nRef];

        for (int i = 0; i < nRef; ++i) {
            int lName = readInt(uncba, p);
            String name = new String(uncba, p + 4, lName - 1);
            int lRef = readInt(uncba, p + lName + 4);

            if (genome != null) {
                name = genome.getCanonicalChrName(name);
            }

            chrToIndex.put(name, i);
            indexToChr[i] = name;

            p = p + 8 + lName;
        }

        this.chrToIndex = chrToIndex;
        this.indexToChr = indexToChr;


    }

    /**
     * Split the alignment record into blocks as specified in the cigarArray.  Each aligned block contains
     * its portion of the read sequence and base quality strings.  A read sequence or base quality string
     * of "*" indicates the value is not recorded.  In all other cases the length of the block sequence (block.seq)
     * and quality string (block.qual) must == the block length.
     * <p>
     * NOTE: Insertions are not yet treated // TODO
     *
     * @returns array of blocks
     */
    void makeBlocks(CigarOperator[] cigarArray, BAMAlignment alignment) {

        List<AlignmentBlock> blocks = new ArrayList<>();
        List<AlignmentBlock> insertions = null;
        int seqOffset = 0;
        int pos = alignment.start;
        int len = cigarArray.length;
        byte[] blockSeq;
        char gapType;
        byte[] blockQuals;
//                    gapType,
//                    minQ = 5,  //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MIN)
//                    maxQ = 20; //prefs.getAsInt(PreferenceManager.SAM_BASE_QUALITY_MAX)

        for (int i = 0; i < len; i++) {

            CigarOperator c = cigarArray[i];

            switch (c.letter) {
                case 'H':
                    break; // ignore hard clips
                case 'P':
                    break; // ignore pads
                case 'S':
                    seqOffset += c.length;
                    gapType = 'S';
                    break; // soft clip read bases
                case 'N':
                    pos += c.length;
                    gapType = 'N';
                    break;  // reference skip
                case 'D':
                    pos += c.length;
                    gapType = 'D';
                    break;
                case 'I':
                    blockSeq = alignment.sequence.length == 1 ? alignment.sequence : Arrays.copyOfRange(alignment.sequence, seqOffset, seqOffset + c.length);
                    blockQuals = alignment.qualities == null ? null : Arrays.copyOfRange(alignment.qualities, seqOffset, seqOffset + c.length);
                    if (insertions == null) {
                        insertions = new ArrayList<>();
                    }
                    insertions.add(new AlignmentBlockImpl(pos, blockSeq, blockQuals));
                    seqOffset += c.length;
                    break;
                case 'M':
                case '=':
                case 'X':
                    blockSeq = alignment.sequence.length == 1 ? alignment.sequence : Arrays.copyOfRange(alignment.sequence, seqOffset, seqOffset + c.length);
                    blockQuals = alignment.qualities == null ? null : Arrays.copyOfRange(alignment.qualities, seqOffset, seqOffset + c.length);
                    blocks.add(new AlignmentBlockImpl(pos, blockSeq, blockQuals));
                    seqOffset += c.length;
                    pos += c.length;
                    break;

                default:
                    log.error("Error processing cigar element: " + c.length + c.letter);
            }
        }

        alignment.alignmentBlocks = blocks.toArray(new AlignmentBlockImpl[]{});
        if (insertions != null) {
            alignment.insertions = insertions.toArray(new AlignmentBlockImpl[]{});
        }

    }


    public String readString(ByteBuffer ba) throws IOException {
        ByteArrayOutputStream bis = new ByteArrayOutputStream(100);

        byte b;
        while ((b = ba.get()) != 0) {
            if (b < 0) {
                throw new EOFException();
            }

            bis.write(b);
        }

        return new String(bis.toByteArray());
    }


    public static int readInt(byte[] ba, int offset) {
        return (ba[offset + 3] << 24) + (ba[offset + 2] << 24 >>> 8) + (ba[offset + 1] << 24 >>> 16) + (ba[offset] << 24 >>> 24);
    }

}


