package org.broad.igv.feature.genome;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.LittleEndianInputStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
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

    public TwoBitSequence(String path) throws IOException {
        this.path = path;
        init();

    }


    private void init() throws IOException {

        SeekableStream is = null;

        is = IGVSeekableStreamFactory.getInstance().getStreamFor(path);

        SeekableStream bis = IGVSeekableStreamFactory.getInstance().getBufferedStream(is, 1000);

        LittleEndianInputStream lis = new LittleEndianInputStream(bis);

        int signature = lis.readInt();

        boolean littleEndian = SIGNATURE_LE == signature;
        // TODO -- handle big endian

        int version = lis.readInt();   // Should be zero

        int seqCount = lis.readInt();

        int reserved = lis.readInt();    // Should be zero

        Map<String, Integer> offsets = new LinkedHashMap<>();

        for (int i = 0; i < seqCount; i++) {

            byte nameSize = lis.readByte();

            byte[] seqNameBytes = new byte[nameSize];
            lis.readFully(seqNameBytes);
            String seqName = new String(seqNameBytes);

            int offset = lis.readInt();

            offsets.put(seqName, offset);

        }

        for (Integer offset : offsets.values()) {

            bis.seek(offset);
            int dnaSize = lis.readInt();

            int nBlockCount = lis.readInt();
            int [] nBlockStarts = new int[nBlockCount];
            for(int i=0; i<nBlockCount; i++) {
                nBlockStarts[i] = lis.readInt();
            }
            int [] nBlockSizes = new int[nBlockCount];
            for(int i=0; i<nBlockCount; i++) {
                nBlockSizes[i] = lis.readInt();
            }

            int maskBlockCount = lis.readInt();
            int [] maskBlockStarts = new int[maskBlockCount];
            for(int i=0; i<maskBlockCount; i++) {
                maskBlockStarts[i] = lis.readInt();
            }
            int [] maskBlockSizes = new int[maskBlockCount];
            for(int i=0; i<maskBlockCount; i++) {
                maskBlockSizes[i] = lis.readInt();
            }

        }


    }

    @Override
    public byte[] getSequence(String chr, int start, int end, boolean useCache) {
        return new byte[0];
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
