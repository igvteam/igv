package org.broad.igv.feature.genome.fasta;

import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broad.igv.util.LittleEndianInputStream;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

/**
 * Created by jrobinso on 6/23/17.
 */
public class FastaBlockCompressedSequence extends FastaIndexedSequence {

    Mapping[] gziMappings;
    Mapping zeroMapping = new Mapping(0, 0);

    public FastaBlockCompressedSequence(String path) throws IOException {

        this(path, null);
    }

    public FastaBlockCompressedSequence(String path, String indexPath) throws IOException {

        super(path);

        if(indexPath == null) indexPath = path + ".gzi";

        readGziMappings(indexPath);
    }

    @Override
    /**
     * Read the bytes between VIRTUAL file position posStart and posEnd
     *
     * @throws IOException
     */
    protected byte[] readBytes(long posStart, long posEnd) throws IOException {

        Mapping m1 = findBlockContaining(posStart);
        int d1 = (int) (posStart - m1.uncompressedOffset);
        long vp1 = m1.compressedOffset << 16 | d1;

        SeekableStream ss = null;
        try {
            int nBytes = (int) (posEnd - posStart);

            int bufferSize = Math.max(512000, nBytes / 8);

            ss = new SeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(path), bufferSize);

            BlockCompressedInputStream bis = new BlockCompressedInputStream(ss);

            byte[] bytes = new byte[nBytes];

            bis.seek(vp1);
            readFully(bytes, bis);

            return bytes;
        } finally {
            if (ss != null) {
                ss.close();
            }
        }
    }

    protected Mapping findBlockContaining(long uoffset) {

        int ilo = 0, ihi = gziMappings.length - 1;
        while (ilo <= ihi) {
            int i = (ilo + ihi) / 2;
            Mapping mapping = gziMappings[i];
            if (uoffset < mapping.uncompressedOffset) ihi = i - 1;
            else if (uoffset >= mapping.uncompressedOffset) ilo = i + 1;
            else break;
        }

        return ilo == 0 ? zeroMapping : gziMappings[ilo - 1];
    }


    private void readGziMappings(String gziPath) throws IOException {

        InputStream is = null;
        LittleEndianInputStream dis = null;

        is = ParsingUtils.openInputStream(gziPath);
        dis = new LittleEndianInputStream(new BufferedInputStream(is));

        int nEntries = (int) dis.readLong();

        gziMappings = new Mapping[nEntries];

        for (int i = 0; i < nEntries; i++) {
            gziMappings[i] = new Mapping(dis.readLong(), dis.readLong());
        }

    }

    private void readFully(byte[] b, InputStream is) throws IOException {
        int len = b.length;
        if (len < 0) {
            throw new IndexOutOfBoundsException();
        } else {
            int count;
            for (int n = 0; n < len; n += count) {
                count = is.read(b, n, len - n);
                if (count < 0) {
                    throw new EOFException();
                }
            }

        }
    }

    public static class Mapping {

        long compressedOffset;
        long uncompressedOffset;

        public Mapping(long compressedOffset, long uncompressedOffset) {
            this.compressedOffset = compressedOffset;
            this.uncompressedOffset = uncompressedOffset;
        }
    }


}
