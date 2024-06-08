package org.broad.igv.track;

import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

public class FileFormatUtils {

    static final byte[] BAM_MAGIC = "BAM\1".getBytes();
    static final byte[] CRAM_MAGIC = "CRAM".getBytes();

    static final long BIGWIG_MAGIC = 2291137574l; // BigWig Magic
    static final long BIGBED_MAGIC = 2273964779l; // BigBed Magic


    public static boolean isGzip(SeekableStream seekableStream) throws IOException {
        seekableStream.seek(0);
        byte[] bytes = new byte[2];
        seekableStream.readFully(bytes);
        return bytes[0] == 31 && bytes[1] == -117;
    }

    public static boolean isBAM(SeekableStream seekableStream) throws IOException {
        if (isGzip(seekableStream)) {
            seekableStream.seek(0);
            BlockCompressedInputStream blockCompressedInputStream = new BlockCompressedInputStream(seekableStream);
            byte[] bytes = new byte[4];
            blockCompressedInputStream.read(bytes);
            return Arrays.equals(bytes, BAM_MAGIC);
        } else {
            return false;
        }
    }

    public static boolean isBAM(String path) throws IOException {
        SeekableStream seekableStream = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        if (isGzip(seekableStream)) {
            seekableStream.seek(0);
            BlockCompressedInputStream blockCompressedInputStream = new BlockCompressedInputStream(new BufferedInputStream(seekableStream));
            byte[] bytes = new byte[4];
            blockCompressedInputStream.read(bytes);
            return Arrays.equals(bytes, BAM_MAGIC);
        } else {
            return false;
        }
    }

    public static String determineFormat(String path) throws IOException {

        byte[] bytes = new byte[1000];

        SeekableStream seekableStream = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        if (isGzip(seekableStream)) {
            seekableStream.seek(0);
            GZIPInputStream inputStream = new GZIPInputStream(seekableStream);
            inputStream.read(bytes);
        } else {
            seekableStream.seek(0);
            seekableStream.readFully(bytes);
        }

        // Try BAM and CRAM
        if (Arrays.equals(bytes, 0, 4, BAM_MAGIC, 0, 4)) {
            return "bam";
        }
        if (Arrays.equals(bytes, 0, 4, CRAM_MAGIC, 0, 4)) {
            return "cram";
        }

        // BIGWIG - BIGBED
        UnsignedByteBuffer byteBuffer = UnsignedByteBuffer.wrap(bytes);
        long m = byteBuffer.getUInt();
        if (m == BIGWIG_MAGIC) {
            return "bigwig";
        }
        if (m == BIGBED_MAGIC) {
            return "bigbed";
        }

        String magicString = new String(Arrays.copyOfRange(bytes, 0, 4));
        if ((magicString.startsWith("TDF") || magicString.startsWith("IBF"))) {
            return "tdf";
        }

        // Knowable text formats
        try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(bytes)));
            String firstLine = reader.readLine();
            if (firstLine.startsWith("##fileformat=VCF")) {
                return "vcf";
            }
            if (firstLine.startsWith("##gff-version 3")) {
                return "gff3";
            }
            if (firstLine.startsWith("##gff-version")) {
                return "gff";
            }
        } catch (IOException e) {
            // Ignore, apparently not text
        }


        // Format unknown
        return null;
    }
}
