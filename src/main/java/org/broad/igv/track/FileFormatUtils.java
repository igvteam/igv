package org.broad.igv.track;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;
import org.broad.igv.ucsc.twobit.UnsignedByteBufferImpl;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.*;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

public class FileFormatUtils {

    static final byte[] BAM_MAGIC = "BAM\1".getBytes();
    static final byte[] CRAM_MAGIC = "CRAM".getBytes();
    static final long BIGWIG_MAGIC = 2291137574l; // BigWig Magic
    static final long BIGBED_MAGIC = 2273964779l; // BigBed Magic

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
            try {
                seekableStream.readFully(bytes);
            } catch (EOFException e) {
                // File is < bytes.length, this is o.k. for this application, but trim zeroes
                int idx = 0;
                while (idx < bytes.length && bytes[idx] != 0) {
                    idx++;
                }
                bytes = Arrays.copyOf(bytes, idx);
            }
        }

        // Try BAM and CRAM
        if (Arrays.equals(bytes, 0, 4, BAM_MAGIC, 0, 4)) {
            return "bam";
        }
        if (Arrays.equals(bytes, 0, 4, CRAM_MAGIC, 0, 4)) {
            return "cram";
        }

        // BIGWIG - BIGBED
        UnsignedByteBuffer byteBuffer = UnsignedByteBufferImpl.wrap(bytes);
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
            if (maybeSampleInfo(bytes)) {
                return "sampleinfo";
            }
        } catch (IOException e) {
            // Ignore, apparently not text
        }


        // Format unknown
        return null;
    }

    /**
     * Minimal validation of a sample info file (1) contents are UTF-8, (2) file contains tabs
     *
     * @param bytes
     * @return
     */
    private static boolean maybeSampleInfo(byte[] bytes) {

        try {
            if(!isUTF8(bytes)) {
                return false;
            }

            String converted = new String(bytes, "UTF-8");

            // Now look for tab characters, there needs to be at least 2
            int firstTab = converted.indexOf('\t');
            if (firstTab < 0) return false;
            int secondTab = converted.indexOf('\t', firstTab + 1);
            if (secondTab < 0) return false;

            //
            boolean nullSeen = false;
            for(int i=0; i<converted.length(); i++) {
                if(converted.charAt(i) == 0) {
                    nullSeen = true;
                } else {
                    if(nullSeen) return false;
                }
            }
        } catch (UnsupportedEncodingException e) {
            return false;
        }
        return true;

        //
    }

    private static boolean isGzip(SeekableStream seekableStream) throws IOException {
        seekableStream.seek(0);
        byte[] bytes = new byte[2];
        seekableStream.readFully(bytes);
        return bytes[0] == 31 && bytes[1] == -117;
    }

    private static boolean isUTF8(byte[] data) {


        int len = data.length;

        // for each value in the data array we have to take
        // the "and" with the masks array
        for (int i = 0; i < len; i++) {
            int curr = data[i];

            // method to check the array if the
            // and with the num and masks array is
            // 0 then return true
            int type = getType(curr);

            if (type == 0) {
                continue;
            } else if (type > 1 && i + type <= len) {
                while (type-- > 1) {
                    if (getType(data[++i]) != 1) {
                        return false;
                    }
                }
            } else {
                return false;
            }
        }
        return true;
    }

    // method to check the type
    private static int getType(int num) {
        int[] masks = {128, 64, 32, 16, 8};
        for (int i = 0; i < 5; i++) {

            // checking the each input
            if ((masks[i] & num) == 0) {
                return i;
            }
        }

        return -1;
    }
}

