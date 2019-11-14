package org.broad.igv.sam.lite;

import htsjdk.samtools.util.StringUtil;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.*;
import org.broad.igv.sam.ReadMate;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * Created by jrobinso on 3/13/17.
 */
public class BAMAlignment extends SAMAlignment {

    static int READ_PAIRED_FLAG = 0x1;
    static int PROPER_PAIR_FLAG = 0x2;
    static int READ_UNMAPPED_FLAG = 0x4;
    static int MATE_UNMAPPED_FLAG = 0x8;
    static int READ_STRAND_FLAG = 0x10;
    static int MATE_STRAND_FLAG = 0x20;
    static int FIRST_OF_PAIR_FLAG = 0x40;
    static int SECOND_OF_PAIR_FLAG = 0x80;
    static int SECONDARY_ALIGNMENT_FLAG = 0x100;
    static int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
    static int DUPLICATE_READ_FLAG = 0x400;
    static int SUPPLEMENTARY_ALIGNMENT_FLAG = 0x800;

    public int flags;
    public int fragmentLength;
    public int lengthOnRef;
    public int start;
    public byte[] sequence;
    public byte[] qualities;
    public int mq;
    public String readName;
    public byte[] tagBytes;
    public byte[] cigarBytes;
    public String chr;
    public ReadMate mate;
    private Map<String, Object> tagDictionary;


    @Override
    public String getReadName() {
        return readName;
    }

    @Override
    public int getMappingQuality() {
        return mq;
    }

    @Override
    public int getInferredInsertSize() {
        return 0;
    }

    @Override
    public String getCigarString() {
        return new String(cigarBytes);
    }

    @Override
    public String getReadLengthString() {
        return String.valueOf(sequence.length);
    }

    @Override
    public String getReadSequence() {
        return new String(sequence);
    }

    @Override
    public boolean isNegativeStrand() {
        return (this.flags & READ_STRAND_FLAG) != 0;
    }

    @Override
    public boolean isFirstOfPair() {
        return (this.flags & FIRST_OF_PAIR_FLAG) != 0;
    }

    @Override
    public boolean isSecondOfPair() {
        return (this.flags & SECOND_OF_PAIR_FLAG) != 0;
    }

    @Override
    public boolean isDuplicate() {
        return (this.flags & DUPLICATE_READ_FLAG) != 0;
    }

    @Override
    public boolean isMapped() {
        return (this.flags & READ_UNMAPPED_FLAG) == 0;
    }

    @Override
    public boolean isPaired() {
        return (this.flags & READ_PAIRED_FLAG) != 0;
    }

    @Override
    public boolean isProperPair() {
        return (this.flags & PROPER_PAIR_FLAG) != 0;
    }

    @Override
    public boolean isSupplementary() {
        return (this.flags & SUPPLEMENTARY_ALIGNMENT_FLAG) != 0;
    }

    @Override
    public boolean isVendorFailedRead() {
        return (this.flags & READ_FAILS_VENDOR_QUALITY_CHECK_FLAG) != 0;
    }

    @Override
    public boolean isPrimary() {
        return (this.flags & SECONDARY_ALIGNMENT_FLAG) == 0;
    }

    @Override
    public int getAlignmentStart() {
        return start;
    }

    @Override
    public int getAlignmentEnd() {
        return start + lengthOnRef;
    }

    @Override
    public Object getAttribute(String key) {
        return null;
    }


    @Override
    public String getSample() {
        return null;
    }

    @Override
    public String getReadGroup() {
        return null;
    }

    @Override
    public String getLibrary() {
        return null;
    }

    @Override
    public String getAttributeString(boolean truncate) {

         IGVPreferences prefMgr = PreferencesManager.getPreferences();  // TODO -- pass hidden tags as argument


        // List of tags to skip.  Some tags, like MD and SA, are both quite verbose and not easily
        // interpreted by a human reader.  It is best to just hide these tags.  The list of tags
        // to hide is set through the SAM_HIDDEN_TAGS preference.
        ArrayList<String> tagsToHide = new ArrayList<String>(),
                tagsHidden = new ArrayList<String>();

        String samHiddenTagsPref = prefMgr.get(Constants.SAM_HIDDEN_TAGS);
        for (String s : (samHiddenTagsPref == null ? "" : samHiddenTagsPref).split("[, ]")) {
            if (!s.equals("")) {
                tagsToHide.add(s);
            }
        }

        StringBuffer buf = new StringBuffer();



        Map<String, Object> attributes = getTagDictionary();

        if (attributes != null && !attributes.isEmpty()) {

            for (Map.Entry<String, Object> entry: attributes.entrySet()) {

                String tag = entry.getKey();
                Object value = entry.getValue();


                if (tagsToHide.contains(tag)) {
                    tagsHidden.add(tag);
                    continue;
                }
                buf.append("<br>" + tag + " = ");

                if (value.getClass().isArray()) { // ignore array types
                    buf.append("[not shown]<br>");
                    continue;
                }

                // Break tag
                final String tagValue = value.toString();
                final int maxLength = 70;
                if (tagValue.length() > maxLength && truncate) {
                    String[] tokens = tagValue.split("<br>");
                    for (String token : tokens) {
                        if (token.length() > maxLength) {
                            // Insert line breaks
                            String remainder = token;
                            while (remainder.length() > maxLength) {
                                String tmp = remainder.substring(0, maxLength);
                                int spaceIndex = tmp.lastIndexOf(' ');
                                int idx = spaceIndex > 30 ? spaceIndex : maxLength;
                                final String substring = remainder.substring(0, idx);
                                buf.append(substring);
                                buf.append("<br>");
                                remainder = remainder.substring(idx);
                            }
                            buf.append(remainder);
                            buf.append("<br>");

                        } else {
                            buf.append(token);
                            buf.append("<br>");
                        }
                    }
                } else {
                    buf.append(tagValue);
                }

            }

            if (tagsHidden.size() > 0) {
                buf.append("<br>Hidden tags: " + String.join(", ", tagsHidden));
            }
        }
        return buf.toString();
    }

    private Map<String, Object> getTagDictionary() {
        if (this.tagDictionary == null) {
            if (this.tagBytes == null) {
                return null;
            } else {
                this.tagDictionary = decodeTags(this.tagBytes);

            }
        }
        return this.tagDictionary;
    }

    /**
      A [!-~] Printable character
     i [-+]?[0-9]+ Signed integer5
     f [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? Single-precision
     oating number
     Z [ !-~]+ Printable string, including space
     H [0-9A-F]+ Byte array in the Hex format6
     B [cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+ Integer or numeric array
     */
    private Map<String, Object> decodeTags(byte[] ba) {

        Map<String, Object> tags = new LinkedHashMap<>();
        ByteBuffer byteBuffer = ByteBuffer.wrap(ba);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        while (byteBuffer.hasRemaining()) {
            int p = byteBuffer.position();
            String tag = new String(ba, p, 2);
            byteBuffer.position(p+2);
            char type = (char) byteBuffer.get();
            Object value;

            switch (type) {
                case 'Z':
                    value = readNullTerminatedString(byteBuffer);
                    break;
                case 'A':
                    value = (char)byteBuffer.get();
                    break;
                case 'I':
                    final long val = byteBuffer.getInt() & 0xffffffffL;
                    if ( val <= Integer.MAX_VALUE ) {
                        value = (int)val;
                    }
                    else {
                        value = val;
                    }
                    break;
                case 'i':
                    value =  byteBuffer.getInt();
                    break;
                case 's':
                    value =  byteBuffer.getShort();
                    break;
                case 'S':
                    // Convert to unsigned short stored in an int
                    value =  (int) byteBuffer.getShort() & 0xffff;
                    break;
                case 'c':
                    value =  byteBuffer.get();
                    break;
                case 'C':
                    // Convert to unsigned byte stored in an int
                    value =  (int)byteBuffer.get() & 0xff;
                    break;
                case 'f':
                    value = byteBuffer.getFloat();
                    break;
                case 'H':
                    final String hexRep = readNullTerminatedString(byteBuffer);
                    value = StringUtil.hexStringToBytes(hexRep);
                    break;
                case 'B':
                    value = readArray(byteBuffer);
                    break;
                default:
                    value = "Unrecognized tag type: " + type;
            }

            tags.put(tag, value);
        }
        return tags;
    }

    private static String readNullTerminatedString(final ByteBuffer byteBuffer) {
        // Count the number of bytes in the string
        byteBuffer.mark();
        final int startPosition = byteBuffer.position();
        while (byteBuffer.get() != 0) {}
        final int endPosition = byteBuffer.position();

        // Don't count null terminator
        final byte[] buf = new byte[endPosition - startPosition - 1];
        // Go back to the start of the string and read out the bytes
        byteBuffer.reset();
        byteBuffer.get(buf);
        // Skip over the null terminator
        byteBuffer.get();
        return StringUtil.bytesToString(buf);
    }

    private static Object readArray(final ByteBuffer byteBuffer) {

        final byte arrayType = byteBuffer.get();
        final boolean isUnsigned = Character.isUpperCase(arrayType);
        final int length = byteBuffer.getInt();
        final Object value;
        switch (Character.toLowerCase(arrayType)) {
            case 'c': {
                final byte[] array = new byte[length];
                value = array;
                byteBuffer.get(array);
                break;
            }
            case 's': {
                final short[] array = new short[length];
                value = array;
                for (int i = 0; i < length; ++i) {
                    array[i] = byteBuffer.getShort();
                }
                break;
            }

            case 'i': {
                final int[] array = new int[length];
                value = array;
                for (int i = 0; i < length; ++i) {
                    array[i] = byteBuffer.getInt();
                }
                break;
            }

            case 'f': {
                final float[] array = new float[length];
                value = array;
                for (int i = 0; i < length; ++i) {
                    array[i] = byteBuffer.getFloat();
                }
                break;
            }

            default:
                throw new RuntimeException("Unrecognized tag array type: " + (char)arrayType);
        }
        return  value;
    }


}
