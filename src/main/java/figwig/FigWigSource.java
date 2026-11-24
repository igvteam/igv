package figwig;

import java.io.Closeable;
import java.io.IOException;

/**
 * The bytes of a figwig file the reader does not open itself: anything but a local file (http(s), s3
 * and gs paths, or a local file read some other way).  The application supplies it to
 * {@link FigWigFile#open(String, FigWigSource)}, and with it its own authentication, signing and redirects.
 */
public interface FigWigSource extends Closeable {

    /**
     * Fill {@code out} with the bytes from {@code offset}, failing when the file ends first.  The reader
     * calls this from several threads at once.
     */
    void read(long offset, byte[] out) throws IOException;

    /** The file's length, or 0 when it is not known. */
    long length() throws IOException;

    /**
     * Whether each read is a round trip to a server.  The reader then asks for longer stretches at a
     * time and sends the requests of one load at once.
     */
    boolean isRemote();
}
