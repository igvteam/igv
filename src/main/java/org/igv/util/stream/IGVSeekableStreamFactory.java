package org.igv.util.stream;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.igv.util.FileUtils;
import org.igv.util.HttpUtils;

import java.io.File;
import java.io.IOException;
import java.net.URL;

/**
 * @author Jim Robinson
 */
public class IGVSeekableStreamFactory implements ISeekableStreamFactory {

    private static IGVSeekableStreamFactory instance;

    static {
        instance = new IGVSeekableStreamFactory();
    }

    private IGVSeekableStreamFactory() {
    }

    public static IGVSeekableStreamFactory getInstance() {
        return instance;
    }

    public SeekableStream getStreamFor(URL url) throws IOException {
        return getStreamFor(url.toExternalForm());
    }

    public SeekableStream getStreamFor(String path) throws IOException {

        if (path.endsWith(".list")) {
            return new SeekableSplitStream(path);
        } else {
            path = mapPath(path);
            SeekableStream is = null;
            if (FileUtils.isRemote(path)) {
                final URL url = HttpUtils.createURL(path);
                is = new IGVSeekableHTTPStream(url);
            } else if (path.toLowerCase().startsWith("ftp:")) {
                final URL url = HttpUtils.createURL(path);
                is = new IGVSeekableFTPStream(url);
            } else {
                is = new SeekableFileStream(new File(path));
            }
            return is;
        }
    }

    public SeekableStream getBufferedStream(SeekableStream stream) {
        return getBufferedStream(stream, IGVSeekableBufferedStream.DEFAULT_BUFFER_SIZE);
    }

    public SeekableStream getBufferedStream(SeekableStream stream, int bufferSize) {
        return new IGVSeekableBufferedStream(stream, bufferSize);
    }

    private String mapPath(String path) {
        if (path.startsWith("ftp://ftp.ncbi.nlm.nih.gov/geo")) {
            return path.replace("ftp://ftp.ncbi.nlm.nih.gov/geo", "https://ftp.ncbi.nlm.nih.gov/geo");
        } else {
            return path;
        }
    }

}
