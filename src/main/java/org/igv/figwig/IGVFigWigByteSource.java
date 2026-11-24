package org.igv.figwig;

import org.figwig.FigWigFile;
import org.figwig.FigWigByteSource;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.igv.util.FileUtils;
import org.igv.util.stream.IGVSeekableStreamFactory;

import java.io.EOFException;
import java.io.File;
import java.io.IOException;

/**
 * A figwig file read through IGV's streams, which carry IGV's authentication, signing and redirects:
 * http(s), s3 and gs paths, and any other path that is not a local file.
 */
public class IGVFigWigByteSource implements FigWigByteSource {

    private static final boolean POSITIONAL = !Boolean.getBoolean("figwig.noPositional");

    private final String path;
    private final SeekableStream stream;
    private final boolean remote;

    /** Open a figwig file as IGV reads it: a local file by the reader itself, anything else through IGV's streams. */
    public static FigWigFile open(String path) throws IOException {
        if (POSITIONAL && new File(path).isFile()) {
            return FigWigFile.open(path);
        }
        return FigWigFile.open(path, new IGVFigWigByteSource(path));
    }

    public IGVFigWigByteSource(String path) throws IOException {
        this.path = path;
        this.stream = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
        this.remote = FileUtils.isRemote(path) && !path.startsWith("ftp:");
    }

    @Override
    public void read(long offset, byte[] out) throws IOException {
        if (remote) {
            // One range request per read on a stream of its own, so IGV's loader threads read concurrently.
            try (SeekableStream s = IGVSeekableStreamFactory.getInstance().getStreamFor(path)) {
                fill(s, offset, out);
            }
        } else {
            synchronized (stream) {
                fill(stream, offset, out);
            }
        }
    }

    private void fill(SeekableStream s, long offset, byte[] out) throws IOException {
        s.seek(offset);
        try {
            s.readFully(out);
        } catch (EOFException e) {
            throw new IOException("figwig: " + path + " ended before " + (offset + out.length));
        }
    }

    @Override
    public long length() {
        return stream.length();
    }

    @Override
    public boolean isRemote() {
        return remote;
    }

    @Override
    public void close() throws IOException {
        stream.close();
    }
}
