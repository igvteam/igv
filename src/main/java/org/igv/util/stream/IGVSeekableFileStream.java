package org.igv.util.stream;

import htsjdk.samtools.seekablestream.SeekableStream;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * Seekable stream class for a local file.  This implementaion differs from the htsjdk by treating the
 * random access f
 */
public class IGVSeekableFileStream extends SeekableStream {


    File file;
    long position;
   // RandomAccessFile fis;

    public IGVSeekableFileStream(final File file) throws FileNotFoundException {
        this.file = file;
    }

    @Override
    public long length() {
        return file.length();
    }

    @Override
    public boolean eof() throws IOException {
        return position == file.length();
    }

    @Override
    public void seek(final long position) throws IOException {
        this.position = position;
    }

    @Override
    public long position() throws IOException {
        return position;
    }

    @Override
    public long skip(long n) throws IOException {
        long initPos = position;
        position = Math.min(file.length(), position + n);
        return position - initPos;
    }

    @Override
    public int read(final byte[] buffer, final int offset, final int length) throws IOException {

        if (length < 0) {
            throw new IndexOutOfBoundsException();
        }

        RandomAccessFile fis = null;
        try {
            fis = new RandomAccessFile(file, "r");
            fis.seek(position);
            int n = 0;
            while (n < length) {
                final int count = fis.read(buffer, offset + n, length - n);
                if (count < 0) {
                    if (n > 0) {
                        return n;
                    } else {
                        return count;
                    }
                }
                n += count;
            }
            return n;
        } finally  {
            fis.close();
        }
    }

    @Override
    public int read() throws IOException {
        RandomAccessFile fis = null;
        fis.seek(position);
        try {
            return fis.read();
        } finally  {
            fis.close();
        }
    }

    @Override
    public int read(byte[] b) throws IOException {
        RandomAccessFile fis = null;
        fis.seek(position);
        try {
        return fis.read(b);
        } finally  {
            fis.close();
        }
    }

    @Override
    public String getSource() {
        return file.getAbsolutePath();
    }


    @Override
    public void close() throws IOException {

    }

}
