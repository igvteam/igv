package org.broad.igv.feature.genome.sequence;

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.FileUtils;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;


public class SequenceFactory {

    private static Logger log = LogManager.getLogger(SequenceFactory.class);

    public static Sequence getSequence(String path) throws IOException {
        return getSequence(path, null, null);
    }


    public static Sequence getSequence(String path, String indexPath) throws IOException {
        return getSequence(path, indexPath, null);
    }


        /**
         * Return the Sequence object for the associated paths.
         * <p>
         * TODO -- this logic relieds on file extensions, won't neccessarily work for signed or DRS urls
         *
         * @param path
         * @param indexPath
         * @param gziIndexPath
         * @return
         * @throws IOException
         */
    public static Sequence getSequence(String path, String indexPath, String gziIndexPath) throws IOException {

        String p;
        if (FileUtils.isRemote(path)) {
            try {
                p = ((new URL(path)).getPath());
            } catch (MalformedURLException e) {
                log.error("Error parsing url: " + path);
                p = path;
            }
        } else {
            p = path;
        }

        if (gziIndexPath != null || p.endsWith(".gz")) {
            return new FastaBlockCompressedSequence(path, gziIndexPath, indexPath);
        } else if (p.endsWith(".2bit")) {
            return new TwoBitSequence(path);
        } else {
            return new FastaIndexedSequence(path, indexPath);
        }
    }

}
