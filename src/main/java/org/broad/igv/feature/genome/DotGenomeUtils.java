package org.broad.igv.feature.genome;

import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.genome.load.GenomeDescriptor;
import org.broad.igv.feature.genome.load.GenomeLoader;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.download.Downloader;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.Utilities;

import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLDecoder;
import java.util.Map;

/**
 * Utilities for the ".genome" format.  These files are not created anymore, having been replaced by the genome json
 * format, but we need to maintain support for reading and managing them.
 */
public class DotGenomeUtils {

    private static Logger log = LogManager.getLogger(DotGenomeUtils.class);


    /**
     * Returns a File of the provided genomePath. If the genomePath is a URL, it will be downloaded
     * and saved in the genome cache directory.
     *
     * @param genomePath
     * @return
     * @throws MalformedURLException
     * @throws UnsupportedEncodingException
     */
    public static File getDotGenomeFile(String genomePath) throws MalformedURLException, UnsupportedEncodingException {

        File archiveFile;

        if (HttpUtils.isRemoteURL(genomePath.toLowerCase())) {
            // We need a local copy, as there is no http zip file reader
            URL genomeArchiveURL = HttpUtils.createURL(genomePath);
            final String tmp = URLDecoder.decode(genomeArchiveURL.getFile(), "UTF-8");
            String cachedFilename = Utilities.getFileNameFromURL(tmp);
            if (!DirectoryManager.getGenomeCacheDirectory().exists()) {
                DirectoryManager.getGenomeCacheDirectory().mkdir();
            }
            archiveFile = new File(DirectoryManager.getGenomeCacheDirectory(), cachedFilename);
            Frame parent = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;
            Downloader.download(genomeArchiveURL, archiveFile, parent);
        } else {
            archiveFile = new File(genomePath);
        }
        return archiveFile;
    }



    public static File getLocalFasta(String id) {
        return GenomeLoader.localSequenceMap.get(id);
    }

    public static void removeLocalFasta(String id) {
        GenomeLoader.localSequenceMap.remove(id);
        updateSequenceMapFile();
    }


    private static void updateSequenceMapFile() {

        PrintWriter pw = null;

        try {
            File sequenceFile = new File(DirectoryManager.getGenomeCacheDirectory(), GenomeDescriptor.SEQUENCE_MAP_FILE);
            pw = new PrintWriter(new BufferedWriter(new FileWriter(sequenceFile)));

            for (Map.Entry<String, File> entry : GenomeLoader.localSequenceMap.entrySet()) {
                pw.println(entry.getKey() + "\t" + entry.getValue());
            }
        } catch (IOException e) {
            log.error("Error writing sequence map", e);
        } finally {
            if (pw != null) pw.close();
        }
    }

}
