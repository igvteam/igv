package org.broad.igv.session;

import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Class to parse a UCSC session file
 *
 * @author Jim Robinson
 * @date 1/12/12
 */
public class UCSCSessionReader implements SessionReader {

    IGV igv;

    public UCSCSessionReader(IGV igv) {
        this.igv = igv;
    }

    /**
     * TODO -- load asynchronously when possible
     *
     * @param inputStream
     * @param session
     * @param sessionName
     * @throws IOException
     */
    public void loadSession(InputStream inputStream, Session session, String sessionName) throws IOException {

        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

        String trackLine = null;
        String nextLine;

        LinkedHashMap<String, List<Track>> loadedTracks = new LinkedHashMap();

        while ((nextLine = reader.readLine()) != null) {
            if (nextLine.startsWith("#")) {
                continue;
            } else if (nextLine.startsWith("browser")) {
                parseBrowserLine(nextLine, session);
            } else if (nextLine.startsWith("track")) {
                trackLine = nextLine;
                TrackProperties props = new TrackProperties();
                ParsingUtils.parseTrackLine(nextLine, props);
                String dataURL = props.getDataURL();
                if (dataURL != null) {
                    ResourceLocator locator = new ResourceLocator(dataURL);
                    locator.setTrackLine(trackLine);
                    trackLine = null;
                    loadedTracks.put(dataURL, igv.load(locator));

                }
            } else {
                ResourceLocator locator = new ResourceLocator(nextLine);
                locator.setTrackLine(trackLine);
                trackLine = null;
                loadedTracks.put(nextLine, igv.load(locator));
            }

        }

        for (Map.Entry<String, List<Track>> entry : loadedTracks.entrySet()) {
            String path = entry.getKey();
            //TrackPanel panel = IGV.getInstance().getPanelFor(new ResourceLocator(path));
            // If loading from UCSC use a single panel
            TrackPanel panel = igv.getTrackPanel(IGV.DATA_PANEL_NAME);
            panel.addTracks(entry.getValue());
        }

    }

    private boolean isAlignmentFile(String path) {
        return path.endsWith(".bam") || path.endsWith(".entries") || path.endsWith(".sam");
    }


    /**
     * Browser lines have 2 or 3 tokens,  e.g.
     * <p/>
     * browser position chr19:11500000-12000000
     *
     * @param line
     */
    private void parseBrowserLine(String line, Session session) {

        String[] tokens = line.split("\\s+");
        if (tokens.length >= 3 && tokens[1].equals("position")) {
            session.setLocus(tokens[2]);
        }

    }
}
