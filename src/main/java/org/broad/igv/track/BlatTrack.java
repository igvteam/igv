package org.broad.igv.track;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broad.igv.event.DataLoadedEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.blat.BlatQueryWindow;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;


public class BlatTrack extends FeatureTrack {

    private static Logger log = LogManager.getLogger(BlatTrack.class);

    private static Random rand = new Random();

    String sequence;
    String db;
    List<PSLRecord> features;

    /**
     * Restore from session
     */
    public BlatTrack() {
        setDisplayMode(Track.DisplayMode.SQUISHED);
        setColor(Color.DARK_GRAY);
    }

    /**
     * Create new blat track
     * @param sequence
     * @param features
     * @param trackLabel
     */
    public BlatTrack(String db, String sequence, List<PSLRecord> features, String trackLabel) {
        super(null, guid(), trackLabel);
        setDisplayMode(Track.DisplayMode.SQUISHED);
        setColor(Color.DARK_GRAY);
        this.db = db;
        this.sequence = sequence;
        this.features = features;
        init();
    }

    private void init() {

        setUseScore(true);
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.source = new FeatureCollectionSource(this.features, genome);
        this.renderer = new IGVFeatureRenderer();
        IGVEventBus.getInstance().subscribe(DataLoadedEvent.class, this);
    }

    private void openTableView() {

        BlatQueryWindow win = new BlatQueryWindow(IGV.getMainFrame(), sequence, features);
        win.setVisible(true);
    }

    public List<PSLRecord> getFeatures() {
        return features;
    }

    @Override
    protected boolean isShowFeatures(ReferenceFrame frame) {
        return true;
    }

    public IGVPopupMenu getPopupMenu(final TrackClickEvent te) {
        IGVPopupMenu menu = TrackMenuUtils.getPopupMenu(Arrays.asList(this), "Menu", te);

        JMenuItem item = new JMenuItem("Open table view");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                openTableView();
            }
        });

        menu.addSeparator();
        menu.add(item);

        return menu;
    }


    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);
        element.setAttribute("db", this.db);
        element.setAttribute("sequence", this.sequence);
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("sequence")) {
            String sequence = element.getAttribute("sequence");
            if (sequence == null || sequence.length() == 0) {
                sequence = element.getAttribute("id");  // Bug in versions < 2.11.5
            }
            String db = element.getAttribute("db");
            if (db == null || db.length() == 0) {
                db = GenomeManager.getInstance().getCurrentGenome().getBlatDB();
            }
            try {
                this.features = BlatClient.blat(db, sequence);
            } catch (Exception e) {
                MessageUtils.showMessage("Error restoring blat track: " + e.getMessage());
                log.error("Error restoring blat track", e);
            }
        } else {
            // Legacy tracks
            String tmp = element.getAttribute("tokensList");
            List<String> tokensList = Arrays.asList(tmp.split("\n"));
            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            PSLCodec codec = new PSLCodec(genome, true);
            this.features = new ArrayList<>(tokensList.size());
            for (String tokens : tokensList) {
                PSLRecord f = codec.decode(tokens);
                if (f != null) {
                    this.features.add(f);
                }
            }
        }

        init();

    }

    private static String guid() {
        return "" + rand.nextInt();
    }

}
