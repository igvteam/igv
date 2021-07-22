package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.event.DataLoadedEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.blat.LegacyBlatClient;
import org.broad.igv.util.blat.BlatQueryWindow;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class BlatTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(BlatTrack.class);

    String sequence;
    List<PSLRecord> features;

    public BlatTrack() {
        setDisplayMode(Track.DisplayMode.SQUISHED);
        setColor(Color.DARK_GRAY);
    }

    public BlatTrack(String sequence, List<PSLRecord> features, String trackLabel) {
        super(sequence, trackLabel);
        setDisplayMode(Track.DisplayMode.SQUISHED);
        setColor(Color.DARK_GRAY);
        this.features = features;
        init();
    }

    private void init() {

        setUseScore(true);


        // Convert tokens to features
        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        if (this.features.isEmpty()) {
            MessageUtils.showMessage("No features found");
        }

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

        element.setAttribute("sequence", this.sequence);
    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("sequence")) {
            // Legacy tracks
            String sequence = element.getAttribute("sequence");
            String species = element.getAttribute("species");
            String db = element.getAttribute("db");
            try {
                this.features = BlatClient.blat(db, sequence);
            } catch (IOException e) {
                MessageUtils.showMessage("Error restoring blat track: " + e.getMessage());
                log.error("Error restoring blat track", e);
            }
        } else {
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


    public static void createBlatTrack(final String userSeq, final String trackLabel) {

        String serverType = PreferencesManager.getPreferences().get(Constants.BLAT_SERVER_TYPE);
        if (serverType.equalsIgnoreCase("web_blat")) {

            LegacyBlatClient.doBlatQuery(userSeq, trackLabel);

        } else {

            LongRunningTask.submit(new NamedRunnable() {

                public String getName() {
                    return "Blat sequence";
                }

                public void run() {
                    try {
                        Genome genome = IGV.hasInstance() ? GenomeManager.getInstance().getCurrentGenome() : null;
                        String db = genome == null ? "hg19" : genome.getUCSCId();
                        List<PSLRecord> features = BlatClient.blat(db, userSeq);
                        BlatTrack newTrack = new BlatTrack(userSeq, features, trackLabel); //species, userSeq, db, genome, trackLabel);
                        IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTrack(newTrack);
                        IGV.getInstance().repaint();
                        BlatQueryWindow win = new BlatQueryWindow(IGV.getMainFrame(), userSeq, newTrack.getFeatures());
                        win.setVisible(true);

                    } catch (Exception e1) {
                        MessageUtils.showErrorMessage("Error running blat", e1);
                    }
                }
            });
        }
    }

    public static void createBlatTrackFromRegion(final String chr, final int start, final int end, Strand strand) {

        if ((end - start) > 8000) {
            MessageUtils.showMessage("BLAT searches are limited to 8kb.  Please try a shorter sequence.");
            return;
        }

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        final byte[] seqBytes = genome.getSequence(chr, start, end);
        String userSeq = new String(seqBytes);

        if (strand == Strand.NEGATIVE) {
            userSeq = SequenceTrack.getReverseComplement(userSeq);
        }

        createBlatTrack(userSeq, "Blat");
    }


}
