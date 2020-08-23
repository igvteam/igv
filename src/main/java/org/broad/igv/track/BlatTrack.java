package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.event.DataLoadedEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.blat.BlatClient;
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

    String species;

    String db;

    String sequence;

    Genome genome;

    List<PSLRecord> features;

    public BlatTrack(String species, String sequence, String db, Genome genome, String trackLabel) {

        super(sequence, trackLabel);

        this.sequence = sequence;
        this.db = db;
        this.species = species;
        this.genome = genome;

        setDisplayMode(Track.DisplayMode.SQUISHED);
        setColor(Color.DARK_GRAY);
        init();
    }


    public BlatTrack() {
        this.genome = GenomeManager.getInstance().getCurrentGenome();
    }

    private void init() {

        setUseScore(true);

        try {
            List<String> tokensList = BlatClient.blat(species, db, sequence);

            // Convert tokens to features
            PSLCodec codec = new PSLCodec(genome, true);

            features = new ArrayList<PSLRecord>(tokensList.size());
            for (String tokens : tokensList) {
                PSLRecord f = codec.decode(tokens);
                if (f != null) {
                    features.add(f);
                }
            }

            if (features.isEmpty()) {
                MessageUtils.showMessage("No features found");
            }

            this.source = new FeatureCollectionSource(features, genome);

            this.renderer = new IGVFeatureRenderer();
        } catch (IOException e) {
            log.error("Error blatting sequence", e);
            MessageUtils.showErrorMessage("Error blatting sequence", e);
        }

        IGVEventBus.getInstance().subscribe(DataLoadedEvent.class, this);
    }

    private void openTableView() {

        BlatQueryWindow win = new BlatQueryWindow(IGV.getMainFrame(), sequence, features);
        win.setVisible(true);
    }

    public List<PSLRecord> getFeatures() {
        return features;
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

        element.setAttribute("db", db);
        element.setAttribute("sequence", sequence);
        element.setAttribute("species", species);

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        this.sequence = element.getAttribute("sequence");
        this.species = element.getAttribute("species");
        this.db = element.getAttribute("db");
        this.genome = GenomeManager.getInstance().getCurrentGenome();

        init();

    }


}
