package org.broad.igv.track;

import org.broad.igv.event.DataLoadedEvent;
import org.broad.igv.event.IGVEventBus;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.renderer.IGVFeatureRenderer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.blat.BlatClient;
import org.broad.igv.util.blat.BlatQueryWindow;
import org.w3c.dom.Element;

import javax.swing.*;
import javax.xml.bind.annotation.XmlType;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.*;
import java.util.List;


public class BlatTrack extends FeatureTrack {


    protected Class clazz = BlatTrack.class;

    String species;

    String db;

    String sequence;

    Genome genome;

    List<PSLRecord> features;

    public BlatTrack(String species, String sequence, String db, Genome genome) throws IOException {

        super(sequence, "Blat");

        this.sequence = sequence;
        this.db = db;
        this.species = species;
        this.genome = genome;

        setUseScore(true);
        setDisplayMode(Track.DisplayMode.SQUISHED);

        init();
    }

    private void init() throws IOException {

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
        ;

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


    public Map<String, String> getPersistentState() {

        Map<String, String> state = new LinkedHashMap<String, String>();

        state.put("name", getName());
        state.put("db", db);
        state.put("sequence", sequence);
        state.put("species", species);
        state.put("displayMode", getDisplayMode().toString());

        Color c = getColor();
        if (c != null) {
            String cs = ColorUtilities.colorToString(c);
            state.put("color", cs);
        }

        return state;
    }


    public BlatTrack(Element element) throws IOException {

        super(element);

        String sequence = element.getAttribute("sequence");
        String species = element.getAttribute("species");
        String db = element.getAttribute("db");
        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        this.id = sequence;
        this.sequence = sequence;
        this.db = db;
        this.species = species;
        this.genome = genome;

        setUseScore(true);

        init();


    }

}
