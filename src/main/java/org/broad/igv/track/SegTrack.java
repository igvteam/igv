package org.broad.igv.track;

import org.broad.igv.jbrowse.CircularViewUtilities;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.HeatmapRenderer;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.VariantMenu;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

public class SegTrack extends CompositeTrack {

    CNFreqTrack freqTrack;

    public SegTrack(ResourceLocator locator, List<Track> tracks) {
        super(locator, tracks);
    }


    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {
        int trackY = 0;
        for (Track t : tracks) {
            Rectangle r = new Rectangle(trackRectangle);
            r.y = trackY;
            r.height = t.getContentHeight();
            t.renderName(g2D, r, visibleRectangle);
            trackY += r.height;
        }
    }

    @Override
    public Renderer getRenderer() {
        // All tracks have the same renderer
        return tracks.isEmpty() ? new HeatmapRenderer() : tracks.get(0).getRenderer();
    }

    @Override
    public void setRenderer(Renderer renderer) {
        for (Track t : tracks) {
            t.setRenderer(renderer);
        }
    }

    @Override
    public Collection<WindowFunction> getAvailableWindowFunctions() {
        return tracks.isEmpty() ? Collections.EMPTY_LIST : tracks.get(0).getAvailableWindowFunctions();
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        List<Track> tmp = Arrays.asList(this);

        IGVPopupMenu popupMenu = new IGVPopupMenu();

        popupMenu.add(TrackMenuUtils.getTrackRenameItem(tmp));
        popupMenu.add(TrackMenuUtils.getChangeTrackHeightItem(tmp));

        JMenuItem rowHeightItem = new JMenuItem("Change Sample Row Height...");
        rowHeightItem.addActionListener(e -> changeTrackHeight());
        popupMenu.add(rowHeightItem);

        popupMenu.addSeparator();
        TrackMenuUtils.addRendererItems(popupMenu, tmp);

        popupMenu.addSeparator();

        popupMenu.add(TrackMenuUtils.getHeatmapScaleItem(tracks));


//        popupMenu.addSeparator();
//        this.addSnpTresholdItem(popupMenu);
//
//        popupMenu.addSeparator();
//        addLoadCoverageDataItem(popupMenu);
//
//        popupMenu.addSeparator();
//        addCopyDetailsItem(popupMenu, te);
//
//        if (alignmentTrack != null) {
//            popupMenu.addSeparator();
//            addShowItems(popupMenu);
//        }

        return popupMenu;
    }

//    JLabel popupTitle = new JLabel("<html><b>" + this.track.getName(), JLabel.LEFT);
//    Font newFont = getFont().deriveFont(Font.BOLD, 12);
//        popupTitle.setFont(newFont);
//    add(popupTitle);
//
//
//        if (PreferencesManager.getPreferences().getAsBoolean(Constants.CIRC_VIEW_ENABLED) && CircularViewUtilities.ping()) {
//        addSeparator();
//        JMenuItem circItem = new JMenuItem("Add SVs to Circular View");
//        circItem.addActionListener(e1 -> track.sendToCircularView(e));
//        add(circItem);
//    }
//
//
//    //Change Track Settings
//    addSeparator();
//
//    List<Track> selectedTracks = Arrays.asList(variantTrack);
//    add(TrackMenuUtils.getTrackRenameItem(selectedTracks));
//    add(TrackMenuUtils.getChangeFontSizeItem(selectedTracks));
//

//    public JMenuItem getGenotypeSortItem(final Variant variant) {
//
//        JMenuItem item = new JMenuItem("Sort By Genotype");
//        if (variant != null) {
//            item.addActionListener(evt -> {
//                VariantMenu.GenotypeComparator compare = new VariantMenu.GenotypeComparator(variant);
//                genotypeSortingDirection = !genotypeSortingDirection;
//                track.sortSamples(compare);
//                IGV.getInstance().getContentPane().repaint();
//            });
//        }
//
//        return item;
//    }


    public void changeTrackHeight() {

        if (tracks.isEmpty()) {
            return;
        }

        final String parameter = "Sample Row height";
        int initialValue = tracks.iterator().next().getHeight();
        Integer value = TrackMenuUtils.getIntegerInput(parameter, initialValue);
        if (value == null) {
            return;
        }

        value = Math.max(0, value);
        for (Track track : tracks) {
            track.setHeight(value);
        }

        IGV.getInstance().getMainPanel().forceTracksLayout();

    }

}
