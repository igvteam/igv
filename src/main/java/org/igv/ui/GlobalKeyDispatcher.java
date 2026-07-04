package org.igv.ui;

import org.igv.logging.*;
import org.igv.charts.ScatterPlotUtils;
import org.igv.feature.Range;
import org.igv.feature.RegionOfInterest;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.alignment.AlignmentTrackUtils;
import org.igv.alignment.SortOption;
import org.igv.ui.panel.FrameManager;
import org.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import static org.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */
public class GlobalKeyDispatcher implements KeyEventDispatcher {

    private static Logger log = LogManager.getLogger(GlobalKeyDispatcher.class);

    private final InputMap inputMap = new InputMap();
    private final ActionMap actionMap = new ActionMap();

    private static GlobalKeyDispatcher theInstance;

    public static synchronized GlobalKeyDispatcher getInstance() {
        if(theInstance == null) {
            theInstance = new GlobalKeyDispatcher();
        }
        return theInstance;
    }

    private GlobalKeyDispatcher() {
        init();
    }

    public boolean dispatchKeyEvent(KeyEvent event) {

        // If the source of this event is a text component don't process it here.
        final Object source = event.getSource();
        if(JTextComponent.class.isInstance(source) ||
                TextComponent.class.isInstance(source)) {
            return false;   // <= important, returning true will prevent further dispatching of event
        }


        KeyStroke ks = KeyStroke.getKeyStrokeForEvent(event);
        String actionKey = (String) inputMap.get(ks);

        // Disable tooltip if any modifier control key is pressed
        if (event.getKeyCode() == KeyEvent.VK_CONTROL || event.getKeyCode() == KeyEvent.VK_ALT) {
            boolean flag = !(event.isControlDown() || event.isAltDown() || event.isMetaDown());
            ToolTipManager.sharedInstance().setEnabled(flag);
        }



        if (actionKey != null) {
            Action action = actionMap.get(actionKey);
            if (action != null && action.isEnabled()) {
                // I'm not sure about the parameters
                action.actionPerformed(
                        new ActionEvent(source, event.getID(),
                                actionKey, ((KeyEvent) event).getModifiers()));
                return true; // consume event
            }
        }

        return false;
    }

    /**
     * Initialize the input and action map.   The indirection here strikes me as odd but it is apparently the standard pattern.
     */
    public void init() {

        final IGV igv = IGV.getInstance();
        final IGVPreferences prefMgr = PreferencesManager.getPreferences();


        // Set the focus to the "search" box
        final KeyStroke searchBoxKeyCtrl = KeyStroke.getKeyStroke(KeyEvent.VK_L, KeyEvent.CTRL_MASK, false);
        final KeyStroke searchBoxKeyMeta = KeyStroke.getKeyStroke(KeyEvent.VK_L, KeyEvent.META_MASK, false);
        final Action searchBoxAction = new EnableWrappedAction(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                igv.focusSearchBox();
            };
        });
        inputMap.put(searchBoxKeyCtrl, "focusSearch");
        inputMap.put(searchBoxKeyMeta, "focusSearch");
        actionMap.put("focusSearch", searchBoxAction);

        // Show extras menu
        final KeyStroke extrasKey = KeyStroke.getKeyStroke(KeyEvent.VK_T, KeyEvent.ALT_MASK, false);
        final Action extrasAction = new EnableWrappedAction(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                igv.enableExtrasMenu();
            }
        });
        inputMap.put(extrasKey, "tools");
        actionMap.put("tools", extrasAction);

        // Create region-on-interest
        final KeyStroke regionKey = KeyStroke.getKeyStroke(KeyEvent.VK_R, KeyEvent.CTRL_MASK, false);
        final Action regionAction = new EnableWrappedAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (FrameManager.isGeneListMode()) {
                    return;
                }
                Range currentRange = FrameManager.getDefaultFrame().getCurrentRange();
                RegionOfInterest regionOfInterest =
                        new RegionOfInterest(
                                currentRange.getChr(),
                                currentRange.getStart(),
                                currentRange.getEnd(),
                                null);
                igv.addRegionOfInterest(regionOfInterest);
            }
        });
        inputMap.put(regionKey, "region");
        actionMap.put("region", regionAction);

        // Create region-of-interest at center of view (1 bp wide)
        final KeyStroke regionCenterKey = KeyStroke.getKeyStroke(KeyEvent.VK_R, KeyEvent.CTRL_MASK + KeyEvent.SHIFT_MASK, false);
        final Action regionCenterAction = new EnableWrappedAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (FrameManager.isGeneListMode()) {
                    return;
                }
                int center = (int) FrameManager.getDefaultFrame().getCenter();
                RegionOfInterest regionOfInterest =
                        new RegionOfInterest(
                                FrameManager.getDefaultFrame().getChrName(),
                                center,
                                center + 1,
                                null);
                igv.addRegionOfInterest(regionOfInterest);
            }
        });
        inputMap.put(regionCenterKey, "regionCenter");
        actionMap.put("regionCenter", regionCenterAction);

        // Sort alignments
        final KeyStroke sortByLastKey = KeyStroke.getKeyStroke(KeyEvent.VK_S, KeyEvent.CTRL_DOWN_MASK, false);
        final Action sorAlignmentTracksAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                sortAlignmentTracks();
            }
        };
        inputMap.put(sortByLastKey, "sortByLast");
        actionMap.put("sortByLast", sorAlignmentTracksAction);

        // Open scatter plot
        final KeyStroke scatterplotKey = KeyStroke.getKeyStroke(KeyEvent.VK_P, KeyEvent.CTRL_MASK, false);
        final Action scatterplotAction = new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (ScatterPlotUtils.hasPlottableTracks()) {
                    ReferenceFrame defaultFrame = FrameManager.getDefaultFrame();
                    String chr = defaultFrame.getChrName();
                    int start = (int) defaultFrame.getOrigin();
                    int end = (int) defaultFrame.getEnd();
                    int zoom = defaultFrame.getZoom();
                    ScatterPlotUtils.openPlot(chr, start, end, zoom);
                }
            }
        };
        inputMap.put(scatterplotKey, "scatterPlot");
        actionMap.put("scatterPlot", scatterplotAction);

        // Back button
        final KeyStroke backKey1 = KeyStroke.getKeyStroke(KeyEvent.VK_CLOSE_BRACKET, KeyEvent.META_DOWN_MASK, false);
        final KeyStroke backKey2 = KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, KeyEvent.ALT_DOWN_MASK, false);
        final Action backAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                igv.getSession().getHistory().back();
            }
        };
        inputMap.put(backKey1, "back");

        inputMap.put(backKey2, "back");
        actionMap.put("back", backAction);

        // Forward button
        final KeyStroke forwardKey1 = KeyStroke.getKeyStroke(KeyEvent.VK_OPEN_BRACKET, KeyEvent.META_DOWN_MASK, false);
        final KeyStroke forwardKey2 = KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, KeyEvent.ALT_DOWN_MASK, false);
        final Action forwardAction = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                igv.getSession().getHistory().forward();
            }
        };
        inputMap.put(forwardKey1, "forward");
        inputMap.put(forwardKey2, "forward");

    }


    /**
     * Sort all alignment tracks according to user preference setting.
     */
    public static void sortAlignmentTracks() {
        IGVPreferences prefMgr = PreferencesManager.getPreferences();
        String sortOptionString = prefMgr.get(SAM_SORT_OPTION);
        if (sortOptionString != null) {
            try {
                SortOption option = SortOption.fromString(sortOptionString);
                String lastSortTag = prefMgr.get(SAM_SORT_BY_TAG);
                AlignmentTrackUtils.sortAlignmentTracks(option, lastSortTag, prefMgr.getAsBoolean(SAM_INVERT_SORT));
            } catch (IllegalArgumentException e1) {
                log.error("Unrecognized sort option: " + sortOptionString);
            }
        }
    }


    /**
     * TODO I'm actually pretty sure this class doesn't do what it's intended to do,
     * but I just refactored it to condense code, there were no functional changes.
     * -JS
     */
    private class EnableWrappedAction extends AbstractAction {

        private Action action;

        private EnableWrappedAction(Action action) {
            this.action = action;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            setEnabled(false); // stop any other events from interfering
            this.action.actionPerformed(e);
            setEnabled(true);

        }
    }

}
