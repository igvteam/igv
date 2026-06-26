package org.igv.circview.ui;

import org.igv.circview.model.Assembly;
import org.igv.circview.model.Chord;
import org.igv.circview.model.Chromosome;
import org.igv.circview.model.Mate;
import org.junit.Test;

import javax.swing.JCheckBox;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Verifies the control panel wiring: it rebuilds in response to the view's
 * structure changes (chords added/cleared, regrouped) and always carries the
 * group-by checkbox row.
 */
public class CircularViewPanelTest {

    private static Assembly assembly() {
        return new Assembly("g", "g", List.of(
                new Chromosome("1", 1000, Color.RED),
                new Chromosome("2", 1000, Color.BLUE)));
    }

    private static List<Chord> chords() {
        return List.of(new Chord("x", "1", 1, 2, new Mate("2", 3, 4), null));
    }

    private static CircularViewPanel panelWithTwoSets() {
        CircularView view = new CircularView(new CircularViewConfig());
        view.setAssembly(assembly());
        view.addChords(chords(), "Set A", Color.BLUE);
        view.addChords(chords(), "Set B", Color.RED);
        return new CircularViewPanel(view);
    }

    @Test
    public void firstRowIsGroupByCheckbox() {
        CircularViewPanel panel = panelWithTwoSets();
        assertTrue("first control row should hold the group-by checkbox",
                containsCheckbox((Container) panel.controlPanelRow(0)));
    }

    @Test
    public void oneRowPerCollectionPlusGroupByRow() {
        CircularViewPanel panel = panelWithTwoSets();
        // group-by row + 2 chord sets
        assertEquals(3, panel.controlPanelRowCount());
    }

    @Test
    public void rebuildsWhenChordsAdded() {
        CircularView view = new CircularView(new CircularViewConfig());
        view.setAssembly(assembly());
        CircularViewPanel panel = new CircularViewPanel(view);
        assertEquals(1, panel.controlPanelRowCount()); // group-by row only

        view.addChords(chords(), "Set A", Color.BLUE);
        assertEquals(2, panel.controlPanelRowCount());
    }

    @Test
    public void rebuildsWhenCleared() {
        CircularViewPanel panel = panelWithTwoSets();
        assertEquals(3, panel.controlPanelRowCount());
        panel.getView().clearChords();
        assertEquals(1, panel.controlPanelRowCount());
    }

    @Test
    public void regroupingByTrackRebuilds() {
        // Two chord sets sharing one track name ("Calls") collapse to a single
        // track row when grouped.
        CircularView view = new CircularView(new CircularViewConfig());
        view.setAssembly(assembly());
        view.addChords(chords(), "Calls region1", Color.BLUE);
        view.addChords(chords(), "Calls region2", Color.RED);
        CircularViewPanel panel = new CircularViewPanel(view);
        assertEquals(3, panel.controlPanelRowCount()); // group-by + 2 sets

        view.setGroupByTrack(true);
        assertEquals(2, panel.controlPanelRowCount()); // group-by + 1 track
    }

    private static boolean containsCheckbox(Container c) {
        for (Component child : c.getComponents()) {
            if (child instanceof JCheckBox) {
                return true;
            }
            if (child instanceof Container && containsCheckbox((Container) child)) {
                return true;
            }
        }
        return false;
    }
}
