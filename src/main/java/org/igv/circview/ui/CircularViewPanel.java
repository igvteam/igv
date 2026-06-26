package org.igv.circview.ui;

import org.igv.circview.model.ChordCollection;
import org.igv.circview.util.ColorUtils;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;

/**
 * A composed Swing UI for the circular view: a toolbar, a collapsible control
 * panel (one row per chord set or track), and the {@link CircularView} itself.
 *
 * <p>Port of createControls()/addToControlPanel() in circularView.js. The control
 * panel rebuilds itself whenever the view's collections change (chords added,
 * cleared, or regrouped), via {@link CircularView#addStructureListener}.
 */
public class CircularViewPanel extends JPanel {

    private static final double EXP5 = Math.exp(5);
    private static final Color PANEL_BG = new Color(216, 230, 234);

    /**
     * Max characters shown for a collection name. Longer names overflow the row's
     * FlowLayout and wrap out of the (single-line) visible height, so they are
     * truncated to their tail — the most meaningful part of a track name — with a
     * leading ellipsis. The full name remains available in the row's tooltip.
     */
    private static final int MAX_NAME_CHARS = 50;

    private final CircularView view;
    private final CircularViewConfig config;

    private final JPanel controlPanel = new JPanel();
    private final JButton showControlsButton = new JButton();

    public CircularViewPanel(CircularView view) {
        super(new java.awt.BorderLayout());
        this.view = view;
        this.config = view.getConfig();

        controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.Y_AXIS));
        controlPanel.setVisible(false); // hidden by default, matching the JS
        controlPanel.setAlignmentX(LEFT_ALIGNMENT);

        JComponent toolbar = createToolbar();
        toolbar.setAlignmentX(LEFT_ALIGNMENT);

        JPanel north = new JPanel();
        north.setLayout(new BoxLayout(north, BoxLayout.Y_AXIS));
        north.add(toolbar);
        north.add(controlPanel);

        add(north, java.awt.BorderLayout.NORTH);
        add(view, java.awt.BorderLayout.CENTER);

        view.addStructureListener(this::rebuildControlPanel);
        rebuildControlPanel();
    }

    /** Convenience: build the view from a config and wrap it. */
    public static CircularViewPanel create(CircularViewConfig config) {
        return new CircularViewPanel(new CircularView(config));
    }

    public CircularView getView() {
        return view;
    }

    /** Show or hide the control panel (the "Show/Hide Controls" toggle). */
    public void setControlsVisible(boolean visible) {
        controlPanel.setVisible(visible);
        updateShowControlsText();
        revalidate();
        repaint();
    }

    public boolean isControlsVisible() {
        return controlPanel.isVisible();
    }

    // ---- Toolbar ------------------------------------------------------------

    private JComponent createToolbar() {
        JPanel toolbar = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 4));

        updateShowControlsText();
        showControlsButton.addActionListener(e -> setControlsVisible(!controlPanel.isVisible()));
        toolbar.add(showControlsButton);

        JButton clearAll = new JButton("Clear All");
        clearAll.addActionListener(e -> view.clearChords());
        toolbar.add(clearAll);

        return toolbar;
    }

    private void updateShowControlsText() {
        showControlsButton.setText(controlPanel.isVisible() ? "Hide Controls" : "Show Controls");
    }

    // ---- Control panel ------------------------------------------------------

    /** Rows in the control panel: the group-by row plus one per collection. Test hook. */
    int controlPanelRowCount() {
        return controlPanel.getComponentCount();
    }

    /** The control panel's row component at the given index. Test hook. */
    Component controlPanelRow(int index) {
        return controlPanel.getComponent(index);
    }

    private void rebuildControlPanel() {
        controlPanel.removeAll();
        controlPanel.add(createGroupByRow());
        for (ChordCollection c : view.getActiveCollections()) {
            controlPanel.add(createCollectionRow(c));
        }
        controlPanel.revalidate();
        controlPanel.repaint();
    }

    private JComponent createGroupByRow() {
        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 4));
        row.setBackground(PANEL_BG);
        row.setAlignmentX(Component.LEFT_ALIGNMENT);
        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, row.getPreferredSize().height));

        JCheckBox groupBy = new JCheckBox("Group by track", view.isGroupByTrack());
        groupBy.setBackground(PANEL_BG);
        groupBy.addActionListener(e -> view.setGroupByTrack(groupBy.isSelected()));
        row.add(groupBy);
        return row;
    }

    private JComponent createCollectionRow(ChordCollection collection) {
        final String name = collection.getName();

        JPanel row = new JPanel(new FlowLayout(FlowLayout.LEFT, 6, 2));
        row.setAlignmentX(Component.LEFT_ALIGNMENT);

        // Hide / Show toggle
        JButton hideShow = new JButton(collection.isVisible() ? "Hide" : "Show");
        hideShow.addActionListener(e -> {
            boolean nowVisible = !collection.isVisible();
            view.setVisible(name, nowVisible);
            hideShow.setText(nowVisible ? "Hide" : "Show");
        });
        row.add(hideShow);

        // Color swatch -> color chooser
        JButton swatch = new JButton();
        swatch.setPreferredSize(new Dimension(28, 18));
        swatch.setBackground(opaque(collection.getColor()));
        swatch.setOpaque(true);
        swatch.setBorder(BorderFactory.createLineBorder(Color.GRAY));
        swatch.setToolTipText("Set arc color");
        row.add(swatch);

        // Alpha (transparency) slider
        JSlider alpha = new JSlider(0, 1000, alphaToValue(ColorUtils.getAlpha(collection.getColor())));
        alpha.setPreferredSize(new Dimension(120, alpha.getPreferredSize().height));
        alpha.setToolTipText("Adjust transparency of arcs");
        alpha.addChangeListener(e -> {
            float a = valueToAlpha(alpha.getValue());
            view.setColor(name, ColorUtils.setAlpha(collection.getColor(), a));
        });
        row.add(alpha);

        swatch.addActionListener(e -> {
            Color chosen = JColorChooser.showDialog(this, "Arc color for " + name,
                    opaque(collection.getColor()));
            if (chosen != null) {
                float a = ColorUtils.getAlpha(collection.getColor());
                view.setColor(name, ColorUtils.setAlpha(chosen, a));
                swatch.setBackground(opaque(collection.getColor()));
            }
        });

        JLabel label = new JLabel(displayName(name));
        label.setToolTipText(name);
        row.add(label);

        row.setMaximumSize(new Dimension(Integer.MAX_VALUE, row.getPreferredSize().height));
        return row;
    }

    /** The label text for a collection name: the tail (with a leading ellipsis) when too long. */
    private static String displayName(String name) {
        if (name == null) {
            return "";
        }
        if (name.length() <= MAX_NAME_CHARS) {
            return name;
        }
        return "…" + name.substring(name.length() - MAX_NAME_CHARS);
    }

    // ---- Alpha <-> slider value mapping (ported from circularView.js) -------

    private static float valueToAlpha(int value) {
        return (float) (Math.exp(value / 200.0) / EXP5);
    }

    private static int alphaToValue(float alpha) {
        double a = Math.max(1e-6, alpha);
        int v = (int) Math.round(200.0 * Math.log(a * EXP5));
        return Math.max(0, Math.min(1000, v));
    }

    private static Color opaque(Color c) {
        return new Color(c.getRed(), c.getGreen(), c.getBlue());
    }
}
