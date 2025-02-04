package org.broad.igv.ucsc.hub;

import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

public class CollapsiblePanel extends JPanel {

    public static final Color HEADER_BG = new Color(180,204,226);


    private final JButton collapseButton;
    private final JComponent content;
    private final JPanel header;
    private ImageIcon openIcon;
    private ImageIcon closeIcon;

    public CollapsiblePanel(String label, JComponent content) {
        this(label, content, false);
    }

    public CollapsiblePanel(String label, JComponent content, boolean isOpen) {

        setLayout(new BorderLayout());
        this.content = content;

        this.openIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.MINUS);
        this.closeIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.PLUS);

        content.setVisible(isOpen);
        this.add(content, BorderLayout.CENTER);

        header = new JPanel();
        header.setLayout(new BorderLayout());
        header.setBackground(HEADER_BG);

        this.collapseButton = new JButton();
        collapseButton.setIcon(isOpen ? openIcon : closeIcon);
        collapseButton.setBorder(new EmptyBorder(10, 5, 10, 0));
        collapseButton.setHorizontalAlignment(SwingConstants.LEFT);

        collapseButton.addActionListener(e -> {
            collapseButton.setIcon(content.isVisible() ? closeIcon : openIcon);
            content.setVisible(!content.isVisible());
            this.getParent().revalidate();
        });
        header.add(collapseButton, BorderLayout.WEST);

        final JLabel jLabel = new JLabel(label);
        jLabel.setFont(FontManager.getFont(14));
        jLabel.setHorizontalAlignment(SwingConstants.CENTER);
        header.add(jLabel, BorderLayout.CENTER);

        this.add(header, BorderLayout.NORTH);

    }

    public void collapse() {
        collapseButton.setIcon(closeIcon);
        content.setVisible(false);
    }

    public void expand() {
        collapseButton.setIcon(openIcon);
        content.setVisible(true);
    }


    /**
     * Constrain the maximum height to prevent BoxLayout from needlessly resizing the panel to fill space.  This is
     * rather hardcoded for the TrackHubSelectionDialog.
     *
     * @return
     */

    @Override
    public Dimension getMaximumSize() {
        Dimension d4 = header.getMinimumSize();
        if(!content.isVisible()) {
            return new Dimension(Integer.MAX_VALUE, d4.height);
        } else {
            Dimension d5 = content.getMinimumSize();
            return new Dimension(Integer.MAX_VALUE, d4.height + (int) (1.2 * d5.height));
        }
    }

    public static void main(String[] args) {

        JComponent content = new JTextArea("alskdfjalskdjflsdkjfsaldkfjsladkfjsalkfj");
        CollapsiblePanel cp = new CollapsiblePanel("Expand/Collapse", content);

        JFrame f = new JFrame("test");
        f.setSize(500, 500);
        f.setContentPane(cp);
        f.setVisible(true);
    }
}


