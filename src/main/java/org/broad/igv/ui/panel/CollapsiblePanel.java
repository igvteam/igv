package org.broad.igv.ui.panel;

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

        //setBorder(BorderFactory.createLineBorder(Color.BLACK));

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

    @Override
    public Dimension getPreferredSize() {
        Dimension d = super.getPreferredSize();
        if(!content.isVisible()) {
            return new Dimension(d.width, header.getHeight());
        } else {
            return d;
        }
    }

    @Override
    public Dimension getMaximumSize() {
        Dimension d = super.getMaximumSize();
        if(!content.isVisible()) {
            return new Dimension(d.width, header.getHeight());
        } else {
            return d;
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
