package org.broad.igv.ui.panel;

import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;

public class CollapsiblePanel extends JPanel {

    public static final Color HEADER_BG = new Color(210, 210, 210);


    private final JButton collapseButton;
    private ImageIcon openIcon;
    private ImageIcon closeIcon;

    public CollapsiblePanel(String label, JComponent content) {
        this(label, content, false);
    }

    public CollapsiblePanel(String label, JComponent content, boolean isOpen) {

        setLayout(new BorderLayout());

        setBorder(BorderFactory.createLineBorder(Color.BLACK));

        this.openIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.MINUS);
        this.closeIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.PLUS);

        content.setVisible(isOpen);
        this.add(content, BorderLayout.CENTER);

        JPanel header = new JPanel();
        header.setLayout(new BorderLayout());
        header.setBackground(HEADER_BG);

        this.collapseButton = new JButton();
        collapseButton.setIcon(isOpen ? openIcon : closeIcon);
        collapseButton.setBorder(new EmptyBorder(10, 5, 10, 0));
        collapseButton.setHorizontalAlignment(SwingConstants.LEFT);

        collapseButton.addActionListener(e -> {
            collapseButton.setIcon(content.isVisible() ? closeIcon : openIcon);
            content.setVisible(!content.isVisible());
        });
        header.add(collapseButton, BorderLayout.WEST);

        final JLabel jLabel = new JLabel(label);
        jLabel.setFont(FontManager.getFont(14));
        jLabel.setHorizontalAlignment(SwingConstants.CENTER);
        header.add(jLabel, BorderLayout.CENTER);

        this.add(header, BorderLayout.NORTH);

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
