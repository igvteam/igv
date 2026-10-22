package org.igv.ucsc.hub;

import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.event.MenuDragMouseEvent;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.net.URI;

/**
 * Menu item for a track hub in the "Tracks" menu.  If the hub has a description url an "info" icon is drawn at the
 * right edge, styled as in the hub selection dialog: the url is its tooltip and clicking it opens the url in a
 * browser instead of selecting the item.
 */
public class HubMenuItem extends JMenuItem {

    private static Logger log = LogManager.getLogger(HubMenuItem.class);

    private static final int ICON_GAP = 8;

    private final String infoUrl;
    private final Icon infoIcon;

    public HubMenuItem(Action action, String infoUrl) {
        super(action);
        this.infoUrl = infoUrl == null || infoUrl.isBlank() ? null : infoUrl;
        this.infoIcon = this.infoUrl == null ? null : IconFactory.getInstance().getIcon(IconFactory.IconID.INFO);
    }

    private Rectangle iconBounds() {
        if (infoIcon == null) {
            return null;
        }
        Insets insets = getInsets();
        return new Rectangle(getWidth() - insets.right - infoIcon.getIconWidth(),
                (getHeight() - infoIcon.getIconHeight()) / 2, infoIcon.getIconWidth(), infoIcon.getIconHeight());
    }

    private boolean overIcon(MouseEvent e) {
        Rectangle bounds = iconBounds();
        return bounds != null && bounds.contains(e.getPoint());
    }

    @Override
    public Dimension getPreferredSize() {
        Dimension d = super.getPreferredSize();
        if (infoIcon != null) {
            d.width += ICON_GAP + infoIcon.getIconWidth();
        }
        return d;
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Rectangle bounds = iconBounds();
        if (bounds != null) {
            infoIcon.paintIcon(this, g, bounds.x, bounds.y);
        }
    }

    @Override
    public String getToolTipText(MouseEvent e) {
        return overIcon(e) ? infoUrl : super.getToolTipText(e);
    }

    @Override
    protected void processMouseMotionEvent(MouseEvent e) {
        setCursor(overIcon(e) ? Cursor.getPredefinedCursor(Cursor.HAND_CURSOR) : Cursor.getDefaultCursor());
        super.processMouseMotionEvent(e);
    }

    /**
     * Mouse clicks on the icon open the description url and close the menu.  They are not passed on, so the item's
     * action is not fired.
     */
    @Override
    protected void processMouseEvent(MouseEvent e) {
        int id = e.getID();
        if (overIcon(e) && (id == MouseEvent.MOUSE_PRESSED || id == MouseEvent.MOUSE_RELEASED || id == MouseEvent.MOUSE_CLICKED)) {
            if (id == MouseEvent.MOUSE_RELEASED) {
                MenuSelectionManager.defaultManager().clearSelectedPath();
                openInfo();
            }
            return;
        }
        super.processMouseEvent(e);
    }

    @Override
    public void processMenuDragMouseEvent(MenuDragMouseEvent e) {
        if (overIcon(e) && e.getID() == MouseEvent.MOUSE_RELEASED) {
            MenuSelectionManager.defaultManager().clearSelectedPath();
            openInfo();
            return;
        }
        super.processMenuDragMouseEvent(e);
    }

    protected void openInfo() {
        try {
            Desktop.getDesktop().browse(new URI(infoUrl));
        } catch (Exception ex) {
            log.error("Error opening hub description url: " + infoUrl, ex);
        }
    }
}
