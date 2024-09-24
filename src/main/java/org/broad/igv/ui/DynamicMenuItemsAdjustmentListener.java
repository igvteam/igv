package org.broad.igv.ui;

import javax.swing.*;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Function;

/**
 * MenuListener which updates the availble menu items based on the contents of a changeable collection.
 *
 * Whenever the menu is selected the relevant JMenuItems are regenerated according to the current values in the backing
 * collection.
 *
 * Items are inserted in a row beneath a given separator.
 *
 * @param <T> element type of the collection
 */
public class DynamicMenuItemsAdjustmentListener<T> implements MenuSelectedListener {
    private final JMenu menu;

    private final JSeparator insertionPoint;
    private final Collection<T> values;
    private final Function<T, JMenuItem> itemConstructor;

    //Store the currently visible components here so they can be removed when necessary
    private final List<JComponent> activeComponents;

    /**
     *
     * @param menu the menu to modify
     * @param insertionPoint a JSeparator which acts as an anchor to always insert elements below.  This element is hidden
     *                       when the collection is empty
     * @param values a collection which is used to generate menu items
     * @param itemConstructor a function to create a JMenuItem from an element in the collection
     */
    public DynamicMenuItemsAdjustmentListener(JMenu menu, JSeparator insertionPoint, Collection<T> values, Function<T, JMenuItem> itemConstructor) {
        this.menu = menu;
        this.insertionPoint = insertionPoint;
        this.values = values;
        this.itemConstructor = itemConstructor;
        this.activeComponents = new ArrayList<>();
    }

    private List<JMenuItem> getCurrentItems() {
        return values.stream().map(itemConstructor).toList();
    }

    @Override
    public void menuSelected(MenuEvent e) {
        List<JMenuItem> newComponents = getCurrentItems();

        // We definitely don't want to be doing this while other things are also changing the menu
        // this should at least protect against multiple of these listeners modifying the same menu at once.
        synchronized (menu) {
            activeComponents.forEach(menu::remove);
            if (newComponents.isEmpty()) {
                insertionPoint.setVisible(false);
            } else {
                insertionPoint.setVisible(true);

                final int componentIndex = Arrays.asList(menu.getMenuComponents()).indexOf(insertionPoint);
                for (int i = 0; i < newComponents.size(); i++) {
                    menu.insert(newComponents.get(i), componentIndex + i + 1);
                }
                activeComponents.addAll(newComponents);
            }
        }

        menu.revalidate();
        menu.repaint();
    }
}
