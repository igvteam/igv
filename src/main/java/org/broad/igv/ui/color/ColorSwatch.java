package org.broad.igv.ui.color;

import org.broad.igv.ui.util.IGVMouseInputAdapter;

import javax.swing.*;
import javax.swing.border.BevelBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

public class ColorSwatch extends JPanel {

    private Color selectedColor;
    private ArrayList<ColorChangeListener> listeners;

    public ColorSwatch(Color selectedColor) {
        this.selectedColor = selectedColor;
        this.listeners = new ArrayList<>();
        this.setPreferredSize(new Dimension(15, 15));
        this.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        this.setBackground(selectedColor);
        this.addMouseListener(new IGVMouseInputAdapter() {
            @Override
            public void igvMouseClicked(MouseEvent mouseEvent) {
                changeColorAction(mouseEvent);
            }
        });
    }

    public void setSelectedColor(Color selectedColor) {
        this.selectedColor = selectedColor;
        this.setBackground(selectedColor);
    }

    public void addColorChangeListener(ColorChangeListener listener) {
        listeners.add(listener);
    }

    private void changeColorAction(MouseEvent mouseEvent) {

        final JColorChooser colorChooser = new JColorChooser(selectedColor);
        JDialog dialog = JColorChooser.createDialog(
                this,
                "Select color",
                true,
                colorChooser,
                new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        Color c = colorChooser.getColor();
                        if (c != null && !c.equals(selectedColor)) {
                            setSelectedColor(c);
                            for (ColorChangeListener l : listeners) {
                                l.colorChanged(c);
                            }
                        }
                    }
                },
                null);
        dialog.setVisible(true);

    }

    public static void main(String [] args) {
        JFrame f= new JFrame("Example");
        f.setLayout(new FlowLayout());
        f.setSize(400,400);

        ColorSwatch colorSwatch = new ColorSwatch(Color.BLUE);
        colorSwatch.addColorChangeListener(c -> {
            //Color c = ((ColorSwatch) e.getSource()).getSelectedColor();
            if(c != null) {
                float [] components = c.getRGBComponents(new float[4]);
                int r = (int) (components[0] * 255);
                int g = (int) (components[1] * 255);
                int b = (int) (components[2] * 255);
                System.out.println(r + "," + g + "," + b);
            }

        });
        f.add(colorSwatch);

        f.setVisible(true);
    }

    public interface ColorChangeListener {
        void colorChanged(Color c);
    }

}
