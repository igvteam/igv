package org.igv.hic;

import org.igv.renderer.ContinuousColorScale;
import org.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.util.List;

public class ContactMatrixView extends JPanel {

    HicFile hicFile;
    ReferenceFrame frame;
    ContinuousColorScale colorScale;
    int binSize;

    public ContactMatrixView(HicFile hicFile, ReferenceFrame frame, ContinuousColorScale colorScale) {
        this.hicFile = hicFile;
        this.frame = frame;
        this.colorScale = colorScale;
        binSize = hicFile.getBinSize(frame.getChrName(), frame.getScale());

        int dim = (int) ((frame.getEnd() - frame.getOrigin()) / binSize);
        Dimension pref = new Dimension(dim, dim);
        setPreferredSize(pref);
        setMinimumSize(pref);
        setMaximumSize(pref);
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2D = (Graphics2D) g.create();
        try {
            renderMap(g2D);
        } catch (IOException e) {
            // Draw a simple error message on the panel if rendering fails
            g2D.setColor(Color.RED);
            String msg = "Error rendering contact matrix: " + e.getMessage();
            g2D.drawString(msg, 10, 20);
        } finally {
            g2D.dispose();
        }
    }

    private void renderMap(Graphics2D g2D) throws IOException {
        // Rendering logic for contact matrix


        final Region region1 = new Region(frame.getChrName(), (int) frame.getOrigin(), (int) frame.getEnd());

        List<ContactRecord> records = hicFile.getContactRecords(
                region1,
                region1,
                "BP",
                binSize,
                true
        );

        int startBin = region1.start() / binSize;

        // create an offscreen image sized to the view
        int width = getWidth();
        int height = getHeight();

        // Guard against zero-sized component (avoid IllegalArgumentException from BufferedImage)
        if (width <= 0 || height <= 0) {
            return;
        }

        java.awt.image.BufferedImage img = new java.awt.image.BufferedImage(width, height, java.awt.image.BufferedImage.TYPE_INT_ARGB);

        // paint pixels directly
        for (ContactRecord record : records) {
            int binX = record.bin1() - startBin;
            int binY = record.bin2() - startBin;
            int x = binX; // rectSize == 1
            int y = binY;
            if (x >= 0 && x < width && y >= 0 && y < height) {
                Color color = colorScale.getColor(record.counts());
                img.setRGB(x, y, color.getRGB());
                img.setRGB(y, x, color.getRGB());
            }
        }

//        // Max pooling
//        int bMax = 0;
//        for (ContactRecord record : records) {
//            bMax = Math.max(bMax, Math.max(record.bin1(), record.bin2()));
//        }
//        int bWidth = bMax - startBin + 1;
//        float scale = (float) width / (float) bWidth;
//
//        float[][] maxValues = new float[width][height];
//        for (ContactRecord record : records) {
//            int x = (int)((record.bin1() - startBin) * scale);
//            int y = (int)((record.bin2() - startBin) * scale);
//            if (x >= 0 && x < width && y >= 0 && y < height) {
//                maxValues[x][y] = Math.max(maxValues[x][y], record.counts());
//            }
//        }
//
//        // paint pixels based on max pooled values
//        for (int x = 0; x < width; x++) {
//            for (int y = 0; y < height; y++) {
//                if (maxValues[x][y] > 0) {
//                    Color color = colorScale.getColor(maxValues[x][y]);
//                    img.setRGB(x, y, color.getRGB());
//                    img.setRGB(y, x, color.getRGB());
//                }
//            }
//        }



        // draw the composed image once
        g2D.drawImage(img, 0, 0, null);
    }

    /**
     * Open the ContactMatrixView in a simple popup JFrame. The frame is created with DISPOSE_ON_CLOSE
     * so closing it won't exit the application. The panel will be initialized with the provided
     * hicFile, frame, and colorScale and will prefer 800x800 pixels as requested.
     */
    public static void showPopup(HicFile hicFile, ReferenceFrame frame, ContinuousColorScale colorScale) {
        SwingUtilities.invokeLater(() -> {
            JFrame jf = new JFrame("Contact Matrix");
            jf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
            ContactMatrixView panel = new ContactMatrixView(hicFile, frame, colorScale);
            jf.getContentPane().add(panel);
            jf.pack();
            jf.setLocationRelativeTo(null);
            jf.setVisible(true);
        });
    }

}
