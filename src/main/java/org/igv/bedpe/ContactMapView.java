package org.igv.bedpe;

import org.igv.Globals;
import org.igv.event.IGVEvent;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
import org.igv.event.ViewChange;
import org.igv.hic.ContactRecord;
import org.igv.hic.HicFile;
import org.igv.hic.Region;
import org.igv.renderer.ContinuousColorScale;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.panel.RulerPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class ContactMapView extends JPanel implements IGVEventObserver {

    private InteractionTrack track;
    private HicFile hicFile;
    private String normalization;
    private Color color;
    private Map<String, ContinuousColorScale> colorScaleCache = new HashMap<>();  // Key: "normalization_binSize"
    private Map<String, String> colorScaleDataKeys = new HashMap<>();  // Tracks which data cache key was used for each color scale
    private ReferenceFrame frame;  // Current genomic region frame (dynamic)
    private int binSize; // Bin size for the current genomic region (dynamic)
    private int startBin; // Start bin for the current genomic region (dynamic)
    private int nBins;  // Number of map bins for the current genomic region (dynamic)
    private RulerPanel rulerPanel;  // Reference to ruler panel for synchronized updates
    private JTextField maxField;
    private JSlider maxSlider;
    private double sliderMinValue = 0;
    private double sliderMaxValue = 100;

    // Cached data for rendering
    private volatile List<ContactRecord> cachedRecords;
    private volatile boolean isLoading = false;
    private volatile String loadingError = null;

    // Background executor for data fetching
    private final ExecutorService dataFetchExecutor = Executors.newSingleThreadExecutor();
    private Future<?> currentFetchTask;

    // Cache key to detect when data needs to be refetched
    private volatile String dataCacheKey = null;

    public ContactMapView(InteractionTrack track, HicFile hicFile, String normalization, ReferenceFrame frame, Color color) {

        this.track = track;
        this.hicFile = hicFile;
        this.normalization = normalization;
        this.frame = frame;
        this.color = color;
        binSize = hicFile.getBinSize(frame.getChrName(), frame.getScale());
        startBin = (int) (frame.getOrigin() / binSize);

        // Adjust frame to align with bin boundaries
        int endBin = (int) (frame.getEnd() / binSize);
        int nBins = endBin - startBin;
        this.nBins = nBins;  // Store the initial number of bins
        frame.setOrigin(startBin * binSize);
        frame.setScale(binSize);
        frame.setWidthInPixels(nBins);

        // Set initial preferred size but don't constrain min/max to allow resizing
        setPreferredSize(new Dimension(nBins, nBins));

        // Add component listener to handle resize events
        this.addComponentListener(new ComponentAdapter() {
            @Override
            public void componentResized(ComponentEvent e) {
                updateFrameForSize();
                if (rulerPanel != null) {
                    rulerPanel.repaint();
                }
                repaint();
            }
        });

        this.addMouseMotionListener(new MouseMotionAdapter() {
            @Override
            public void mouseMoved(MouseEvent e) {
                super.mouseMoved(e);
                // Calculate coordinates based on current scale
                double scaleFactor = (double) getWidth() / ContactMapView.this.nBins;
                int binX = (int) (e.getX() / scaleFactor) + startBin;
                int binY = (int) (e.getY() / scaleFactor) + startBin;
                int coordX = binX * binSize + binSize / 2;
                int coordY = binY * binSize + binSize / 2;
                track.setMarkerBounds(new int[]{coordX, coordY});

                //String tooltip = String.format("Bin1: %d (Coord: %d), Bin2: %d (Coord: %d)", binX, coordX, binY, coordY);
                //setToolTipText(tooltip);
            }
        });

        this.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseExited(MouseEvent e) {
                super.mouseExited(e);
                track.setMarkerBounds(null);
            }
        });

        track.setContactMapView(this);

        // Initial data fetch
        fetchDataAsync();

        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
    }

    /**
     * Generate a cache key based on current view parameters.
     * If this key changes, we need to refetch data.
     */
    private String generateCacheKey() {
        return String.format("%s:%d:%d:%d:%s",
                frame.getChrName(),
                (int) frame.getOrigin(),
                (int) frame.getEnd(),
                binSize,
                normalization);
    }

    /**
     * Generate a compound key for the color scale cache.
     * Combines normalization and binSize so each combination has its own cached color scale.
     */
    private String getColorScaleCacheKey() {
        return normalization + "_" + binSize;
    }

    /**
     * Fetch contact records asynchronously in a background thread.
     * When complete, triggers a repaint on the EDT.
     */
    private void fetchDataAsync() {
        String newCacheKey = generateCacheKey();

        // Check if we already have valid cached data
        if (newCacheKey.equals(dataCacheKey) && cachedRecords != null && !isLoading) {
            return; // Data is already cached
        }

        // Cancel any pending fetch task
        if (currentFetchTask != null && !currentFetchTask.isDone()) {
            currentFetchTask.cancel(true);
        }

        isLoading = true;
        loadingError = null;
        // DON'T update dataCacheKey yet - only update it when data actually arrives
        // This prevents the race condition where dataMatchesCurrentView returns true
        // but we're still using stale cachedRecords

        // Trigger repaint to show loading indicator
        repaint();

        // Capture current state for the background thread
        final String chrName = frame.getChrName();
        final int origin = (int) frame.getOrigin();
        final int end = (int) frame.getEnd();
        final int currentBinSize = binSize;
        final String fetchCacheKey = newCacheKey;

        currentFetchTask = dataFetchExecutor.submit(() -> {
            try {
                Region region = new Region(chrName, origin, end);
                List<ContactRecord> records = hicFile.getContactRecords(
                        region,
                        region,
                        "BP",
                        currentBinSize,
                        normalization,
                        true
                );

                // Update cache on EDT to avoid race conditions
                SwingUtilities.invokeLater(() -> {
                    // Always update - we have fresh data for this fetch
                    cachedRecords = records;
                    dataCacheKey = fetchCacheKey;  // NOW update dataCacheKey with the fresh data
                    isLoading = false;
                    loadingError = null;
                    repaint();
                });

            } catch (IOException e) {
                SwingUtilities.invokeLater(() -> {
                    dataCacheKey = fetchCacheKey;  // Update cache key even on error
                    isLoading = false;
                    loadingError = e.getMessage();
                    repaint();
                });
            } catch (Exception e) {
                SwingUtilities.invokeLater(() -> {
                    dataCacheKey = fetchCacheKey;  // Update cache key even on error
                    isLoading = false;
                    loadingError = "Unexpected error: " + e.getMessage();
                    repaint();
                });
            }
        });
    }

    /**
     * Shutdown the executor service. Should be called when the view is disposed.
     */
    public void dispose() {
        if (currentFetchTask != null) {
            currentFetchTask.cancel(true);
        }
        dataFetchExecutor.shutdownNow();
    }

    /**
     * Change the genomic position by supplying a new ReferenceFrame.
     * The frame will be adjusted to align with bin boundaries of the hicFile.
     * If the bin size changes, the frame scale is adjusted to fit the specified
     * genomic range within the existing component size (no resize).
     *
     * @param newFrame the new ReferenceFrame specifying the genomic region to display
     */
    public void setReferenceFrame(ReferenceFrame newFrame) {

        String currentChr = frame.getChrName();

        // Recalculate bin size for the new frame's scale
        int newBinSize = hicFile.getBinSize(newFrame.getChrName(), newFrame.getScale());

        if (newBinSize == binSize && newFrame.getChrName().equals(currentChr)) {
            // No change in bin size or chromosome, just update frame's origin
            // Update startBin to match the new origin
            startBin = (int) (newFrame.getOrigin() / binSize);
            frame.setOrigin(startBin * binSize);

            // Update frame end based on initialNBins
            int endBin = startBin + nBins;
            int genomicEnd = endBin * binSize;

            // Don't update frame.setEnd() as it's calculated, just ensure scale is correct
            updateFrameForSize();

            fetchDataAsync();
            return;
        }

        // Update to the new frame
        this.frame = new ReferenceFrame(newFrame);
        binSize = newBinSize;

        // Calculate how many bins are needed to span the requested genomic range
        int genomicSpan = (int) (newFrame.getEnd() - newFrame.getOrigin());
        int requestedNBins = (int) Math.ceil((double) genomicSpan / binSize);

        // Use the requested number of bins (this represents the genomic span)
        nBins = requestedNBins;

        // Align to bin boundaries - startBin is the first bin we'll display
        startBin = (int) (newFrame.getOrigin() / binSize);

        // Set frame origin to align with bin boundaries
        frame.setOrigin(startBin * binSize);

        // The end is determined by the requested genomic span
        int endBin = startBin + nBins;

        // Set the frame scale and width based on the requested region
        frame.setScale(binSize);
        frame.setWidthInPixels(nBins);

        // Update the frame scale to fit the requested genomic span in the current component width
        // This will adjust the scale so the genomic region fits in the available pixels
        updateFrameForSize();

        // Trigger async data fetch
        fetchDataAsync();

        // Trigger repaint to show the new region
        if (rulerPanel != null) {
            rulerPanel.repaint();
        }
        repaint();
    }

    /**
     * Update the ReferenceFrame scale and widthInPixels based on the current component size.
     * This allows the view to be responsive to window resizing while maintaining the same genomic region.
     */
    private void updateFrameForSize() {
        int currentWidth = getWidth();
        if (currentWidth <= 0) return;

        // Calculate the scale factor: how many pixels per bin
        double pixelsPerBin = (double) currentWidth / nBins;

        // Update the frame's scale to reflect the new pixels-per-base-pair ratio
        // Scale is base pairs per pixel, so we need binSize / pixelsPerBin
        double newScale = binSize / pixelsPerBin;
        frame.setScale(newScale);

        // Update width in pixels to match the current component width
        frame.setWidthInPixels(currentWidth);
    }

    @Override
    public Dimension getPreferredSize() {
        // Maintain square aspect ratio
        Dimension pref = super.getPreferredSize();
        int size = Math.min(pref.width, pref.height);
        return new Dimension(size, size);
    }

    @Override
    public Dimension getMinimumSize() {
        // Maintain square aspect ratio with a reasonable minimum
        int minSize = Math.max(100, nBins / 10);
        return new Dimension(minSize, minSize);
    }

    @Override
    public Dimension getMaximumSize() {
        // Maintain square aspect ratio
        Dimension max = super.getMaximumSize();
        int size = Math.min(max.width, max.height);
        return new Dimension(size, size);
    }

    @Override
    public void setBounds(int x, int y, int width, int height) {
        // Force square dimensions by using the smaller of width/height
        int size = Math.min(width, height);
        super.setBounds(x, y, size, size);
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2D = (Graphics2D) g.create();
        try {
            if (isLoading) {
                // Show loading indicator
                drawLoadingIndicator(g2D);
            } else if (loadingError != null) {
                // Show error message
                g2D.setColor(Color.RED);
                g2D.drawString("Error: " + loadingError, 10, 20);
            } else if (cachedRecords != null) {
                // Render from cached data
                Image img = renderMapFromCache();
                if (img != null) {
                    g2D.drawImage(img, 0, 0, null);
                }
            }
        } finally {
            g2D.dispose();
        }
    }

    /**
     * Draw a loading indicator while data is being fetched.
     */
    private void drawLoadingIndicator(Graphics2D g2D) {
        int width = getWidth();
        int height = getHeight();

        // Draw semi-transparent overlay if we have cached data
        if (cachedRecords != null) {
            Image img = renderMapFromCache();
            if (img != null) {
                g2D.drawImage(img, 0, 0, null);
            }
            g2D.setColor(new Color(255, 255, 255, 128));
            g2D.fillRect(0, 0, width, height);
        }

        // Draw loading text
        g2D.setColor(Color.DARK_GRAY);
        String message = "Loading...";
        FontMetrics fm = g2D.getFontMetrics();
        int textWidth = fm.stringWidth(message);
        int textHeight = fm.getHeight();
        g2D.drawString(message, (width - textWidth) / 2, (height + textHeight) / 2);
    }

    /**
     * Render the contact map from cached data.
     * This is called from paintComponent and only uses cached data.
     */
    private Image renderMapFromCache() {
        if (cachedRecords == null) {
            return null;
        }

        int width = getWidth();
        int height = getHeight();

        if (width <= 0 || height <= 0) {
            return null;
        }

        BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

        // Calculate scale factor for rendering
        double scaleFactor = (double) width / nBins;

        String colorScaleKey = getColorScaleCacheKey();
        ContinuousColorScale colorScale = colorScaleCache.get(colorScaleKey);

        // Check if the cached records match the current view parameters
        // If they don't, we have stale data and shouldn't compute a new color scale yet
        String currentCacheKey = generateCacheKey();
        boolean dataMatchesCurrentView = currentCacheKey.equals(dataCacheKey);

        // Check if this color scale was computed from the current data
        String colorScaleDataKey = colorScaleDataKeys.get(colorScaleKey);
        boolean colorScaleNeedsUpdate = colorScale == null || !currentCacheKey.equals(colorScaleDataKey);

        // We need to compute/recompute the color scale if:
        // 1. No color scale exists for this key, OR
        // 2. The color scale exists but was computed from different data
        if (dataMatchesCurrentView && colorScaleNeedsUpdate) {
            // We have fresh data and need to compute/update color scale from it
            // First, compute the percentile and update slider range
            double upper = computePercentileOnly(cachedRecords, .5);

            if (colorScale == null) {
                // Create new color scale
                Color lowerColor = Globals.isDarkMode() ? Color.BLACK : Color.WHITE;
                colorScale = new ContinuousColorScale(0, upper, lowerColor, color);
                colorScaleCache.put(colorScaleKey, colorScale);
            } else {
                // Update existing color scale with new max value
                colorScale.setPosEnd(upper);
            }

            // Track that this color scale was computed from the current data
            colorScaleDataKeys.put(colorScaleKey, currentCacheKey);

            // Now update text field and slider position based on the new color scale
            maxField.setText(String.format("%.2f", colorScale.getMaximum()));
            updateSliderFromValue(colorScale.getMaximum());
        } else if (colorScale == null) {
            // We don't have fresh data and no color scale exists
            // Just return null and wait for fresh data
            return null;
        }


        // paint pixels directly with scaling
        for (ContactRecord record : cachedRecords) {
            int binX = record.bin1() - startBin;
            int binY = record.bin2() - startBin;

            // Scale the bin positions to pixel positions
            int x = (int) (binX * scaleFactor);
            int y = (int) (binY * scaleFactor);

            if (x >= 0 && y >= 0) {
                Color color = colorScale.getColor(record.counts());

                // Calculate the size of each scaled bin
                int pixelSize = Math.max(1, (int) Math.ceil(scaleFactor));

                // Fill the scaled pixel area
                for (int dx = 0; dx < pixelSize && (x + dx) < width; dx++) {
                    for (int dy = 0; dy < pixelSize && (y + dy) < height; dy++) {
                        if ((x + dx) < width && (y + dy) < height) {
                            img.setRGB(x + dx, y + dy, color.getRGB());
                        }
                        // Mirror for symmetry
                        if ((y + dx) < width && (x + dy) < height) {
                            img.setRGB(y + dx, x + dy, color.getRGB());
                        }
                    }
                }
            }
        }

        return img;
    }

    /**
     * Open the ContactMatrixView in a simple popup JFrame. The frame is created with DISPOSE_ON_CLOSE
     * so closing it won't exit the application.
     */
    public static void showPopup(InteractionTrack track, HicFile hicFile, String normalization, ReferenceFrame frame, Color color) {

        ReferenceFrame refFrame = new ReferenceFrame(frame);

        SwingUtilities.invokeLater(() -> {
            JFrame jf = new JFrame(track.getDisplayName());
            jf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

            // Main panel with a border layout
            JPanel mainPanel = new JPanel(new BorderLayout());

            // ContactMatrixView
            ContactMapView contactMapView = new ContactMapView(track, hicFile, normalization, refFrame, color);
            mainPanel.add(contactMapView, BorderLayout.CENTER);

            // Create a panel to hold both the control panel and ruler panel
            JPanel topPanel = new JPanel(new BorderLayout());

            // Control panel with color scale max adjustment
            JPanel controlPanel = contactMapView.createControlPanel();
            topPanel.add(controlPanel, BorderLayout.NORTH);

            // RulerPanel
            RulerPanel rulerPanel = new RulerPanel(refFrame);
            Dimension rulerSize = new Dimension(contactMapView.getPreferredSize().width, 50);
            rulerPanel.setPreferredSize(rulerSize);
            topPanel.add(rulerPanel, BorderLayout.CENTER);

            mainPanel.add(topPanel, BorderLayout.NORTH);

            // Store ruler panel reference for synchronized updates
            contactMapView.rulerPanel = rulerPanel;

            // Add component listener to keep ruler panel width in sync with contact matrix view
            contactMapView.addComponentListener(new ComponentAdapter() {
                @Override
                public void componentResized(ComponentEvent e) {
                    int width = contactMapView.getWidth();
                    rulerPanel.setPreferredSize(new Dimension(width, 50));
                    rulerPanel.revalidate();
                }
            });

            // Add window listener to unsubscribe from event bus when window closes
            jf.addWindowListener(new WindowAdapter() {
                @Override
                public void windowClosing(WindowEvent e) {
                    IGVEventBus.getInstance().unsubscribe(contactMapView);
                    contactMapView.dispose();
                    track.setContactMapView(null);
                }
            });

            jf.getContentPane().add(mainPanel);
            jf.pack();
            jf.setLocationRelativeTo(null);
            jf.setVisible(true);
        });
    }

    /**
     * Create a control panel with input controls for adjusting visualization parameters.
     * Currently includes a max value control for the color scale with both slider and text field.
     * This panel can be extended with additional controls in the future.
     */
    private JPanel createControlPanel() {
        JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

        // Label for max value
        JLabel maxLabel = new JLabel("Color Scale Max:");
        controlPanel.add(maxLabel);

        // Slider for max value - instance variable for this view
        // Initial range 0-100, will be updated dynamically
        this.maxSlider = new JSlider(0, 100, 50);
        this.maxSlider.setPreferredSize(new Dimension(200, 25));
        this.maxSlider.setToolTipText(String.format("Range: %.2f - %.2f", sliderMinValue, sliderMaxValue));
        this.maxSlider.addChangeListener(e -> {
            if (!this.maxSlider.getValueIsAdjusting()) {
                // Map slider value (0-100) to actual range (sliderMinValue to sliderMaxValue)
                double value = sliderMinValue + (sliderMaxValue - sliderMinValue) * this.maxSlider.getValue() / 100.0;
                this.maxField.setText(String.format("%.2f", value));

                // Update tooltip to show current value
                this.maxSlider.setToolTipText(String.format("Value: %.2f (Range: %.2f - %.2f)",
                        value, sliderMinValue, sliderMaxValue));

                // Update the color scale
                String colorScaleKey = getColorScaleCacheKey();
                ContinuousColorScale colorScale = this.colorScaleCache.get(colorScaleKey);
                if (colorScale != null) {
                    colorScale.setPosEnd(value);
                    this.repaint();
                }
            }
        });
        controlPanel.add(this.maxSlider);

        // Text field for max value - instance variable for this view
        this.maxField = new JTextField("10", 10);
        this.maxField.addActionListener(e -> {
            try {
                double newMax = Double.parseDouble(this.maxField.getText().trim());
                String colorScaleKey = getColorScaleCacheKey();
                ContinuousColorScale colorScale = this.colorScaleCache.get(colorScaleKey);
                if (colorScale != null) {
                    colorScale.setPosEnd(newMax);
                    this.maxField.setText(String.format("%.2f", newMax));

                    // Update slider position based on text field value
                    updateSliderFromValue(newMax);

                    this.repaint();
                }
            } catch (NumberFormatException ex) {
                JOptionPane.showMessageDialog(controlPanel,
                        "Invalid number format. Please enter a valid number.",
                        "Input Error",
                        JOptionPane.ERROR_MESSAGE);
            }
        });
        controlPanel.add(this.maxField);

        // Apply button
        JButton applyButton = new JButton("Apply");
        applyButton.addActionListener(e -> {
            try {
                String colorScaleKey = getColorScaleCacheKey();
                ContinuousColorScale colorScale = this.colorScaleCache.get(colorScaleKey);
                double newMax = Double.parseDouble(this.maxField.getText().trim());
                if (colorScale != null) {
                    colorScale.setPosEnd(newMax);
                    // Update the maxField value with the newMax
                    this.maxField.setText(String.format("%.2f", newMax));

                    // Update slider position based on text field value
                    updateSliderFromValue(newMax);

                    this.repaint();
                }
            } catch (NumberFormatException ex) {
                JOptionPane.showMessageDialog(controlPanel,
                        "Invalid number format. Please enter a valid number.",
                        "Input Error",
                        JOptionPane.ERROR_MESSAGE);
            }
        });
        controlPanel.add(applyButton);

        return controlPanel;
    }

    /**
     * Update slider position based on a value within the current range.
     * If value exceeds sliderMaxValue, expand the range.
     */
    private void updateSliderFromValue(double value) {
        // Expand slider range if value exceeds current maximum
        if (value > sliderMaxValue) {
            sliderMaxValue = value;
        }

        if (sliderMaxValue > sliderMinValue) {
            double normalized = (value - sliderMinValue) / (sliderMaxValue - sliderMinValue);
            int sliderPos = (int) Math.round(normalized * 100);
            sliderPos = Math.max(0, Math.min(100, sliderPos)); // Clamp to 0-100
            this.maxSlider.setValue(sliderPos);

            // Update tooltip to show current value and range
            this.maxSlider.setToolTipText(String.format("Value: %.2f (Range: %.2f - %.2f)",
                    value, sliderMinValue, sliderMaxValue));
        }
    }

    /**
     * Update the slider range based on min and max values from the filtered data.
     */
    private void updateSliderRange(double minValue, double maxValue) {
        this.sliderMinValue = minValue;
        this.sliderMaxValue = maxValue;

        // Update slider position to maintain current value if possible
        String colorScaleKey = getColorScaleCacheKey();
        ContinuousColorScale colorScale = this.colorScaleCache.get(colorScaleKey);
        if (colorScale != null) {
            updateSliderFromValue(colorScale.getMaximum());
        }
    }

    /**
     * Format bin size for display (e.g., 1000 -> "1 kb", 1000000 -> "1 Mb")
     */
    private static String formatBinSize(int binSize) {
        if (binSize >= 1000000) {
            return String.format("%.1f Mb", binSize / 1000000.0);
        } else if (binSize >= 1000) {
            return String.format("%.1f kb", binSize / 1000.0);
        } else {
            return binSize + " bp";
        }
    }

    @Override
    public void receiveEvent(IGVEvent event) {
        if (event instanceof ViewChange) {
            ViewChange vc = (ViewChange) event;
            if (!vc.panning) {
                setReferenceFrame(vc.referenceFrame);
            }
        }
    }

    /**
     * Compute a percentile of the "counts" value in contactRecords,
     * filtering records where Math.abs(bin2 - bin1) < 5.
     * Also updates the slider range based on the min/max of filtered counts.
     * Does NOT update the slider position - caller should do that after updating the color scale.
     *
     * @return the percentile value, or -1 if contactRecords is empty
     */
    private double computePercentileOnly(List<ContactRecord> contactRecords, double percentile) {
        if (contactRecords == null || contactRecords.isEmpty()) {
            return -1; // Return -1 if there are no records
        }

        // Filter records where Math.abs(bin2 - bin1) < 5
        List<Float> counts = contactRecords.stream()
                .filter(record -> Math.abs(record.bin2() - record.bin1()) < 5)
                .map(ContactRecord::counts)
                .sorted()
                .collect(Collectors.toList());

        if (counts.isEmpty()) {
            return -1; // Return -1 if no records pass the filter
        }

        int index = (int) Math.ceil(percentile * counts.size()) - 1;
        float midCount = counts.get(index);

        // Update slider range based on min and max of filtered counts
        float minCount = counts.get(0);
        float maxCount = 2 * midCount;
        this.sliderMinValue = minCount;
        this.sliderMaxValue = maxCount;


        return midCount;
    }

    public void setNormalization(String type) {
        if (this.normalization != null && this.normalization.equals(type)) {
            return; // No change
        }

        this.normalization = type;

        // Clear cached data since normalization affects the contact record values
        this.cachedRecords = null;
        this.dataCacheKey = null;

        // No need to clear colorScaleCache since it uses compound keys (normalization_binSize)
        // Each normalization + binSize combination maintains its own color scale

        // Trigger async data fetch with new normalization
        fetchDataAsync();
    }
}
