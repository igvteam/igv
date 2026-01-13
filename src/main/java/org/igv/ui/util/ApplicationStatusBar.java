package org.igv.ui.util;

import com.jidesoft.swing.JideBoxLayout;
import org.igv.Globals;
import org.igv.event.IGVEvent;
import org.igv.event.StopEvent;
import org.igv.logging.*;
import org.igv.oauth.OAuthProvider;
import org.igv.oauth.OAuthUtils;
import org.igv.ui.FontManager;
import org.igv.event.IGVEventBus;
import org.igv.event.IGVEventObserver;
//import org.igv.event.StopEvent;

import javax.swing.*;
import java.awt.*;
import java.text.NumberFormat;
import java.util.TimerTask;


/**
 * @author eflakes
 */
public class ApplicationStatusBar extends JPanel implements IGVEventObserver { //StatusBar {

    static Logger log = LogManager.getLogger(ApplicationStatusBar.class);
    public JButton stopButton;
    private JLabel messageBox;
    private JLabel messageBox2;
    private JLabel messageBox3;
    private JLabel memoryStatus;

    java.util.Timer timer;

    public ApplicationStatusBar() {
        initialize();
    }

    private void initialize() {

        setBackground(Globals.isDarkMode() ? Color.darkGray : new Color(240, 240, 240));
        Color messageBG = Globals.isDarkMode() ? Color.darkGray :  new Color(230, 230, 230);

        Font messageFont = FontManager.getFont(11);

        setMinimumSize(new Dimension(200, 20));
        setPreferredSize(new Dimension(800, 20));

        JideBoxLayout layout = new JideBoxLayout(this, JideBoxLayout.X_AXIS);
        layout.setGap(3);
        setLayout(layout);

        messageBox = createMessageField(messageBG, messageFont);
        messageBox.setMinimumSize(new Dimension(135, 10));
        messageBox.setPreferredSize(new Dimension(135, 20));
        add(messageBox, JideBoxLayout.FIX);

        stopButton = new JButton();
        stopButton.setBorder(BorderFactory.createLineBorder(Globals.isDarkMode() ? Color.white : Color.black));
        stopButton.addActionListener(e -> IGVEventBus.getInstance().post(new StopEvent()));
        stopButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/toolbarButtonGraphics/general/Stop16.gif")));
        stopButton.setMaximumSize(new java.awt.Dimension(16, 16));
        stopButton.setMinimumSize(new java.awt.Dimension(16, 16));
        stopButton.setPreferredSize(new java.awt.Dimension(16, 16));
        stopButton.setSize(new java.awt.Dimension(16, 16));
        stopButton.setToolTipText("Cancel load");
        stopButton.setEnabled(false);
        add(stopButton, JideBoxLayout.FIX);

        messageBox2 = createMessageField(messageBG, messageFont);
        messageBox2.setMinimumSize(new Dimension(150, 10));
        messageBox2.setPreferredSize(new Dimension(150, 20));
        add(messageBox2, JideBoxLayout.FIX);


        messageBox3 = createMessageField(messageBG, messageFont);
        messageBox3.setMinimumSize(new Dimension(165, 10));
        messageBox3.setPreferredSize(new Dimension(165, 20));
        add(messageBox3, JideBoxLayout.VARY);

        memoryStatus = createMessageField(messageBG, messageFont);
        memoryStatus.setPreferredSize(new Dimension(100, 20));
        memoryStatus.setMinimumSize(new Dimension(100, 10));
        memoryStatus.setBackground(messageBG);
        add(memoryStatus, JideBoxLayout.FIX);

        MemoryUpdateTask updateTask = new MemoryUpdateTask(memoryStatus);
        timer = new java.util.Timer();
        timer.schedule(updateTask, 0, 1000);

        // Once we get positive result from the auth flow, perhaps fetch user full name and/or profile avatar?
        IGVEventBus.getInstance().subscribe(OAuthProvider.AuthStateEvent.class, this);
    }

    public void setMessage(final String message) {
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            messageBox.setText(message);
            messageBox.validate();
            messageBox.paintImmediately(messageBox.getBounds());
        });
    }

    public void setMessage2(final String message) {
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            messageBox2.setText(message);
            messageBox.validate();
            messageBox2.paintImmediately(messageBox2.getBounds());
        });
    }

    public void setMessage3(final String message) {
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            messageBox3.setText(message);
            messageBox.validate();
            messageBox3.paintImmediately(messageBox2.getBounds());
        });
    }

    private JLabel createMessageField(Color bg, Font font) {
        JLabel messageField = new JLabel("");
        messageField.setBackground(bg);
        messageField.setFont(font);
        messageField.setBorder(BorderFactory.createLineBorder(Globals.isDarkMode() ? Color.white : Color.black));
        if(Globals.isDarkMode()) {
            messageField.setForeground(Color.white);
        }
        return messageField;

    }

    public void enableStopButton(boolean enable) {
        stopButton.setEnabled(enable);
    }

    @Override
    public void receiveEvent(IGVEvent event) {
        if (event instanceof OAuthProvider.AuthStateEvent) {

            String msg = "";
            for (OAuthProvider provider : OAuthUtils.getInstance().getAllProviders()) {
                if (provider.isLoggedIn()) {
                    if(msg.length() == 0) {
                        msg += "Logged in as: ";
                    } else {
                        msg += ";  ";
                    }
                    String user = provider.getCurrentUserEmail();
                    if(user == null || user.length() == 0) {
                        user = provider.getCurrentUserName();
                    }
                    msg += user + " via " + provider.getAuthProvider();
                }
                setMessage3(msg);
            }
        }
    }

    class MemoryUpdateTask extends TimerTask {

        JLabel textField;
        NumberFormat format;

        public MemoryUpdateTask(JLabel textField) {
            this.textField = textField;
            format = NumberFormat.getIntegerInstance();
        }

        @Override
        public void run() {
            Runtime runtime = Runtime.getRuntime();
            int freeMemory = (int) (runtime.freeMemory() / 1000000);
            int totalMemory = (int) (runtime.totalMemory() / 1000000);
            int usedMemory = (totalMemory - freeMemory);
            String um = format.format(usedMemory);
            String tm = format.format(totalMemory);
            UIUtilities.invokeOnEventThread(() -> textField.setText(um + "M of " + tm + "M"));
        }

    }
}
