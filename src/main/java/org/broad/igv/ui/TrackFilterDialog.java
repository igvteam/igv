package org.broad.igv.ui;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.util.TrackFilter;
import org.broad.igv.util.FilterElement;

import javax.swing.*;
import java.awt.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class TrackFilterDialog extends JDialog {

    JPanel trackFilterPane;
    JScrollPane scrollPane;

    private JCheckBox showAllTracksFilterCheckBox;
    private JCheckBox matchAllCheckBox;
    private JCheckBox matchAnyCheckBox;

    private boolean cancelled = true;
    private List<String> attributeKeys;
    private final int vgap = 5;

    public TrackFilterDialog(Frame owner, String title, TrackFilter trackFilter) {
        super(owner, title, true);
        setLocationRelativeTo(owner);
        setResizable(true);
        init(trackFilter);
        setSize(700, 600);
    }

    private void init(TrackFilter trackFilter) {

        if (trackFilter == null) {
            trackFilter = new TrackFilter();
        }

        setLayout(new BorderLayout());

        JPanel filterHeaderPanel = new JPanel();
        filterHeaderPanel.setLayout(new GridLayout(0, 1));

        showAllTracksFilterCheckBox = new JCheckBox("Show all tracks");
        showAllTracksFilterCheckBox.setSelected(trackFilter.isShowAll());
        showAllTracksFilterCheckBox.addActionListener(e -> {
            matchAllCheckBox.setEnabled(!showAllTracksFilterCheckBox.isSelected());
            matchAnyCheckBox.setEnabled(!showAllTracksFilterCheckBox.isSelected());
        });
        filterHeaderPanel.add(showAllTracksFilterCheckBox);

        matchAllCheckBox = new JCheckBox("match all of the following");
        matchAnyCheckBox = new JCheckBox("match any of the following");
        matchAllCheckBox.setSelected(trackFilter.isMatchAll());
        matchAnyCheckBox.setSelected(!trackFilter.isMatchAll());
        ButtonGroup booleanButtonGroup = new ButtonGroup();
        booleanButtonGroup.add(matchAllCheckBox);
        booleanButtonGroup.add(matchAnyCheckBox);


        if (trackFilter.isShowAll()) {
            matchAllCheckBox.setEnabled(false);
            matchAnyCheckBox.setEnabled(false);
        } else {
            matchAllCheckBox.setEnabled(true);
            matchAnyCheckBox.setEnabled(true);
        }

        JPanel controls = new JPanel();
        FlowLayout layoutManager = new FlowLayout();
        layoutManager.setAlignment(FlowLayout.LEFT);
        controls.setLayout(layoutManager);
        controls.add(new JLabel("Show tracks that: "));
        controls.add(matchAllCheckBox);
        controls.add(matchAnyCheckBox);

        controls.add(Box.createHorizontalStrut(10));
        JButton addButton = new JButton("+");
        addButton.addActionListener(e -> {
            addComponent();
        });
        controls.add(addButton);

        filterHeaderPanel.add(controls);

        add(filterHeaderPanel, BorderLayout.NORTH);

        List<String> uniqueAttributeKeys = AttributeManager.getInstance().getAttributeNames();
        this.attributeKeys = uniqueAttributeKeys;

        trackFilterPane = new JPanel();
        trackFilterPane.setLayout(new FlowLayout(FlowLayout.LEFT));

        scrollPane = new JScrollPane(trackFilterPane);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        //scrollPane.setPreferredSize(new Dimension(700, 500));

        if(trackFilter.getFilterElements().hasNext()) {
            Iterator iterator = trackFilter.getFilterElements();
            while (iterator.hasNext()) {
                FilterElement element = (FilterElement) iterator.next();
                trackFilterPane.add(new FilterComponent(this, uniqueAttributeKeys, element));
            }
        } else {
        // If there are no filters, add a blank one
            trackFilterPane.add(new FilterComponent(this, uniqueAttributeKeys, null));
        }

        add(scrollPane, BorderLayout.CENTER);

        // Button bar with OK and Cancel buttons
        JPanel buttonBar = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        JButton okButton = new JButton("OK");
        JButton cancelButton = new JButton("Cancel");

        okButton.addActionListener(e -> {
            cancelled = false;
            setVisible(false);
            dispose();
        });

        cancelButton.addActionListener(e -> {
            cancelled = true;
            setVisible(false);
            dispose();
        });

        buttonBar.add(okButton);
        buttonBar.add(cancelButton);

        add(buttonBar, BorderLayout.SOUTH);

    }

    /**
     * Add another filter
     */
    public void addComponent() {

        trackFilterPane.add(new FilterComponent(this, attributeKeys, null));
        resizeFilterPane();
    }

    public void removeComponent(FilterComponent filterComponent) {

        if (trackFilterPane.getComponents().length < 2) {
            // Can not leave less than one filter element on the screen
            return;
        }

        trackFilterPane.remove(filterComponent);
        resizeFilterPane();

    }

    private void resizeFilterPane() {
        int h = 0;
        int w = 0;
        for (Component c : trackFilterPane.getComponents()) {
            h += c.getHeight() + 2*vgap;
            w = c.getWidth();
        }
        trackFilterPane.setPreferredSize(new Dimension(w, h));
        scrollPane.getViewport().revalidate();

    }




    public boolean isCancelled() {
        return cancelled;
    }

    public TrackFilter getFilter() {
        boolean matchAll = matchAllCheckBox.isSelected();
        boolean showAllTracks = showAllTracksFilterCheckBox.isSelected();
        List<FilterElement> filterElements = new ArrayList<>();
        for (Component component : trackFilterPane.getComponents()) {
            if (component instanceof FilterComponent) {
                if (((FilterComponent) component).isComplete()) {
                    filterElements.add(((FilterComponent) component).getFilterElement(matchAll));
                }
            }
        }
        return new TrackFilter(showAllTracks, matchAll, filterElements);
    }

}

 class FilterComponent extends JPanel {


     private TrackFilterDialog filterDialog;

     private JComboBox comparisonOperatorComboBox;
     private JComboBox itemComboBox;
     private JButton removeButton;
     private JTextField valueTextField;


     public FilterComponent(TrackFilterDialog filterDialog, List<String> items, FilterElement element) {

         initComponents();

         // Load the Item ComboBox
         this.filterDialog = filterDialog;
         itemComboBox.setModel(new DefaultComboBoxModel(items.toArray()));

         // Load available operators into combobox
         List<String> textForOperators = new ArrayList<String>();
         FilterElement.Operator[] operators = FilterElement.Operator.values();
         for (int i = 0; i < operators.length; i++) {

             // List the operators to skip
             if (operators[i].equals(FilterElement.Operator.GREATER_THAN_OR_EQUAL) ||
                     operators[i].equals(FilterElement.Operator.LESS_THAN_OR_EQUAL)) {
                 continue;
             }

             textForOperators.add(operators[i].getValue());
         }
         Collections.sort(textForOperators);
         comparisonOperatorComboBox.setModel(
                 new javax.swing.DefaultComboBoxModel(textForOperators.toArray()));


         // Initialize the UI with the FilterElement if present
         if (element != null) {
             itemComboBox.setSelectedItem(element.getAttributeKey());
             comparisonOperatorComboBox.setSelectedItem(element.getComparisonOperator().getValue());
             valueTextField.setText(element.getValue());
         }

     }


     /**
      * Helper method to convert the string representation of an operator
      * to the appropriate object representation.
      */
     private FilterElement.Operator getOperatorForText(String operatorText) {

         FilterElement.Operator selected = null;

         FilterElement.Operator[] operators = FilterElement.Operator.values();
         for (FilterElement.Operator operator : operators) {

             if (operatorText.equals(operator.getValue())) {
                 selected = operator;
                 break;
             }
         }

         return selected;
     }

     /**
      * Save the UI content into a non-UI version of the FilterElement
      */
     public FilterElement getFilterElement(boolean matchAll) {
         String attributeKey = itemComboBox.getSelectedItem().toString();
         FilterElement.Operator operator = getOperatorForText((String) comparisonOperatorComboBox.getSelectedItem());
         String expectedValue = valueTextField.getText();
         return new FilterElement(attributeKey, operator, expectedValue);
     }

     public boolean isComplete() {
         return valueTextField.getText() != null && !valueTextField.getText().isEmpty();
     }

     protected void remove() {
         filterDialog.removeComponent(this);
     }


     private void initComponents() {


         setRequestFocusEnabled(false);
         setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 2, 0));

         itemComboBox = new JComboBox();
         itemComboBox.setMinimumSize(new java.awt.Dimension(50, 27));
         itemComboBox.setPreferredSize(new java.awt.Dimension(150, 27));

         comparisonOperatorComboBox = new JComboBox();
         comparisonOperatorComboBox.setActionCommand("comparisonOperatorComboBoxChanged");
         comparisonOperatorComboBox.setMinimumSize(new java.awt.Dimension(50, 27));
         comparisonOperatorComboBox.setPreferredSize(new java.awt.Dimension(150, 27));

         valueTextField = new JTextField();
         valueTextField.setMaximumSize(new Dimension(32767, 20));
         valueTextField.setMinimumSize(new Dimension(50, 27));
         valueTextField.setPreferredSize(new Dimension(150, 27));

         add(itemComboBox);
         add(comparisonOperatorComboBox);
         add(valueTextField);

         removeButton = new JButton("-");
         removeButton.setFont(new Font("Arial", 0, 14));
         removeButton.setContentAreaFilled(false);
         removeButton.setHorizontalTextPosition(SwingConstants.CENTER);
         removeButton.setPreferredSize(new Dimension(45, 27));
         removeButton.addActionListener(evt -> remove());
         add(removeButton);

     }
 }
