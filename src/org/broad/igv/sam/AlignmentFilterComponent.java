package org.broad.igv.sam;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import com.jidesoft.swing.JideButton;

import net.sf.samtools.SAMRecord.SAMTagAndValue;

import org.broad.igv.util.FilterElement.Operator;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class AlignmentFilterComponent extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 3328021930786326253L;
	private JPanel jPanel1;
	private javax.swing.JButton moreButton;
	private javax.swing.JButton removeButton;

	final static String m_number = "Number";
	final static String m_string = "String";
	final static String m_missing = "missing value";
	static final String m_hasWildCards = "has wild cards";
	static final String m_isRegExprs = "is regular experssion";
	static final String m_asString = "interprete as string";
	static final String m_include = "include alignments";
	static final String m_exclude = "exclude alignments";
	static final String m_min = "minimum value";
	static final String m_max = "maximum value";
	static final String m_caseSensitive = "case sensitive";
	static final String m_filterExpression = "Filter expression";

	JDialog dialog;
	JPanel allFilterPanel;
	JFormattedTextField minField, maxField;
	JTextField expressionString;
	JRadioButton hasWildCardsButton, isRegExprButton, asStringButton;
	JCheckBox caseSensitiveBox;
	JRadioButton includeButton, excludeButton;
	private JPanel newP;
	private JTextField tagValue;
	JRadioButton numberButton, stringButton, missingButton;
	JComboBox comparisonOperatorComboBox;

	public boolean getNumberButton() {
		return numberButton.isSelected();
	}

	public boolean getStringButton() {
		return stringButton.isSelected();
	}

	public boolean getMissingButton() {
		return missingButton.isSelected();
	}

	public Double getMinField() {
		return Double.valueOf(minField.getText());
	}

	public Double getMaxField() {
		return Double.valueOf(maxField.getText());
	}

	public String getExpressionString() {
		return expressionString.getText();
	}

	public boolean getHasWildCardsButton() {
		return hasWildCardsButton.isSelected();
	}

	public boolean getIsRegExprButton() {
		return isRegExprButton.isSelected();
	}

	public boolean getAsStringButton() {
		return asStringButton.isSelected();
	}

	public boolean getCaseSensitiveBox() {
		return caseSensitiveBox.isSelected();
	}

	public boolean getIncludeButton() {
		return includeButton.isSelected();
	}

	public boolean getExcludeButton() {
		return excludeButton.isSelected();
	}

	public String getTagValue() {
		return tagValue.getText();
	}

	public Integer getTagValueWidth(){
		return tagValue.getWidth();
	}
	public void setTagValueWidth(int witdth, int height){
		tagValue.setSize(witdth, height);
	}

	public AlignmentFilterComponent(final AlignmentFilter alf, String index,
			ActionListener actListener) {
		super();
		jPanel1 = new javax.swing.JPanel();

		setBackground(new java.awt.Color(255, 255, 255));
		setMinimumSize(new java.awt.Dimension(530, 40));
		setPreferredSize(new java.awt.Dimension(700, 40));
		setRequestFocusEnabled(false);
		setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT, 2, 5));

		jPanel1.setBackground(new java.awt.Color(255, 255, 255));
		jPanel1.setMinimumSize(new java.awt.Dimension(470, 31));
		jPanel1.setLayout(new java.awt.FlowLayout(java.awt.FlowLayout.LEFT));

		tagValue = new JTextField(4);
		tagValue.setActionCommand("Tag" + "_" + index);
		//tagValue.addPropertyChangeListener("value",propListener);
		// JLabel tagFieldLabel = new JLabel("tag to be used:");
		// tagFieldLabel.setLabelFor(tagValue);
		tagValue.setText(alf.getTag());
		tagValue.setMinimumSize(new java.awt.Dimension(50, 27));
		tagValue.setPreferredSize(new java.awt.Dimension(150, 27));
		tagValue.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void changedUpdate(DocumentEvent e) {
				    warn();
				  }
			  @Override
				  public void removeUpdate(DocumentEvent e) {
				    warn();
				  }
				  @Override
				  public void insertUpdate(DocumentEvent e) {
				    warn();
				  }

				  public void warn() {
				     alf.setTag(tagValue.getText());
				  }
				
				});
		
		//tagValue.addActionListener(actListener);
		System.out.println("tag bounds:" + tagValue.getBounds());
		jPanel1.add(tagValue);

		List<String> textForOperators = new ArrayList<String>();
		AlignmentFilter.Operator[] operators = AlignmentFilter.Operator
				.values();
		for (int i = 0; i < operators.length; i++) {

			// List the operators to skip
			if (operators[i].equals(Operator.GREATER_THAN_OR_EQUAL)
					|| operators[i].equals(Operator.LESS_THAN_OR_EQUAL)) {
				continue;
			}

			textForOperators.add(operators[i].getValue());
		}
		Collections.sort(textForOperators);
		comparisonOperatorComboBox = new JComboBox();
		comparisonOperatorComboBox.setModel(new DefaultComboBoxModel(
				textForOperators.toArray()));

		comparisonOperatorComboBox
				.setActionCommand("compOp" + "_" + index);
		comparisonOperatorComboBox
				.setMinimumSize(new java.awt.Dimension(50, 27));
		comparisonOperatorComboBox.setPreferredSize(new java.awt.Dimension(150,
				27));
		comparisonOperatorComboBox.addActionListener(actListener);
		comparisonOperatorComboBox.setSelectedItem(alf.getOperation());
		jPanel1.add(comparisonOperatorComboBox);

		// String parameters
		expressionString = new JTextField(20);
		expressionString.setText(alf.getExpression());
		//expressionString.addPropertyChangeListener(propListener);
		
		expressionString.getDocument().addDocumentListener(new DocumentListener() {
			@Override
			public void changedUpdate(DocumentEvent e) {
				    warn();
				  }
			  @Override
				  public void removeUpdate(DocumentEvent e) {
				    warn();
				  }
				  @Override
				  public void insertUpdate(DocumentEvent e) {
				    warn();
				  }

				  public void warn() {
				     alf.setExpression(expressionString.getText());
				  }
				
				});
		jPanel1.add(expressionString);
		
		final JCheckBox inclExcl = new JCheckBox("");
		inclExcl.setActionCommand("incl" + "_" + index);
		inclExcl.setSelected(alf.isExclude());
		inclExcl.setToolTipText("if selected matching reads will be displayed\\nif not selected matching reads will not be displayed");
		inclExcl.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				alf.setExclude(inclExcl.isSelected());
				System.out.println(e);
				
			}
			
		});
		jPanel1.add(inclExcl);
		
		add(jPanel1);
		
		moreButton = new JideButton("+" + "_" + index);
		//moreButton.setActionCommand();
		moreButton.setFont(new java.awt.Font("Arial", 0, 14));
		moreButton.setText("+");
		moreButton.setContentAreaFilled(false);
		moreButton.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
		moreButton.setPreferredSize(new java.awt.Dimension(45, 27));
		moreButton.addActionListener(actListener);
		add(moreButton);
		
		removeButton = new JideButton("-" + "_" + index);
		removeButton.setText("-");
		removeButton.setActionCommand("-" + "_" + index);
		
		removeButton.setFont(new java.awt.Font("Arial", 0, 14));
		removeButton.setContentAreaFilled(false);
		removeButton
				.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
		removeButton.setPreferredSize(new java.awt.Dimension(45, 27));
		removeButton.addActionListener(actListener);
		add(removeButton);

	}

	public JPanel getDialog() {
		return newP;
	}

}
