package org.broad.igv.sam;

import javax.swing.*;

import org.broad.igv.sam.AlignmentFilter.FilterType;

import java.awt.*;
import java.awt.event.ActionListener;
import java.text.NumberFormat;

public class AlignmentFilterPanel extends JPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 3328021930786326253L;

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

	public FilterType getFilterType() {
		if (missingButton.isSelected())
			return FilterType.missing;
		if (stringButton.isSelected())
			return FilterType.pattern;
		if (numberButton.isSelected())
			return FilterType.range;
		return null;
	}

	public AlignmentFilterPanel(AlignmentFilter alf, String index,
			ActionListener actListener) {
		super();
		Color bg=null;
		if (index.equals("new")) {
			bg = Color.LIGHT_GRAY;
		} else {
			bg = super.getBackground();
		}

		tagValue = new JTextField(4);
		JLabel tagFieldLabel = new JLabel("tag to be used:");
		tagFieldLabel.setLabelFor(tagValue);
		tagValue.setText(alf.getTag());
		

		// String or Number Filter
		numberButton = new JRadioButton(m_number + "_" + index);
		numberButton.setText(m_number);
		numberButton.addActionListener(actListener);
		numberButton
				.setSelected(alf.getFilterType() == AlignmentFilter.FilterType.range);
		numberButton.setBackground(bg);

		stringButton = new JRadioButton(m_string + "_" + index);
		stringButton.setText(m_string);
		stringButton.addActionListener(actListener);
		stringButton
				.setSelected(alf.getFilterType() == AlignmentFilter.FilterType.pattern);
		stringButton.setBackground(bg);

		missingButton = new JRadioButton(m_missing + "_" + index);
		missingButton.setText(m_missing);
		missingButton.addActionListener(actListener);
		missingButton
				.setSelected(alf.getFilterType() == AlignmentFilter.FilterType.missing);
		missingButton.setBackground(bg);

		ButtonGroup group = new ButtonGroup();
		group.add(numberButton);
		group.add(stringButton);
		group.add(missingButton);

		// Number parameters
		NumberFormat numFormat = NumberFormat.getNumberInstance();
		;

		minField = new JFormattedTextField(numFormat);
		minField.setText("Minimum Value");
		minField.setColumns(10);
		minField.addActionListener(actListener);
		
		JLabel minLabel = new JLabel(m_min);
		minLabel.setLabelFor(minField);
		minField.setValue(alf.getMin());

		maxField = new JFormattedTextField(numFormat);
		maxField.setColumns(10);
		maxField.setText("Maximum Value");
		maxField.setValue(0.0);
		JLabel maxLabel = new JLabel(m_max);
		maxLabel.setLabelFor(maxField);
		maxField.setValue(alf.getMax());
		

		JPanel numberButtonGroup = new JPanel(new FlowLayout(FlowLayout.LEFT));
		numberButtonGroup.add(minLabel);
		numberButtonGroup.add(minField);
		numberButtonGroup.add(maxLabel);
		numberButtonGroup.add(maxField);
		numberButtonGroup.setBackground(bg);
		
		// TODO up / down functionality

		// String parameters
		boolean stringEnabled = !alf.isHasWildCards() && !alf.isRegExpr();
		expressionString = new JTextField(20);
		JLabel expressionLabel = new JLabel(m_filterExpression + "_" + index);
		expressionLabel.setBackground(bg);
		expressionLabel.setLabelFor(expressionString);
		expressionString.setText(alf.getExpression());
		expressionString.setBackground(bg);
		
		hasWildCardsButton = new JRadioButton(m_hasWildCards + "_new");
		hasWildCardsButton.setText(m_hasWildCards);
		hasWildCardsButton.setEnabled(stringEnabled);
		hasWildCardsButton.addActionListener(actListener);
		hasWildCardsButton.setSelected(alf.isHasWildCards());
		hasWildCardsButton.setBackground(bg);

		isRegExprButton = new JRadioButton(m_isRegExprs + "_new");
		isRegExprButton.setText(m_isRegExprs);
		isRegExprButton.setEnabled(stringEnabled);
		isRegExprButton.addActionListener(actListener);
		isRegExprButton.setSelected(alf.isRegExpr());
		isRegExprButton.setBackground(bg);

		asStringButton = new JRadioButton(m_asString + "_" + index);
		asStringButton.setText(m_asString);
		asStringButton.setEnabled(stringEnabled);
		asStringButton.addActionListener(actListener);
		asStringButton.setSelected(stringEnabled);
		asStringButton.setBackground(bg);
		
		ButtonGroup stringGroup = new ButtonGroup();
		stringGroup.add(hasWildCardsButton);
		stringGroup.add(isRegExprButton);
		stringGroup.add(asStringButton);
		
		caseSensitiveBox = new JCheckBox(m_caseSensitive + "_" + index);
		caseSensitiveBox.setEnabled(false);
		caseSensitiveBox.setText(m_caseSensitive);
		caseSensitiveBox.addActionListener(actListener);
		caseSensitiveBox.setSelected(alf.isCaseSensitive());
		caseSensitiveBox.setBackground(bg);
		
		JButton button;
		JButton up = null;
		JButton down = null;
		JButton apply = null;
		GridLayout buttonLayout = new GridLayout(0, 2);
		JPanel applyRMGroup = new JPanel();
		applyRMGroup.setBackground(bg);
		applyRMGroup.setLayout(buttonLayout);
		if (index.equals("new")) {
			button = new JButton("create");
			applyRMGroup.add(button);
		} else {
			button = new JButton("remove");
			button.setName("remove" + "_" + index);
			up = new JButton("up");
			up.setName("up" + "_" + index);
			up.addActionListener(actListener);
			down = new JButton("down");
			down.setName("down" + "_" + index);
			down.addActionListener(actListener);
			apply = new JButton("apply");
			apply.setName("apply" + "_" + index);
			apply.addActionListener(actListener);
			applyRMGroup.add(up, 0);
			applyRMGroup.add(down, 1);
			applyRMGroup.add(apply, 2);
			applyRMGroup.add(button, 3);
		}
		button.addActionListener(actListener);

		// include exclude
		includeButton = new JRadioButton(m_include + "_" + index);
		includeButton.setText(m_include);
		includeButton.setSelected(!alf.isExclude());
		includeButton.setBackground(bg);
		excludeButton = new JRadioButton(m_exclude + "_" + index);
		excludeButton.setText(m_exclude);
		excludeButton.setSelected(alf.isExclude());
		excludeButton.setBackground(bg);

		JPanel inclGroup = new JPanel(new FlowLayout(FlowLayout.LEFT));
		inclGroup.add(includeButton);
		inclGroup.add(excludeButton);
		inclGroup.setBackground(bg);

		ButtonGroup inclButtonGroup = new ButtonGroup();
		inclButtonGroup.add(includeButton);
		inclButtonGroup.add(excludeButton);

		// JPanel inclExclPanel = new JPanel(new
		// FlowLayout(FlowLayout.TRAILING));
		// inclExclPanel.add(includeButton);
		// inclExclPanel.add(excludeButton);

		JPanel stringParameterGroup = new JPanel(new GridLayout(0, 1));
		JPanel stringButtonGroup = new JPanel(new FlowLayout(FlowLayout.LEFT));
		stringButtonGroup.add(hasWildCardsButton);
		stringButtonGroup.add(isRegExprButton);
		stringButtonGroup.add(asStringButton);
		stringButtonGroup.setBackground(bg);
		
		stringParameterGroup.add(expressionString);
		stringParameterGroup.add(stringButtonGroup);
		stringParameterGroup.add(caseSensitiveBox);
		stringParameterGroup.setBackground(bg);
		
		FlowLayout alignmentFilterLayout = new FlowLayout();
		newP = new JPanel();
		newP.setLayout(new BoxLayout(newP, BoxLayout.Y_AXIS));

		alignmentFilterLayout.setAlignment(FlowLayout.LEFT);

		newP.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createTitledBorder(index),
				BorderFactory.createEmptyBorder(1, 1, 1, 1)));

		JPanel tagPanel = new JPanel(new FlowLayout());
		tagPanel.setAlignmentX(FlowLayout.LEFT);
		tagPanel.add(tagFieldLabel);
		tagPanel.add(tagValue);
		tagPanel.setBackground(bg);
		newP.add(tagPanel);

		// JPanel choicePanel = new JPanel();
		// choicePanel.setLayout(new BoxLayout(choicePanel,
		// BoxLayout.Y_AXIS));

		JPanel numberChoice = new JPanel();
		numberChoice.setBackground(bg);
		GroupLayout glayout = new GroupLayout(numberChoice);
		numberChoice.setLayout(glayout);
		// glayout.setAutoCreateContainerGaps(true);
		// glayout.setAutoCreateGaps(true);
		glayout.setHorizontalGroup(glayout.createParallelGroup(
				GroupLayout.Alignment.LEADING).addGroup(
				glayout.createParallelGroup(GroupLayout.Alignment.LEADING)
						.addComponent(numberButton, 0,
								GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
						.addComponent(numberButtonGroup, 0,
								GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
						.addComponent(stringButton, 0,
								GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
						.addComponent(stringButtonGroup, 0,
								GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
						.addComponent(missingButton, 0,
								GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)));
		glayout.setVerticalGroup(glayout
				.createSequentialGroup()
				.addGroup(
						glayout.createSequentialGroup()
								.addComponent(numberButton, 0,
										GroupLayout.DEFAULT_SIZE,
										Short.MAX_VALUE)
								.addPreferredGap(
										LayoutStyle.ComponentPlacement.RELATED)
								.addComponent(numberButtonGroup, 0,
										GroupLayout.DEFAULT_SIZE,
										Short.MAX_VALUE))
				.addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
				.addGroup(
						glayout.createSequentialGroup()
								.addComponent(stringButton, 0,
										GroupLayout.DEFAULT_SIZE,
										Short.MAX_VALUE)
								.addComponent(stringButtonGroup, 0,
										GroupLayout.DEFAULT_SIZE,
										Short.MAX_VALUE))
				.addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
				.addComponent(missingButton));

		newP.add(numberChoice);
		JPanel lastGroup = new JPanel();
		lastGroup.setBackground(bg);
		lastGroup.setLayout(new BoxLayout(lastGroup, BoxLayout.X_AXIS));
		button.setAlignmentX(Component.RIGHT_ALIGNMENT);
		inclGroup.setAlignmentX(Component.LEFT_ALIGNMENT);
		inclGroup.setBackground(bg);
		lastGroup.add(inclGroup);
		lastGroup.add(Box.createHorizontalGlue());
		lastGroup.add(applyRMGroup);
		newP.add(lastGroup);
		if (index.equals("new")) {
			newP.setBackground(Color.LIGHT_GRAY);
		}
		this.add(newP);

	}

	public JPanel getDialog() {
		return newP;
	}

}
