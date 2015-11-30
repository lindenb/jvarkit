package com.github.lindenb.jvarkit.util.swing;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JTextField;

 @SuppressWarnings("serial")
public abstract class AbstractOneFileChooser extends AbstractFileChooser
	{
	public enum SelectType {SELECT_FILE,SELECT_DIRECTORY};
	public static final String EDITABLE="central.file.editable";
	public static final String TYPEOFFILE="central.file.typeoffile";
	private JTextField textField;
	private AbstractAction setAction;
	private AbstractAction resetAction;
	private SelectType selectType = SelectType.SELECT_FILE;
	protected  AbstractOneFileChooser()
		{
		setLayout(new FlowLayout(FlowLayout.LEADING));
		this.textField = new JTextField(20);
		this.add(this.textField);
		this.textField.setEditable(false);
		this.setAction = new AbstractAction("[+]")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				File f = chooseFile();
				if(f==null) return;
				textField.setText(f.getPath());
				}
			};
		this.setAction.putValue(AbstractAction.LONG_DESCRIPTION, "Set File...");
		this.add(new JButton(this.setAction));
		this.resetAction = new AbstractAction("[-]")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				textField.setText("");
				}
			};
		this.resetAction.putValue(AbstractAction.LONG_DESCRIPTION, "Delete File...");
		this.add(new JButton(this.resetAction));
		}
	
	public boolean isEmpty()
		{
		return getText().trim().isEmpty();
		}
	
	public String getText()
		{
		return getTextField().getText();
		}
	
	public JTextField getTextField()
		{
		return textField;
		}
	
	public void setSelectType(SelectType selectType) {
		SelectType old = this.selectType;
		this.selectType = selectType;
		firePropertyChange(TYPEOFFILE, old, selectType);
		}
	
	public SelectType getSelectType() {
		return selectType;
	}

	
	
	private File chooseFile()
		{
		
		final JFileChooser fc = new JFileChooser(PreferredDirectory.get(getClass()));
		switch(getSelectType())
			{
			case SELECT_DIRECTORY: fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			default:  fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
			}
		fc.setFileFilter(getFileFilter());
		if(showFileChooser(fc)!=JFileChooser.APPROVE_OPTION) return null;
		File f= fc.getSelectedFile();
		PreferredDirectory.update(f.getParentFile());
		return f;
		}
	protected abstract int showFileChooser(final JFileChooser fc);
	}
