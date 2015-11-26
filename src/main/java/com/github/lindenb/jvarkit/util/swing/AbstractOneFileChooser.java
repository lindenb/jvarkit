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
	public static final String EDITABLE="central.file.editable";
	private JTextField textField;
	private AbstractAction setAction;
	private AbstractAction resetAction;
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
	private File chooseFile()
		{
		
		final JFileChooser fc = new JFileChooser(PreferredDirectory.get(getClass()));
		fc.setFileFilter(getFileFilter());
		if(showFileChooser(fc)!=JFileChooser.APPROVE_OPTION) return null;
		File f= fc.getSelectedFile();
		PreferredDirectory.update(f.getParentFile());
		return f;
		}
	protected abstract int showFileChooser(final JFileChooser fc);
	}
