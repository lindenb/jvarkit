package com.github.lindenb.jvarkit.util.swing;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

@SuppressWarnings("serial")
public class OutputChooser extends AbstractOneFileChooser {

	public OutputChooser()
		{
		}
	
	@Override
	protected int showFileChooser(final JFileChooser fc)
		{
		final int ret= fc.showSaveDialog(this);
		if(ret==JFileChooser.APPROVE_OPTION)
			{
			final File f=fc.getSelectedFile();
			if(f.exists() && JOptionPane.showConfirmDialog(this, "File "+f.getName()+" exists. Overwite ?", "Overwite ?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
				{
				return 	JFileChooser.CANCEL_OPTION ;
				}
				
			}
		return ret;
		}
	}
