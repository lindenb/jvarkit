package com.github.lindenb.jvarkit.util.swing;

import javax.swing.JFileChooser;

@SuppressWarnings("serial")
public class InputChooser extends AbstractOneFileChooser {

	public InputChooser()
		{
		}
	
	@Override
	protected int showFileChooser(final JFileChooser fc)
		{
		return fc.showOpenDialog(this);
		}
	}
