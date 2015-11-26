package com.github.lindenb.jvarkit.util.swing;

import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.filechooser.FileFilter;

@SuppressWarnings("serial")
public abstract class AbstractFileChooser extends JPanel{
	public static final String FILTER_CHANGED="central.file.filter.changed";
	private FileFilter fileFilter;
	protected AbstractFileChooser()
		{
		super(new BorderLayout());
		}
	public void setFileFilter(final FileFilter fileFilter)
		{
		FileFilter old= this.fileFilter;
		this.fileFilter = fileFilter;
		firePropertyChange(FILTER_CHANGED, old, fileFilter);
		}
	public FileFilter getFileFilter()
		{
		return fileFilter;
		}
	}
