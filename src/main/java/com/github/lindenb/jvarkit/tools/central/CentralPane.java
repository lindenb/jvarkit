package com.github.lindenb.jvarkit.tools.central;

import java.awt.BorderLayout;
import java.awt.Component;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

@SuppressWarnings("serial")
public abstract class CentralPane extends JPanel
	{
	private Map<String, String> externalResourceMap = Collections.emptyMap();
	
	public CentralPane()
		{
		super(new BorderLayout(5,5));
		this.setBorder(new EmptyBorder(5, 5, 5, 5));
		}
	public abstract Class<?> getMainClass();
	public abstract String getLabel();
	public abstract String getDescription();
	public String getValidationMessage()
		{
		return null;
		}
	
	public void fillCommandLine(List<String> command)
		{
		
		}
	
	protected String makeLabel(String label)
		{
		return label;
		}
	protected JLabel makeLeftLabel(String label)
		{
		JLabel lbl= new JLabel(makeLabel(label)+":",JLabel.RIGHT);
		
		return lbl;
		}
	
	public void setExternalResourceMap(Map<String, String> externalResourceMap) {
		this.externalResourceMap = externalResourceMap;
		}
	public Map<String, String> getExternalResourceMap() {
		return externalResourceMap;
		}
	}
