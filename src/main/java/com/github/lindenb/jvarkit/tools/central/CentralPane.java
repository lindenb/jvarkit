package com.github.lindenb.jvarkit.tools.central;

import java.awt.BorderLayout;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;


@SuppressWarnings("serial")
public abstract class CentralPane extends JPanel
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(CentralPane.class);
	private Map<String, String> externalResourceMap = Collections.emptyMap();
	
	public CentralPane()
		{
		super(new BorderLayout(5,5));
		this.setBorder(new EmptyBorder(5, 5, 5, 5));
		}
	public abstract Class<? extends com.github.lindenb.jvarkit.util.command.Command> getMainClass();
	public abstract String getLabel();
	public abstract String getDescription();
	public String getValidationMessage()
		{
		return null;
		}
	
	public void fillCommandLine(List<String> command)
		{
		
		}
	protected String makeLabel(String s)
		{
		while(s.startsWith("-")) s=s.substring(1);
		for(int i=0;i+1<s.length();++i)
			{
			if( Character.isLetter(s.charAt(i)) &&
				Character.isLetter(s.charAt(i+1)) &&
				Character.isLowerCase(s.charAt(i)) &&
				Character.isUpperCase(s.charAt(i+1)))
				{
				return makeLabel(s.substring(0,i+1)+" "+
						Character.toLowerCase(s.charAt(i+1))
						+(i+2<s.length()?s.substring(i+2):"")
						);
				}
			}
		return s.replace('_', ' ');
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
