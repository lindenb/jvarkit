package com.github.lindenb.jvarkit.tools.ws.rmi.client;

import java.util.List;

import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.SwingUtilities;
import javax.swing.table.AbstractTableModel;

import com.github.lindenb.jvarkit.tools.ws.WSProject;
import com.github.lindenb.jvarkit.tools.ws.rmi.NGSService;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.swing.AbstractGenericTable;

class ProjectsTable extends AbstractGenericTable<WSProject>
	{
	private static String COLS[]={"Label","Description"};
	@Override
	public int getColumnCount()
		{
		return COLS.length;
		}
	@Override
	public Object getValueOf(WSProject o, int columnIndex)
		{
		switch (columnIndex)
			{
			case 0: return o.getLabel();
			case 1: return o.getDescription();
			}
		return null;
		}
	}
class ProjectFrame extends AbstractInternalFrame
	{
	private WSProject project;
	ProjectFrame(NGSServiceClientFrame owner,WSProject project)
		{
		super(owner);
		this.project=project;
		setTitle("Project: "+project.getLabel());
		for(String s:project.getBamIds())
			{
			
			}
		for(String s:project.getVcfIds())
			{
			
			}
		}
	}

class ProjectsFrame extends AbstractTableFrame<WSProject>
	{
	ProjectsFrame(NGSServiceClientFrame owner)
		{
		super(owner);
		}
	
	
	private void x()
		{
		new Thread()
			{
			public void run()
					{
					try
						{
						final List<? extends WSProject> L=getService().getProjects();
						SwingUtilities.invokeAndWait(new Runnable()
							{
							@Override
							public void run()
								{
								
								}
							});
						}
					catch (Exception e) {
						
						}
					}
			}.start();
		}
	}

abstract class AbstractTableFrame<T> extends AbstractInternalFrame
	{
	AbstractTableFrame(NGSServiceClientFrame owner)
		{
		super(owner);
		}
	}

abstract class AbstractInternalFrame
extends JInternalFrame
	{
	private NGSServiceClientFrame owner;
	
	AbstractInternalFrame(NGSServiceClientFrame owner)
		{
		this.owner=owner;
		}
	
	NGSServiceClientFrame getOwnerWindow()
		{
		return this.owner;
		}
	
	NGSService getService()
		{
		return getOwnerWindow().getService();
		}
	}


class NGSServiceClientFrame extends JFrame
	{
	private NGSService service;
	JDesktopPane desktop;
	
	
	NGSServiceClientFrame(NGSService service)
		{
		this.service=service;
		
		}
	NGSService getService()
		{
		return this.service;
		}
	}



public class NGSServiceClient extends AbstractCommandLineProgram
	{
	}
