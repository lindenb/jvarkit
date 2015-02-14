package com.github.lindenb.jvarkit.knime;

import org.knime.core.node.port.PortType;

public abstract class AbstractJVarkitNodeModel extends AbstractNodeModel
	{
	protected AbstractJVarkitNodeModel(int inport,int outport)
		{
		/* super(inport,outport) */
		super(inport,outport);
		}
	
	protected AbstractJVarkitNodeModel(PortType[] inPortTypes, PortType[] outPortTypes)
		{
		super(inPortTypes, outPortTypes);
		
		}
	}
