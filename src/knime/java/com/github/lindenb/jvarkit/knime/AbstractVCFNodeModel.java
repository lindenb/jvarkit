package com.github.lindenb.jvarkit.knime;

import org.knime.core.node.port.PortType;

public abstract class AbstractVCFNodeModel extends AbstractJVarkitNodeModel
	{
	protected AbstractVCFNodeModel(int inport,int outport)
		{
		/* super(inport,outport) */
		super(inport,outport);
		}
	
	protected AbstractVCFNodeModel(PortType[] inPortTypes, PortType[] outPortTypes)
		{
		super(inPortTypes, outPortTypes);
		}
	
	}
