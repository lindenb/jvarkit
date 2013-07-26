package com.github.lindenb.jvarkit.tools.ws;

import java.io.Serializable;

public interface WSReference extends Serializable
	{
	public String getId();
	public String getLabel();
	public String getDescription();
	}
