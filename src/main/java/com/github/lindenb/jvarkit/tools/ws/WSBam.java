package com.github.lindenb.jvarkit.tools.ws;

import java.io.Serializable;


public interface WSBam extends Serializable
	{
	public String getId();
	public String getReferenceId();
	public String getPath();
	}
