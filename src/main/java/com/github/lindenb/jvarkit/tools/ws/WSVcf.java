package com.github.lindenb.jvarkit.tools.ws;

import java.io.Serializable;

public interface WSVcf extends Serializable
	{
	public String getId();
	public String getReferenceId();
	public String getPath();
	}
