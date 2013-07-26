package com.github.lindenb.jvarkit.tools.ws;

import java.io.Serializable;
import java.util.Set;

public interface WSProject extends Serializable
	{
	public String getId();
	public String getLabel();
	public String getDescription();
	public Set<String> getBamIds();
	public Set<String> getVcfIds();
	}
