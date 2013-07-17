package com.github.lindenb.jvarkit.lang;


public interface Function<PARAM, RETURN>
	{
	public RETURN apply(PARAM param);
	}
