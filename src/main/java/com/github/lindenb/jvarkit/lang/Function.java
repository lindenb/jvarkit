package com.github.lindenb.jvarkit.lang;

/** replace with java 8 functional please */
@Deprecated
public interface Function<PARAM, RETURN>
	{
	public RETURN apply(PARAM param);
	}
