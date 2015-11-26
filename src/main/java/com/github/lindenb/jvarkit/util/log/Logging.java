package com.github.lindenb.jvarkit.util.log;

import org.slf4j.Logger;

public class Logging
	{
	public static Logger getLog(Class<?> clazz)
		{
		return getLog("jvarkit");
		}
	public static Logger getLog(final String name)
		{
		return org.slf4j.LoggerFactory.getLogger("jvarkit");
		}
	}
