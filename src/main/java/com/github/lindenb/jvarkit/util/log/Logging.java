package com.github.lindenb.jvarkit.util.log;

import org.slf4j.Logger;

@Deprecated
public class Logging
	{
	public static Logger getLog(Class<?> clazz)
		{
		return getLog("jvarkit");
		}
	public static Logger getLog(final String name)
		{
		String logName="jvarkit";
		try {
			final String s = System.getProperty("jvarkit.log.name",logName);
			logName=s;
		} catch (final java.security.AccessControlException e) {
			
			}
		
		return org.slf4j.LoggerFactory.getLogger(logName);
		}
	}
