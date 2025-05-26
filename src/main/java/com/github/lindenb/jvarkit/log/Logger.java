/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.log;

import java.io.PrintStream;
import java.util.Objects;


/** Generic & simple Logging class */
public interface Logger  {
public enum Level {DEBUG,INFO,WARN,SEVERE,FATAL}


public static Logger of(Class<?> clazz) {
	final DefaultLogger logger = new DefaultLogger(clazz);
	return logger;
	}

public default Logger fine(Object s) {
	return info(s);
	}


public default Logger info(Object s) {
	return log(Level.INFO,s);
	}
public default Logger debug(Object s) {
	return log(Level.DEBUG,s);
	}
public default Logger severe(final Object err) {
	return log(Level.SEVERE,err);
	}

public default Logger severe(final Object err,Throwable t) {
	return log(Level.SEVERE,err,t);
	}


public default Logger error(final Object err) {
	return severe(err);
	}
public default Logger error(final Object err,Throwable t) {
	return severe(err,t);
	}


public default Logger warning(final Object err) {
	return warn(err);
	}

public default Logger warn(final Object err) {
	return log(Level.WARN,err);
	}

public default Logger fatal(final Object err) {
	return log(Level.FATAL,err);
	}
public default Logger log(Logger.Level level,final Object o) {
	if(o==null)
		{
		return this.log(level, "null",null);
		}
	else if(o instanceof Throwable)
		{
		final Throwable t=Throwable.class.cast(o);
		return this.log(level, String.valueOf(t.getMessage()), t);
		}
	else
		{
		return this.log(level, String.valueOf(o), null);
		}
	}	
public  Logger log(Logger.Level level,final Object msg,final Throwable err);	

public Logger setLevel(Level level);
public default Logger setLevel(String level) {
	return setLevel(Level.valueOf(level));
}
public Level getLevel();

public default Logger setDebug() {
	return setLevel(Level.DEBUG);
}

public default boolean isDebug() {
	return getLevel().compareTo(Level.DEBUG)<=0;
	}
public default boolean isInfo() {
	return getLevel().compareTo(Level.INFO)<=0;
	}
public default boolean isWarning() {
	return getLevel().compareTo(Level.WARN)<=0;
	}
public default boolean isFatal() {
	return getLevel().compareTo(Level.FATAL)<=0;
	}
static class DefaultLogger implements Logger
	{
	private PrintStream out=System.err;
	private final Class<?> clazz;
	private Level level = Level.INFO;
	DefaultLogger(Class<?> clazz) {
		this.clazz = Objects.requireNonNull(clazz);
		final String property = "log."+clazz.getName();
		final String value=java.lang.System.getProperty(property);
		if(value!=null) {
			try {
				this.level = Level.valueOf(value);
				}
			catch(Throwable err) {
				System.err.println("[WARN] not a valid value for -D"+property+"="+value);
				}
			}
		}
	@Override
	public Logger setLevel(Level level) {
		this.level=level;
		return this;
		}
	@Override
	public Level getLevel() {
		return this.level;
		}
	
	
	@Override
	public Logger log(final Logger.Level level,final Object msg,final Throwable err) {
		//if(level.compareTo(this.getLevel()) <=0) {
			this.out.print("[");
			this.out.print(clazz.getName());
			this.out.print(":");
			this.out.print(this.getLevel().name());
			this.out.print("]");
			if(msg!=null)
				{
				this.out.print(String.valueOf(msg));
				}
			this.out.println();
			if(err!=null)
				{
				if( err instanceof java.net.ConnectException) {
					this.log(level,"It could be a proxy error. Did you set the java proxy ? See http://docs.oracle.com/javase/6/docs/technotes/guides/net/proxies.html. "+
							"Something like `java -Dhttp.proxyHost=myproxyhost  -Dhttp.proxyPort=myproxyport -Dhttps.proxyHost=myproxyhost  -Dhttps.proxyPort=myproxyport  `",
							null);
					}
				err.printStackTrace(this.out);
				}
		//	}
		return this;
		}	
	

	}


}
