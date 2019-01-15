/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.log;

import java.io.PrintStream;
import java.lang.reflect.Method;

/** Generic & simple Logging class */
public abstract class Logger  {
public enum Level {DEBUG,INFO,WARN,SEVERE,FATAL}
protected PrintStream out=System.err;
protected String prefix=null;
private static boolean TRY_KNIME = true;

public static class Builder
	{
	private Class<?> clazz=null;
	private PrintStream out=System.err;
	private String prefix="[LOG]";
	public Builder output(final PrintStream out) {
		this.out = out;
		return this;
		}
	public Builder prefix(final String prefix) {
		this.prefix = prefix;
		return this;
		}
	public Builder setClass(final Class<?> clazz) {
		this.clazz = clazz;
		return this;
		}
	public Logger make() {
		if(this.clazz!=null && TRY_KNIME)
			{
			try
				{
				final KnimeLogger K=new KnimeLogger(this.clazz);
				K.prefix = prefix;
				return K;
				}
			catch(final Throwable err) {
				TRY_KNIME=false;
				}
			}
		
		
		final DefaultLogger L= new DefaultLogger();
		L.out=this.out==null?System.err:this.out;
		L.prefix=this.prefix;

		if(this.clazz!=null) {
			final String s= this.clazz.getSimpleName();
			if(!(s==null || s.isEmpty() || s.contains("[")))
				{
				L.prefix = s;
				}
			}
		
		return L;
		}
	}
private Logger() {
}

public Logger fine(Object s) {
	return info(s);
	}


public Logger info(Object s) {
	return log(Level.INFO,s);
	}
public Logger debug(Object s) {
	return log(Level.DEBUG,s);
	}
public Logger severe(final Object err) {
	return log(Level.SEVERE,err);
	}

public Logger severe(final Object err,Throwable t) {
	return log(Level.SEVERE,err,t);
	}


public Logger error(final Object err) {
	return severe(err);
	}
public Logger error(final Object err,Throwable t) {
	return severe(err,t);
	}


public Logger warning(final Object err) {
	return warn(err);
	}

public Logger warn(final Object err) {
	return log(Level.WARN,err);
	}

public Logger fatal(final Object err) {
	return log(Level.FATAL,err);
	}
public Logger log(Logger.Level level,final Object o) {
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
public abstract Logger log(Logger.Level level,final Object msg,final Throwable err);	

public static Builder build() {
	return new Builder();
	}

public static Builder build(final Class<?> c) {
	final Builder b=build().setClass(c);
	return b;
	}

protected void andThen(final Logger.Level level,Throwable err) {
	if(err==null) return;
	if( err instanceof java.net.ConnectException) {
		this.log(level,"It could be a proxy error. Did you set the java proxy ? See http://docs.oracle.com/javase/6/docs/technotes/guides/net/proxies.html. "+
				"Something like `java -Dhttp.proxyHost=myproxyhost  -Dhttp.proxyPort=myproxyport -Dhttps.proxyHost=myproxyhost  -Dhttps.proxyPort=myproxyport  `",
				null);
		}
}


private static class DefaultLogger
	extends Logger
	{
	@Override
	public Logger log(final Logger.Level level,final Object msg,final Throwable err) {
		this.out.print("[");
		this.out.print(level.name());
		this.out.print("]");
		this.out.print("[");
		this.out.print(prefix);
		this.out.print("]");
		if(msg!=null)
			{
			this.out.print(String.valueOf(msg));
			}
		this.out.println();
		if(err!=null)
			{
			andThen(level,err);
			err.printStackTrace(this.out);
			}
		
		return this;
		}	

	}

private static class KnimeLogger
	extends Logger
	{
	private final Object nodeLogger;
	
	private final Method info;
	private final Method debug;
	private final Method warn;
	private final Method error;
	private final Method fatal;
	
	
	private final Method info_t;
	private final Method debug_t;
	private final Method warn_t;
	private final Method error_t;
	private final Method fatal_t;
	
	KnimeLogger(final Class<?> clazz) throws Exception {
		final Class<?> org_knime_core_node_NodeLogger_clazz = Class.forName("org.knime.core.node.NodeLogger"); 
	    final Method getLogger = org_knime_core_node_NodeLogger_clazz.getMethod("getLogger",Class.class);
	    this.nodeLogger = getLogger.invoke(null, clazz);
	    this.info = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class);
	    this.debug = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class);
	    this.warn = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class);
	    this.error = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class);
	    this.fatal = org_knime_core_node_NodeLogger_clazz.getMethod("fatal", Object.class);

	    this.info_t = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class,Throwable.class);
	    this.debug_t = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class,Throwable.class);
	    this.warn_t = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class,Throwable.class);
	    this.error_t = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class,Throwable.class);
	    this.fatal_t = org_knime_core_node_NodeLogger_clazz.getMethod("warn", Object.class,Throwable.class);
		}
	@Override
	public Logger log(final Level level,final Object msg,final Throwable err) {
		try {
			if(err!=null)
				{
				switch(level)
					{
					case SEVERE: this.error_t.invoke(this.nodeLogger,msg, err);break;
					case FATAL: this.fatal_t.invoke(this.nodeLogger,msg, err);break;
					case WARN : this.warn_t.invoke(this.nodeLogger,msg, err);break;
					case DEBUG : this.debug_t.invoke(this.nodeLogger,msg, err);break;
					default : this.info_t.invoke(this.nodeLogger,msg, err);break;	
					}
				} 
			else
				{
				switch(level)
					{
					case SEVERE: this.error.invoke(this.nodeLogger,msg);break;
					case FATAL: this.fatal.invoke(this.nodeLogger,msg);break;
					case WARN : this.warn.invoke(this.nodeLogger,msg);break;
					case DEBUG : this.debug.invoke(this.nodeLogger,msg);break;
					default : this.info.invoke(this.nodeLogger,msg);break;	
					}
				}
			andThen(level,err);
		} catch(final Exception exc)
			{
			System.err.println("KNIME Logger broken "+exc.getMessage());
			System.err.println(msg);
			}
		return this;
		}
	}

}
