package com.github.lindenb.jvarkit.util.log;

import java.io.PrintStream;

public class Logger  {
public enum Level {DEBUG,INFO,WARN,SEVERE,FATAL}
private PrintStream out=System.err;
private String prefix=null;

public static class Builder
	{
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
	public Logger make() {
		final Logger L= new Logger();
		L.out=this.out==null?System.err:this.out;
		L.prefix=this.prefix;
		return L;
		}
	}
private Logger() {
}
public Logger info(Object s) {
	return log(Level.INFO,s);
	}

public Logger severe(final Object err) {
	return log(Level.FATAL,err);
	}

public Logger fatal(final Throwable err) {
	return log(Level.FATAL,err.getMessage(),err);
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
public Logger log(Logger.Level level,String msg,Throwable err) {
	this.out.print("[");
	this.out.print(level.name());
	this.out.print("]");
	this.out.print("[");
	this.out.print(prefix);
	this.out.print("]");
	if(msg!=null)
		{
		this.out.print(msg);
		}
	this.out.println();
	if(err!=null)
		{
		err.printStackTrace(this.out);
		}
	
	return this;
	}	

public static Builder build() {
	return new Builder();
	}
}
