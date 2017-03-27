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
public Logger info(String s) {
	return log(Level.INFO,s,null);
	}
public Logger fatal(final Throwable err) {
	return log(Level.FATAL,err.getMessage(),err);
	}

public Logger log(Logger.Level level,String msg,Throwable err) {
	return this;
	}	

public static Builder build() {
	return new Builder();
	}
}
