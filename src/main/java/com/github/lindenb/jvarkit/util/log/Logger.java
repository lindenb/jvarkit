/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
	return severe(err,null);
	}

public Logger severe(final Object err,Throwable t) {
	return log(Level.SEVERE,err,t);
	}


public Logger error(final Object err) {
	return severe(err);
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
public Logger log(Logger.Level level,final Object msg,final Throwable err) {
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
		err.printStackTrace(this.out);
		}
	
	return this;
	}	

public static Builder build() {
	return new Builder();
	}

public static Builder build(final Class<?> c) {
	final Builder b=build();
	if(c!=null) {
		final String s=c.getSimpleName();
		if(!(s==null || s.isEmpty() || s.contains("[")))
			{
			b.prefix(s);
			}
		}
	return b;
	}

}
