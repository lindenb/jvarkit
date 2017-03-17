/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
* 2016 creation

*/
package com.github.lindenb.jvarkit.util.command;


public abstract class Argument<T>
	{
	protected String shortopt=null;
	protected String longOpt=null;
	protected String description=null;
	
	public abstract static class AbstractBuilder<Y,A extends Argument<Y>,B extends AbstractBuilder<Y,A,B>>
		{
		protected final A instance;
		protected AbstractBuilder(final A instance)
			{
			this.instance=instance;
			}
		// http://stackoverflow.com/questions/5818504/
		protected abstract B getThis();
		
		protected void validate()
			{
			
			}
		protected B desc(final String s) { this.instance.description=s; return getThis();}
		protected B shortOpt(final char c) { return this.shortOpt(String.valueOf(c));}
		protected B shortOpt(final String s) { this.instance.shortopt = s; return getThis();}
		protected B longOpt(final String s) { this.instance.longOpt = s; return getThis();}
		public A make() { validate();return this.instance;}
		}
	
	}
