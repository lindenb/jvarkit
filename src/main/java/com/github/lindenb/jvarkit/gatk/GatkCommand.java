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

*/
package com.github.lindenb.jvarkit.gatk;

import java.util.List;

public abstract class GatkCommand {
	protected static class OptionArg
		{
		final String name;
		final String summary;
		public OptionArg(final String name,final String summary) {
			this.name = name;
			this.summary = summary;
			}
		}
	protected abstract class GatkArg
		{
		protected abstract String name();
		protected abstract String synonyms();
		protected abstract String defaultValue();
		protected abstract String type();
		protected abstract String minValue();	
		protected abstract String maxValue();		
		protected abstract String minRecValue();	
		protected abstract String maxRecValue();		
		protected abstract String rodTypes();						
		protected abstract String kind();
		protected abstract String summary();
		protected abstract String required();
		protected abstract OptionArg[] options();
		
		public String getId() {
			return getGatkCommand().getName()+"."+getName();
		}
		
		public String getName() {
			return name().replaceAll("\\-", "");
			}
		public String getSummary() {
			return summary();
			}
		public String getSynonyms() {
			return synonyms();
			}
		
		public boolean isRequired() {
			return !required().equals("no");
			}
		public GatkCommand getGatkCommand() { return GatkCommand.this;}
		}
	protected GatkArg _wrapArg(final GatkArg arg) {
		return arg;
	}
	protected boolean _acceptArg(final GatkArg arg) {
		return true;
		}
	protected List<GatkArg> fillArguments(final List<GatkArg> args) {
		return args;
		}
	public abstract String getName();
	public abstract String getSummary();
	protected abstract List<GatkArg> getArguments();
}
