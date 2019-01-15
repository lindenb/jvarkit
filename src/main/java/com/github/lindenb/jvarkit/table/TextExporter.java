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
package com.github.lindenb.jvarkit.table;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;


import htsjdk.samtools.util.RuntimeIOException;

public class TextExporter extends AbstractTableExporter {
private String margin="";
private boolean allowUnicode = false;

public void setAllowUnicode(boolean allowUnicode) {
	this.allowUnicode = allowUnicode;
	}

public boolean isAllowUnicode() {
	return allowUnicode;
	}

public void print(final Table table,final Appendable out) {
	if(isRemoveEmptyColumns()) table.removeEmptyColumns();
	if(!isPrintEmptyTable() && table.isEmpty()) return;
	//https://en.wikipedia.org/wiki/Box-drawing_character
	final int lengths[]= table.getColumns().
			stream().
			mapToInt(C->C.getName().length()).
			toArray();
	
	for(int x=0;x<table.getColumnCount();++x) {
		final Column c = table.getColumn(x);
		for(int y=0;y<table.getRowCount();++y) {
			final int L = c.toString(table.at(y, x)).length();
			lengths[x] = Math.max(lengths[x], L);
			}
		}
	
	
	try {
		//print header
		
		// line 1 of header
		for(int x=0;x<table.getColumnCount();++x) {
			//final Column col = table.getColumn(x);
			out.append(isAllowUnicode()?(x==0?'\u250C':'\u252C'):'+');
			char c = isAllowUnicode()?'\u2500':'-';
			out.append(c);
			repeat(out,lengths[x],c);
			out.append(c);
			}
		out.append(isAllowUnicode()?'\u2510':'+');
		out.append('\n');
		
		//line 2 of header
		for(int x=0;x<table.getColumnCount();++x) {
			//final Column col = table.getColumn(x);
			out.append(isAllowUnicode()?'\u2502':'|');
			out.append(' ');
			out.append(table.getColumn(x).getName());
			repeat(out,lengths[x]-table.getColumn(x).getName().length(),' ');
			out.append(' ');
			}
		out.append(isAllowUnicode()?'\u2502':'|');
		out.append('\n');
		
		//line 3 of header
		for(int x=0;x<table.getColumnCount();++x) {
			if(!isAllowUnicode())
				{
				out.append('+');
				}
			else if(x==0 && table.getRowCount()==0)
				{	
				out.append('\u2514');
				}
			else if(x==0)
				{
				out.append('\u251C');
				}
			else
				{
				out.append('\u2500');
				}
			char c = isAllowUnicode()?'\u2500':'-';
			out.append(c);
			repeat(out,lengths[x],c);
			out.append(c);
			}
		if(!isAllowUnicode())
			{
			out.append('+');
			}
		else if(table.getRowCount()==0)
			{	
			out.append('\u2518');
			}
		else
			{
			out.append('\u2524');
			}
		out.append('\n');
		
		//print body
		for(int y=0;y< table.getRowCount();++y) {
			final List<Object> row = table.getRow(y);
			//line  of data
			for(int x=0;x<table.getColumnCount();++x) {
				final String str = table.getColumn(x).toString(row.get(x));
				//final Column col = table.getColumn(x);
				out.append(isAllowUnicode()?'\u2502':'|');
				out.append(' ');
				out.append(str);
				repeat(out,lengths[x]-str.length(),' ');
				out.append(' ');
				}
			out.append(isAllowUnicode()?'\u2502':'|');
			out.append('\n');
			}
		//last line
		if(table.getRowCount()>0)
			{
			for(int x=0;x<table.getColumnCount();++x) {
				if(!isAllowUnicode())
					{
					out.append('+');
					}
				else if(x==0) {
					out.append('\u2514');
					}
				else 
					{
					out.append('\u2534');
					}
				char c = isAllowUnicode()?'\u2500':'-';
				out.append(c);
				repeat(out,lengths[x],c);
				out.append(c);
				}
			out.append(isAllowUnicode()?'\u2518':'+');
			out.append('\n');
			}
		
		}
	catch(final IOException err) {
		throw new RuntimeIOException(err);
		}
	}


private void repeat(final Appendable a,int n,char c) throws IOException {
	while(n>0) {
		a.append(c);
		--n;
		}
	}

@Override
public  void saveTableTo(final Table table,final PrintWriter p) throws IOException {
	try(final PrintWriter pw=new PrintWriter(p)) {
		this.print(table,pw);
		pw.flush();
		pw.close();
		}
	}
}
