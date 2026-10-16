/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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

import java.io.Closeable;
import java.io.Flushable;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.util.RuntimeIOException;

public class TableWriterFactory  {

	public interface TableWriter extends Closeable,Flushable {
	public  List<String> getColumns();
	public void write(Map<String,Object> hash);
	}

 static abstract class AbstractTableWriter implements TableWriter {
	private final List<String> columns;
	AbstractTableWriter(final List<String> columns) {
		this.columns = Collections.unmodifiableList(columns);
		}
	
	String objectToString(Object o) {
		return String.valueOf(o);
		}
	
	@Override
	public List<String> getColumns() {
		return this.columns;
		}	
	}
 
 static class JsonTableWriter extends AbstractTableWriter {
	final Writer out;
	final JsonWriter w;
	JsonTableWriter(Path p,final List<String> columns) throws IOException {
		super(columns);
		this.out=IOUtils.openPathForPrintWriter(p);
		w= new JsonWriter(out);
		w.beginArray();
		}
	@Override
	public void write(Map<String, Object> hash) {
		try {
			w.beginObject();
			for(String col:super.columns) {
				w.name(col);
				final Object o = hash.get(col);
				if(o==null) {
					w.nullValue();
					}
				else if(o instanceof Boolean) {
					w.value(Boolean.class.cast(o));
					}
				else if(o instanceof Number) {
					w.value(Number.class.cast(o));
					}
				else {
					w.value(objectToString(o));
					}
				}
			w.endObject();
			} 
		catch (IOException e) {
			throw new RuntimeIOException(e);
			}
		}
	@Override
	public void flush() throws IOException {
		w.flush();
		}
	@Override
	public void close() throws IOException {
		w.endArray();
		w.close();
		out.close();
		}
	}
 
 static class DelimTextWriter extends AbstractTableWriter {
	final PrintWriter out;
	final char sep;
	DelimTextWriter(Path p,final List<String> columns,char sep) throws IOException {
		super(columns);
		this.out=IOUtils.openPathForPrintWriter(p);
		this.sep=sep;
		this.out.println(String.join( String.valueOf(sep), columns));
		}
	DelimTextWriter(Writer p,final List<String> columns,char sep) throws IOException {
		super(columns);
		this.out = (p instanceof PrintWriter?(PrintWriter)p: new PrintWriter(p));
		this.sep=sep;
		this.out.println(String.join( String.valueOf(sep), columns));
		}
	@Override
	public void write(Map<String, Object> hash) {
		for(int i=0;i< getColumns().size();i++) {
			if(i>0) this.out.print(this.sep);
			final String col = getColumns().get(i);
			final Object o = hash.get(col);
			if(o!=null) {
				this.out.print(objectToString(o));
				}
			}
		this.out.println();
		}
	@Override
	public void flush() throws IOException {
		out.flush();
		}
	@Override
	public void close() throws IOException {
		out.close();
		}
	}

 private class TSVWriter extends DelimTextWriter {
	 TSVWriter(Path p,final List<String> columns) throws IOException {
		super(p,columns,'\t'); 
	 	}
	 TSVWriter(Writer os,final List<String> columns) throws IOException {
		super(os,columns,'\t'); 
	 	}
 	}

 private class CSVWriter extends DelimTextWriter {
	 CSVWriter(Path p,final List<String> columns) throws IOException {
		super(p,columns,','); 
	 	}
 	}
 
 
public TableWriter open(final Path path,List<String> columns) throws IOException {
	String fname= path.getFileName().toString().toLowerCase();
	if(fname.endsWith(".json.gz") || fname.endsWith(".json")) {
		return new JsonTableWriter(path,columns);
		}
	if(fname.endsWith(".csv.gz") || fname.endsWith(".csv")) {
		return new CSVWriter(path,columns);
		}
	if(fname.endsWith(".tsv.gz") || fname.endsWith(".tsv")) {
		return new TSVWriter(path,columns);
		}
	return new TSVWriter(path,columns);
	}
public TableWriter open(final Writer os,List<String> columns) throws IOException {
	return new TSVWriter(os,columns);
	}


}
