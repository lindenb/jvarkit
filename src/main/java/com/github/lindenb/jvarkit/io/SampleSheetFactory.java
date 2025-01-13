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

*/
package com.github.lindenb.jvarkit.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.github.lindenb.jvarkit.lang.CharSplitter;

public class SampleSheetFactory {
	private Function<String,List<String>> _splitter = null;
	
	public SampleSheetFactory() {
		splitter(CharSplitter.TAB);
	}
	
	public SampleSheetFactory setSplitter(Function<String, List<String>> splitter) {
		this._splitter = splitter;
		return this;
		}
	public SampleSheetFactory splitter(CharSplitter cs) {
		return setSplitter(S->cs.splitAsStringList(S));
		}
	
	public SampleSheet of(final Path path) throws IOException {
		try(BufferedReader br=IOUtils.openPathForBufferedReading(path)) {
			String line=br.readLine();
			if(line==null) throw new IOException("cannot read first line of "+path);
			final SampleSheetImpl ss = new SampleSheetImpl();
			ss.source = path.toString();
			ss.header =  new FileHeader(line,this._splitter);
			ss.addAll(ss.header.readAll(br));
			return ss;
			}
		catch(IOException err) {
			throw err;
			}
		catch(Throwable err2) {
			throw new IOException(err2);
			}
		}
	
	
	@SuppressWarnings("serial")
	private static class SampleSheetImpl extends ArrayList<FileHeader.RowMap> implements SampleSheet {
		FileHeader header;
		Object source;
		@Override
		public FileHeader getHeader() {
			return header;
			}
		@Override
		public String getSource() {
			return source.toString();
			}
		@Override
		public String toString() {
			StringBuilder sb=new StringBuilder();
			this.writeMarkdown(sb);
			return sb.toString();
			}
		}
}
