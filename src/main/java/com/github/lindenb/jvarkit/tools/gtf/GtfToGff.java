/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gtf;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.gff.Gff3Constants;

/**

BEGIN_DOC

```
```

END_DOC
*/
@Program(
		name="gtf2gff",
		description="Convert GTF to gff",
		creationDate="20220703",
		modificationDate="20250328",
		keywords= {"gtf","gff","gff3"}
		)
public class GtfToGff
	extends Launcher {
	private static final Logger LOG = Logger.build(GtfToGff.class).make();
	
	
	 
	private class Record  implements Locatable {
		final List<String> tokens=new ArrayList<>(8);
		final Map<String,List<String>> attributes = new LinkedHashMap<>();
		
		Record(final String[] tokens) {
			for(int i=0;i< 8;i++) this.tokens.add(tokens[i]);
			
			final String s= tokens[8];
			int i=0;
		
			while(i< s.length()) {
				/* skip ws */
				while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
				final StringBuilder key=new StringBuilder();
				while(i< s.length()) {
					if(Character.isWhitespace(s.charAt(i))) {
						i++;
						break;
						}
					key.append(s.charAt(i));
					i++;
					}
				/* skip ws */
				while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
				if(i>=s.length()) throw new IllegalArgumentException("END of string after "+s.substring(0,i));
				boolean got_quote=false;
				if(s.charAt(i)=='\"') {
					got_quote=true;
					i++;
					}
				
				final StringBuilder value=new StringBuilder();
				while(i< s.length() ) {
					if(got_quote && s.charAt(i)=='\"') break;
					else if(!got_quote && s.charAt(i)==';') break;
					else if(!got_quote && Character.isWhitespace(s.charAt(i))) break;
					
					if(s.charAt(i)=='\\') {
						if(i+1>=s.length()) throw new IllegalArgumentException("unknown escape sequence after "+s.substring(0,i));
						i++;
						switch(s.charAt(i)) {
							case '"': value.append("\"");break;
							case '\'': value.append("\'");break;
							case '\\': value.append("\\");break;
							case 't': value.append("\t");break;
							case 'n': value.append("\n");break;
							default: throw new IllegalArgumentException("unknown escape sequence after "+s.substring(0,i));
							}
						}
					else
						{
						value.append(s.charAt(i));
						}
					i++;
					}
				if(got_quote) {
					if(i>=s.length()) throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
					if(s.charAt(i)!='\"')  throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
					i++;
					}
				final String keyStr= key.toString();
				List<String> L = this.attributes.get(keyStr);
				if(L==null) {
					L=new ArrayList<>(2);
					this.attributes.put(keyStr,L);
					}
				L.add(value.toString());
				
				/* skip ws */
				while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
				if(i<s.length()) {
					if(s.charAt(i)!=';')  throw new IllegalArgumentException("expected semicolon after "+s.substring(0,i));
					i++;
					}
				/* skip ws */
				while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
				}
			}
		
		@Override
		public String getContig()
			{
			return tokens.get(0);
			}
		@Override
		public int getStart()
			{
			return Integer.parseInt(tokens.get(3));
			}
		@Override
		public int getEnd()
			{
			return Integer.parseInt(tokens.get(4));
			}
		private String escape(String value) {
			return value;
			}
		void print(PrintWriter pw) {
			for(int i=0;i< tokens.size();i++) {
				if(i>0) pw.print(Gff3Constants.FIELD_DELIMITER);
				pw.print(tokens.get(i));
				}
			for(String key:this.attributes.keySet()) {
				for(String value:this.attributes.get(key)) {
					pw.print(key);
					pw.print(Gff3Constants.KEY_VALUE_SEPARATOR);
					pw.print(escape(value));
					}
				}
			pw.print(Gff3Constants.END_OF_LINE_CHARACTER);
			}
		}
	
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;


	@Override
	public int doWork(final List<String> args)  {
		try {

			final String input = oneFileOrNull(args);
			try(InputStream in1 = input==null?System.in:Files.newInputStream(Paths.get(input))) {
				final InputStream in = IOUtils.uncompress(in1);
				final BufferedReader br = new BufferedReader(new InputStreamReader(in, "UTF-8"));
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					String line;
					while((line=br.readLine())!=null) {
						if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
						final String[] tokens = CharSplitter.TAB.split(line);
						if(tokens.length!=9) throw new IllegalArgumentException("expected 9 tokens in "+String.join("<\\t>", tokens));
						final Record record=  new Record(tokens);
						record.print(out);
						}
					out.flush();
					}
				br.close();
				in.close();
				}
		
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args)
		{
		new GtfToGff().instanceMainWithExit(args);
		}
	}
