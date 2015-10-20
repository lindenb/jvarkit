/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.net.URLDecoder;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

public class Gtf2Xml extends AbstractGtf2Xml{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Gtf2Xml.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractGtf2Xml.AbstractGtf2XmlCommand
	 	{	
	private abstract class GffCodec
		{
		Map<String,Long> seqdict=new LinkedHashMap<>();
		Set<String> att_keys=new HashSet<>();
		Set<String> sources=new HashSet<>();
		Set<String> types=new HashSet<>();
		protected Pattern tab=Pattern.compile("[\t]");
		protected Pattern semicolon=Pattern.compile("[;]");
		void write(XMLStreamWriter w,String line) throws XMLStreamException,IOException
			{
			String tokens[] = this.tab.split(line);
			if(tokens.length<8) throw new IOException("Expected at least 8 columns in "+line+" got "+tokens.length);
			w.writeStartElement("feature");
			
			w.writeAttribute("chrom", tokens[0]);
			w.writeAttribute("start", tokens[3]);
			w.writeAttribute("end", tokens[4]);
			
			Long contifLength  = seqdict.get(tokens[0]);
			if(contifLength==null) contifLength=0L;
			seqdict.put(tokens[0],Math.max(Long.parseLong(tokens[4]), contifLength));
			
			if(!tokens[5].equals("."))  w.writeAttribute("score", tokens[5]);
			if(!tokens[6].equals("."))  w.writeAttribute("strand", tokens[6]);
			if(!tokens[1].equals("."))
				{
				w.writeAttribute("source", tokens[1]);
				sources.add(tokens[1]);
				}
			if(!tokens[2].equals("."))
				{
				w.writeAttribute("type", tokens[2]);
				types.add(tokens[2]);
				}
			if(!tokens[7].equals(".")) w.writeAttribute("phase", tokens[7]);
			if(!tokens[8].equals("."))
				{
				w.writeStartElement("attributes");
				writeAttributes(w,tokens[8]);
				w.writeEndElement();
				}
			w.writeEndElement();
			w.writeCharacters("\n");
			}
		abstract void writeAttributes(XMLStreamWriter w,String attrString)  throws XMLStreamException,IOException;
			
		
		}
	
	private class DefaultCodec extends GffCodec
		{
		@Override
		void writeAttributes(XMLStreamWriter w,String attrString) throws XMLStreamException,IOException
			{
			StreamTokenizer st=new StreamTokenizer(new StringReader(attrString));
			st.wordChars('_', '_');
			String key=null;
			while(st.nextToken() != StreamTokenizer.TT_EOF)
				{
				String s=null;
				switch(st.ttype)
					{
					case StreamTokenizer.TT_NUMBER: s=String.valueOf(st.nval);break;
					case '"': case '\'' : case StreamTokenizer.TT_WORD: s=st.sval;break;
					case ';':break;
					default:break;
					}
				if(s==null) continue;
				if(key==null)
					{
					key=s;
					this.att_keys.add(key);
					}
				else 
					{
					w.writeStartElement(key);
					w.writeCharacters(s);
					w.writeEndElement();
					key=null;
					}
				}
			
			}
		}
	
	private class Gtf3Codec extends GffCodec
		{
		private String unescape(String s) throws IOException
			{
			return URLDecoder.decode(s, "UTF-8");
			}
		@Override
		void writeAttributes(XMLStreamWriter w,String attrString) throws XMLStreamException,IOException
			{
			String atts[]  = semicolon.split(attrString);
			for(String att: atts)
				{
				if(att.isEmpty()) continue;
				int eq=att.indexOf("=");
				if(eq<=0) throw new IOException("bad att "+att+" in "+attrString);
				String key = att.substring(0,eq);
				String value = att.substring(eq+1);
				w.writeStartElement(key);
				w.writeCharacters(unescape(value) );
				w.writeEndElement();
				this.att_keys.add(key);
				}
			}
		}
	
	
	
	@Override
		protected Collection<Throwable> call(String inputName) throws Exception
			{
			BufferedReader r=null;
			XMLStreamWriter w=null;
			FileWriter fw=null;
			try {
				if(inputName==null)
					{
					r = new BufferedReader(new InputStreamReader(stdin()));
					}
				else 
					{
					r = IOUtils.openURIForBufferedReading(inputName);
					}
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				if(this.outputFile==null)
					{
					w = xof.createXMLStreamWriter(stdout(), "UTF-8");
					}
				else
					{
					w = xof.createXMLStreamWriter((fw=new FileWriter(getOutputFile())));
					}
				GffCodec codec = new DefaultCodec();
				String headerLine=null;
				w.writeStartDocument("UTF-8","1.0");
				w.writeStartElement("gtf");
				while((headerLine=r.readLine())!=null)
					{
					if(!headerLine.startsWith("#")) break;
					int ws = headerLine.indexOf(' ');
					if(ws==-1) continue; //??
					
					if(headerLine.startsWith("##gff-version "))
						{
						String version=headerLine.substring(ws+1).trim();
						LOG.info("version "+version);
						if(version.equals("3"))
							{
							codec = new Gtf3Codec();
							}
						w.writeAttribute("gff-version", version);
						}
					else if(headerLine.startsWith("#!"))
						{
						w.writeAttribute(
								headerLine.substring(2, ws),
								headerLine.substring(ws+1).trim());
						}
					else
						{
						LOG.warn("ignoring "+headerLine);
						}
					}
				w.writeCharacters("\n");
				
				
			
				String line=null;
				for(;;)
					{
					line=(headerLine!=null?headerLine:r.readLine());
					headerLine=null;
					if(line==null) break;
					if(line.isEmpty()) continue;
					codec.write(w,line);
					}
				
				
				w.writeStartElement("attributes");
				for(String k : codec.att_keys)
					{
					w.writeStartElement("attribute");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				
				w.writeStartElement("types");
				for(String k : codec.types)
					{
					w.writeStartElement("type");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				
				
				w.writeStartElement("sources");
				for(String k : codec.sources)
					{
					w.writeStartElement("source");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
	
				w.writeStartElement("dict");
				for(String k : codec.seqdict.keySet())
					{
					w.writeEmptyElement("chrom");
					w.writeAttribute("name",k);
					w.writeAttribute("length", String.valueOf(codec.seqdict.get(k)));
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
	
				
				w.writeEndElement();
				w.writeEndDocument();
				w.flush();
				return RETURN_OK;
				}
			catch (Exception e) {
				LOG.error(e);
				return wrapException(e);
				}
			finally
				{
				CloserUtil.close(r);
				CloserUtil.close(fw);
				}
			}
	 	}
	
	public static void main(String[] args) throws IOException
		{
		new Gtf2Xml().instanceMain(args);
		}
}
