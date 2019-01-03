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
package com.github.lindenb.jvarkit.tools.vcfserver;

import java.io.Closeable;
import java.io.File;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.AbstractHandler;
import org.eclipse.jetty.server.handler.HandlerList;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.misc.VcfToTable;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.JavascriptVariantFilter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;


/**

BEGIN_DOC

## Input

Input is a set of indexed VCF file (tabix or tribble) or a file containing the path to the VCFs.


## Screenshot

https://twitter.com/yokofakun/status/923870331659485184

![twitter](https://pbs.twimg.com/media/DNI-7GZX0AA41qF.jpg "Screenshot")

https://twitter.com/yokofakun/status/924307836968079361

![twitter](https://pbs.twimg.com/media/DNPNBdQWsAAF-5w.jpg "Screenshot")



## Example 

```
$ java -jar dist/vcfserver.jar input.vcf.gz

2017-10-27 23:53:04.140:INFO::main: Logging initialized @510ms
[INFO][VcfServer]Starting com.github.lindenb.jvarkit.tools.vcfserver.VcfServer on http://localhost:8080
2017-10-27 23:53:04.223:INFO:oejs.Server:main: jetty-9.3.7.v20160115
2017-10-27 23:53:04.336:INFO:oejs.ServerConnector:main: Started ServerConnector@9a8472{HTTP/1.1,[http/1.1]}{0.0.0.0:8080}
2017-10-27 23:53:04.337:INFO:oejs.Server:main: Started @717ms

```


END_DOC

**/

@Program(name="vcfserver",
description="Web Server displaying VCF file. A web interface for vcf2table",
keywords={"vcf","table","visualization","server","web"})
public class VcfServer extends Launcher{
private static final Logger LOG=Logger.build(VcfServer.class).make();
private static final int DEFAULT_LIMIT=100;
private static final String REGION_PARAM="rgn";
private static final String VCFIDX_PARAM="vcf";
private static final String LIMIT_PARAM="limit";
private static final String JAVASCRIPT_PARAM="js";
private static final String JEXL_PARAM="jexl";
private static final String SHOW_HEADER_PARAM="header";
private static final String HIDE_NOCALL_PARAM="nc";
private static final String HIDE_HOMREF_PARAM="hr";
private static final String HIDE_GENOTYPES_PARAM="gt";
private static final String TEXT_FORMAT_PARAM="txt";

@Parameter(names={"-p","--ped","--pedigree"},description="Optional Pedigree file:"+Pedigree.OPT_DESCRIPTION)
private File pedigreeFile=null;
@Parameter(names={"-P","--port","-port"},description="Server listening port")
private int port=8080;	
@Parameter(names={"-timeout","--timeout"},description="query timeout in seconds")
private long timeout_seconds =60;	
@Parameter(names={"-jexl","--jexl"},description="Use/Show JEXL filter instead of Javascript filter (which is not filesystem-safe).")
private boolean use_jexl = false;
@Parameter(names={"--url"},description=Launcher.USER_CUSTOM_INTERVAL_URL_DESC)
private String userCustomUrl=null;

/** used to escape the text output in pre tag */
private static class EscapeXmlOutputStream
	extends FilterOutputStream
	{
	final boolean escape;
	EscapeXmlOutputStream(final OutputStream delegate, boolean escape)
		{
		super(delegate);
		this.escape=escape;
		}
	private void _write(final String s) throws IOException
		{
		super.out.write(s.getBytes());
		}
	@Override
	public void write(int b) throws IOException {
		if(!this.escape)
			{
			super.write(b);
			}
		else
			{
			switch(b) {
				case '>': _write("&gt;");break; 
				case '<': _write("&lt;");break; 
				case '&': _write("&amp;");break; 
				case '\'': _write("&apos;");break; 
				case '\"': _write("&quot;");break;
				default:super.write(b); break;
				}
			}
		}
	@Override
	public void close() throws IOException {
		flush();
		//nothing do  not close
		}
	}


private class ViewVcfHandler extends AbstractHandler
	{
	private final List<File> vcfFiles;
	
	
	ViewVcfHandler(final List<File> vcfFiles)
		{
		this.vcfFiles = vcfFiles;
		}
	
	private abstract class DelegateHandler implements Closeable
		{
		private String title="";
		protected final HttpServletRequest request;
		protected final HttpServletResponse response;
		protected XMLStreamWriter writer = null;
		ViewVcfHandler getOwner() { return ViewVcfHandler.this;}
		DelegateHandler(final HttpServletRequest request,final HttpServletResponse response)
			{ 
			this.request = request ;
			this.response = response ;
			}
		
		void flush() {
			try
				{
				this.response.getOutputStream().flush();
				if(this.writer!=null) this.writer.flush();
				this.response.getOutputStream().flush();
				}
			catch(IOException|XMLStreamException err)
				{
				LOG.error(err);
				}
			}
		
		String getTitle() {
			return this.title;
			}
		
		void writeForm()throws XMLStreamException
			{
			this.writer.writeComment("BEGIN FORM");
			this.writer.writeStartElement("div");
			this.writer.writeStartElement("form");
			writeSelectVcf();
			
			//no call
			{
			this.writer.writeStartElement("span");
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id",HIDE_NOCALL_PARAM);
			this.writer.writeAttribute("type","checkbox");
			this.writer.writeAttribute("name",HIDE_NOCALL_PARAM);
			if("true".equals(request.getParameter(HIDE_NOCALL_PARAM))) 	this.writer.writeAttribute("checked","true");
			this.writer.writeAttribute("value","true");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for",HIDE_NOCALL_PARAM);
			this.writer.writeCharacters("Hide NO_CALL");
			this.writer.writeEndElement();//label
			this.writer.writeEndElement();//span
			}
			
			
			//no Hom_ref
			{
			this.writer.writeStartElement("span");
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id",HIDE_HOMREF_PARAM);
			this.writer.writeAttribute("type","checkbox");
			this.writer.writeAttribute("name",HIDE_HOMREF_PARAM);
			if("true".equals(request.getParameter(HIDE_HOMREF_PARAM))) 	this.writer.writeAttribute("checked","true");
			this.writer.writeAttribute("value","true");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for",HIDE_HOMREF_PARAM);
			this.writer.writeCharacters("Hide HOM_REF");
			this.writer.writeEndElement();//label
			this.writer.writeEndElement();//span
			}

			//no GT
			{
			this.writer.writeStartElement("span");
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id",HIDE_GENOTYPES_PARAM);
			this.writer.writeAttribute("type","checkbox");
			this.writer.writeAttribute("name",HIDE_GENOTYPES_PARAM);
			if("true".equals(request.getParameter(HIDE_GENOTYPES_PARAM))) 	this.writer.writeAttribute("checked","true");
			this.writer.writeAttribute("value","true");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for",HIDE_GENOTYPES_PARAM);
			this.writer.writeCharacters("Hide Genotypes");
			this.writer.writeEndElement();//label
			this.writer.writeEndElement();//span
			}

			//text output
			{
			this.writer.writeStartElement("span");
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id",TEXT_FORMAT_PARAM);
			this.writer.writeAttribute("type","checkbox");
			this.writer.writeAttribute("name",TEXT_FORMAT_PARAM);
			if("true".equals(request.getParameter(TEXT_FORMAT_PARAM))) 	this.writer.writeAttribute("checked","true");
			this.writer.writeAttribute("value","true");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for",TEXT_FORMAT_PARAM);
			this.writer.writeCharacters("Text Format");
			this.writer.writeEndElement();//label
			this.writer.writeEndElement();//span
			}
			
			//header
			final String show_header_str = request.getParameter(SHOW_HEADER_PARAM);
			this.writer.writeStartElement("span");
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id",SHOW_HEADER_PARAM);
			this.writer.writeAttribute("type","checkbox");
			this.writer.writeAttribute("name",SHOW_HEADER_PARAM);
			if("true".equals(show_header_str)) 	this.writer.writeAttribute("checked","true");
			this.writer.writeAttribute("value","true");
			
			
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for",SHOW_HEADER_PARAM);
			this.writer.writeCharacters("Show VCF Header");
			this.writer.writeEndElement();//label
			
			this.writer.writeEndElement();//span
			
			//region
			final String rgn_str = request.getParameter(REGION_PARAM);
			this.writer.writeStartElement("span");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for","selrgn");
			this.writer.writeCharacters("Interval");
			this.writer.writeEndElement();
			this.writer.writeCharacters(" : ");
			
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id","selrgn");
			this.writer.writeAttribute("type","text");
			this.writer.writeAttribute("name",REGION_PARAM);
			this.writer.writeAttribute("placeholder","chrom, chrom:pos+extend, chrom:start-end");
			this.writer.writeAttribute("value",StringUtil.isBlank(rgn_str)?"":rgn_str);
			this.writer.writeEndElement();//span
			//
			//limit
			final String limit_str = request.getParameter(LIMIT_PARAM);
			this.writer.writeStartElement("span");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for","limitctx");
			this.writer.writeCharacters("Limit number of variants");
			this.writer.writeEndElement();
			this.writer.writeCharacters(" : ");
			
			this.writer.writeStartElement("select");
			this.writer.writeAttribute("id","limitctx");
			this.writer.writeAttribute("name",LIMIT_PARAM);
			for(int L : new int[] {1,10,20,50,100,1000})
				{
				this.writer.writeStartElement("option");
				this.writer.writeAttribute("value",String.valueOf(L));
				if((!StringUtil.isBlank(limit_str) && limit_str.equals(String.valueOf(L))) || L==10)
					{
						this.writer.writeAttribute("selected","true");
					}
				this.writer.writeCharacters(String.valueOf(L));
				this.writer.writeEndElement();
				}
			this.writer.writeEndElement();
			
			this.writer.writeEndElement();//span
			//
			if(VcfServer.this.use_jexl)
				{
				//jexl
				final String jex = request.getParameter(JEXL_PARAM);
				this.writer.writeStartElement("span");
				this.writer.writeStartElement("label");
				this.writer.writeAttribute("for",JEXL_PARAM);
				this.writer.writeAttribute("title","JEXL expression");
				this.writer.writeCharacters("JEXL Expression");
				this.writer.writeEndElement();
				this.writer.writeCharacters(" : ");
				
				this.writer.writeEmptyElement("input");
				this.writer.writeAttribute("id",JEXL_PARAM);
				this.writer.writeAttribute("type","text");
				this.writer.writeAttribute("placeholder","see hhttps://software.broadinstitute.org/gatk/documentation/article.php?id=1255");
				this.writer.writeAttribute("name",JEXL_PARAM);
				this.writer.writeAttribute("value",StringUtil.isBlank(jex)?"":jex);
				this.writer.writeEndElement();//span
				}
			else
				{
				//javascript
				final String js = request.getParameter(JAVASCRIPT_PARAM);
				this.writer.writeStartElement("span");
				this.writer.writeStartElement("label");
				this.writer.writeAttribute("for",JAVASCRIPT_PARAM);
				this.writer.writeAttribute("title","Javascript Filtering Expression");
				this.writer.writeCharacters("Javascript Expression");
				this.writer.writeEndElement();
				this.writer.writeCharacters(" : ");
				
				this.writer.writeEmptyElement("input");
				this.writer.writeAttribute("id",JAVASCRIPT_PARAM);
				this.writer.writeAttribute("type","text");
				this.writer.writeAttribute("placeholder","see https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/filter/JavascriptVariantFilter.html");
				this.writer.writeAttribute("name",JAVASCRIPT_PARAM);
				this.writer.writeAttribute("value",StringUtil.isBlank(js)?"":js);
				this.writer.writeEndElement();//span
				}
			//
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("type","submit");
			
			this.writer.writeEndElement();//form
			this.writer.writeEndElement();//div
			this.writer.writeComment("END FORM");
			}
		
		void writeSelectVcf()throws XMLStreamException
			{
			if(getOwner().vcfFiles.size()==1) {
				this.writer.writeEmptyElement("input");
				this.writer.writeAttribute("type","hidden");
				this.writer.writeAttribute("name",VCFIDX_PARAM);
				this.writer.writeAttribute("value","0");
				}
			else
				{
				final File selVcf = getOwner().getVcfFile(this.request);
				this.writer.writeStartElement("span");
				this.writer.writeStartElement("label");
				this.writer.writeAttribute("for","selvcf");
				this.writer.writeCharacters("VCF");
				this.writer.writeEndElement();
				this.writer.writeCharacters(" : ");
				this.writer.writeStartElement("select");
				this.writer.writeAttribute("name",VCFIDX_PARAM);
				this.writer.writeAttribute("id","selvcf");
				for(int i=0;i< getOwner().vcfFiles.size();++i)
					{
					this.writer.writeStartElement("option");
					if((selVcf!=null && selVcf.equals(getOwner().vcfFiles.get(i))) ||i==0)
						{
						this.writer.writeAttribute("selected","1");
						}
					this.writer.writeAttribute("name",VCFIDX_PARAM);
					this.writer.writeAttribute("value",String.valueOf(i));
					this.writer.writeCharacters(getOwner().vcfFiles.get(i).getPath());
					this.writer.writeEndElement();
					}
				this.writer.writeEndElement();
				this.writer.writeEndElement();
				}
			}
		
		void writeHtmlTitle() throws XMLStreamException
			{
			this.writer.writeStartElement("title");
			this.writer.writeCharacters(getTitle());
			this.writer.writeEndElement();
			}
		
		void writeHtmlHead() throws XMLStreamException
			{
			writeHtmlTitle();
			writer.writeStartElement("style");
			writer.writeCharacters(VcfToTable.DEFAULT_CSS_STYLE);
			writer.writeEndElement();//style
			writer.writeStartElement("style");
			writer.writeCharacters(
					"label {color: #B4886B; font-weight: bold;}\n"+
					"input[type=text] {border: 1px dotted #999; border-radius: 0;}");

			writer.writeCharacters(".error {background-color:yellow;color:red;");
			writer.writeEndElement();//style
			}
		
		void writeHtmlFooter()  throws XMLStreamException
			{
			flush();
			writer.writeStartElement("div");
			writer.writeCharacters("Author: Pierre Lindenbaum. made with ");
			writer.writeStartElement("a");
			writer.writeAttribute("href", "https://github.com/lindenb/jvarkit");
			writer.writeAttribute("title", "https://github.com/lindenb/jvarkit");
			writer.writeAttribute("target", "_blank");
			writer.writeCharacters("jvarkit");
			writer.writeEndElement();//a
			writer.writeCharacters(".");
			writer.writeEndElement();//div
			this.writer.flush();
			}
		
		
		void writeHtmlBody() throws XMLStreamException
			{
		
			}
		
		void writeException(final Throwable err) throws XMLStreamException
			{
			LOG.error(err);
			final StringWriter sw=new StringWriter();
			final PrintWriter pw=new PrintWriter(sw);
			err.printStackTrace(pw);
			pw.close();
			this.writer.writeStartElement("pre");
			this.writer.writeAttribute("class","error");
			this.writer.writeCharacters(sw.toString());
			this.writer.writeEndElement();
			this.writer.flush();
			}
		void writeError(final String err) throws XMLStreamException
			{
			LOG.error(err);
			this.writer.writeStartElement("span");
			this.writer.writeAttribute("class","error");
			
			this.writer.writeCharacters(err);
			this.writer.writeEndElement();
			this.writer.flush();
			}
		
		void run() {
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			try
				{
				String encoding = this.request.getCharacterEncoding();
				if(StringUtil.isBlank(encoding)) encoding="UTF-8";
				this.request.setCharacterEncoding(encoding);
				this.response.setContentType("text/html");
				this.writer = xof.createXMLStreamWriter(this.response.getOutputStream(),encoding);
				this.writer.writeStartElement("html");
				this.writer.writeComment("BEGIN-VCF-SERVER");
				this.writer.writeStartElement("head");
				this.writer.writeEmptyElement("meta");
				this.writer.writeAttribute("charset", encoding);
				this.writer.writeEmptyElement("meta");
					this.writer.writeAttribute("name","Description");
					this.writer.writeAttribute("content","VCF Server");
				this.writer.writeEmptyElement("meta");
					this.writer.writeAttribute("name","keywords");
					this.writer.writeAttribute("content","VCF,variant,snp,server,bioinformatics");
				writeHtmlHead();
				this.writer.writeEndElement();//head
				this.writer.writeStartElement("body");
				this.writer.flush();
				writeHtmlBody();
				writeHtmlFooter();
				this.writer.writeEndElement();//BODY
				this.writer.writeComment("END-VCF-SERVER");
				this.writer.writeEndElement();
				this.writer.flush();
				flush();
				this.writer.close();
				}
			catch(final Exception err)
				{
				LOG.error(err);
				}
			finally
				{
				CloserUtil.close(this.writer);
				this.writer=null;
				}
			}
		@Override
		public void close() throws IOException {
			
			}
		}
		
	private class WelcomeHandler extends DelegateHandler
		{
		WelcomeHandler(
				final HttpServletRequest request,
				final HttpServletResponse response
				)
			{ 
			super(request,response);
			}
		@Override
		String getTitle() { return VcfServer.class.getSimpleName();}
		@Override
		void writeHtmlBody() throws XMLStreamException {
			writeForm();
			}
		
		}
	
	private class ShowVcfHandler extends DelegateHandler
		{
		ShowVcfHandler(
				final HttpServletRequest request,
				final HttpServletResponse response
				)
			{ 
			super(request,response);
			}
		@Override
		String getTitle() {
			final String rgn=request.getParameter(REGION_PARAM);
			return StringUtil.isBlank(rgn)?"No Region specified":rgn;
			}
		@Override
		void writeHtmlBody() throws XMLStreamException {
			final String rgn_str=request.getParameter(REGION_PARAM);
			final File vcfFile= getOwner().getVcfFile(this.request);
			writeForm();
			writer.writeEmptyElement("hr");
			flush();
			
			
			VCFFileReader reader=null;
			CloseableIterator<VariantContext> iter=null;
			try
				{
				reader = new VCFFileReader(vcfFile, true);
				final VCFHeader header = reader.getFileHeader();
				if(header==null)
					{
					writeError("Cannot get header of "+header);
					return;
					}
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				if(dict==null || dict.isEmpty())
					{
					writeError("Empty or null Dictionary in "+vcfFile);
					return;
					}
				
				final Interval interval;
				
				if(!StringUtil.isBlank(rgn_str)) 
					{
					final IntervalParser parser= new IntervalParser(dict);
					parser.setFixContigName(true);
					parser.setContigNameIsWholeContig(true);
					parser.setRaiseExceptionOnError(false);
					interval = parser.parse(rgn_str);	
					}
				else
					{
					interval = null;
					}
				final Predicate<VariantContext> variantPredicate;
				final String js_str= (VcfServer.this.use_jexl ?null:this.request.getParameter(JAVASCRIPT_PARAM));
				final String jexl_str= (VcfServer.this.use_jexl ?this.request.getParameter(JEXL_PARAM):null);

				if(!StringUtil.isBlank(jexl_str) && VcfServer.this.use_jexl)
					{
				
					try
						{
						final List<JexlVCMatchExp> exps= VariantContextUtils.initializeMatchExps(
								Collections.singletonList("CUSTOM_JEXL_FILTER"),
								Collections.singletonList(jexl_str)
								);
						variantPredicate = (V)-> VariantContextUtils.match(V,exps.get(0));
						}
					catch(final Exception err)
						{
						writeException(err);
						return;
						}
					}
				else if(!StringUtil.isBlank(js_str) && !VcfServer.this.use_jexl)
					{
					try
						{
						final StringReader scriptReader = new StringReader(js_str);
						final JavascriptVariantFilter jsFilter=new JavascriptVariantFilter(scriptReader, header);
						scriptReader.close();
						variantPredicate = (V)-> jsFilter.test(V);
						}
					catch(final Exception err)
						{
						writeException(err);
						return;
						}
					}
				else
					{
					variantPredicate =  (V)->true;
					}
				this.writer.writeComment("BEGIN-TABLE");
				this.writer.writeCharacters("");
				
				final boolean text_output= "true".equals(this.request.getParameter(TEXT_FORMAT_PARAM));

				if(text_output) {
					this.writer.writeStartElement("pre");
					this.writer.writeCharacters("");
					}
				this.flush();
				
				
				final VcfToTable.VcfToTableViewer vcfToTable=new VcfToTable.VcfToTableViewer();
				vcfToTable.setOutputFormat(text_output?
						VcfToTable.OutputFormat.text:
						VcfToTable.OutputFormat.html
						);
				final PrintStream newOut= new PrintStream(
						new EscapeXmlOutputStream(
						IOUtils.uncloseableOutputStream(this.response.getOutputStream()),
						text_output
						));
				vcfToTable.setOutputStream(newOut);
				vcfToTable.setHideHtmlHeader(true);//always
				vcfToTable.setPrintHeader("true".equals(this.request.getParameter(SHOW_HEADER_PARAM)));
				vcfToTable.setHideGenotypes("true".equals(this.request.getParameter(HIDE_GENOTYPES_PARAM)));
				vcfToTable.setHideHomRefGenotypes("true".equals(this.request.getParameter(HIDE_HOMREF_PARAM)));
				vcfToTable.setHideNoCallGenotypes("true".equals(this.request.getParameter(HIDE_NOCALL_PARAM)));
				vcfToTable.setUseANSIColors(!text_output);
				vcfToTable.setUserCustomUrl(VcfServer.this.userCustomUrl);
				
				vcfToTable.writeHeader(header);
				if(VcfServer.this.pedigreeFile!=null)
					{
					vcfToTable.setPedigreeFile(VcfServer.this.pedigreeFile);
					}
				if(interval==null)
					{
					iter = reader.iterator();
					}
				else
					{
					iter = reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
					}
				int limit=DEFAULT_LIMIT;
				final String limit_str = request.getParameter(LIMIT_PARAM);
				try {
					limit=StringUtil.isBlank(limit_str)?
						DEFAULT_LIMIT:
						Integer.parseInt(limit_str)
						;
					}
				catch(NumberFormatException err)
					{
					limit=DEFAULT_LIMIT;
					}
				final long start_millisec = System.currentTimeMillis();
				boolean timeout_flag = false;
				while(iter!=null && iter.hasNext() && limit>0)
					{
					final VariantContext ctx = iter.next();
					if(!variantPredicate.test(ctx)) continue;
					
					final long now_millisec = System.currentTimeMillis();
					if( now_millisec - start_millisec > VcfServer.this.timeout_seconds * 1000L)
						{
						timeout_flag=true;
						break;
						}	
					
					vcfToTable.add(ctx);
					--limit;
					}
				
				vcfToTable.close();
				
				newOut.flush();
				newOut.close();
				
				if(text_output)
					{
					this.writer.writeCharacters("");
					this.writer.writeEndElement();//pre
					}
				if(iter!=null && iter.hasNext())
					{
					this.writer.writeStartElement("p");
					this.writer.writeAttribute("class", "error");
					this.writer.writeCharacters("WARNING: there are more variants but limit was reached.");
					this.writer.writeEndElement();
					this.writer.flush();
					}
				
				if(timeout_flag)
					{
					this.writer.writeStartElement("p");
					this.writer.writeAttribute("class", "error");
					this.writer.writeCharacters("Time out reached!");
					this.writer.writeEndElement();
					this.writer.flush();
					}
				
				
				this.writer.writeComment("END-TABLE");
				this.flush();
				}
			catch(final Exception err)
				{
				super.writeException(err);
				}
			finally
				{
				CloserUtil.close(iter);
				CloserUtil.close(reader);
				}
			}
		}

	
	
	private File getVcfFile(final HttpServletRequest req) {
		
		if(this.vcfFiles.size()==1) return this.vcfFiles.get(0);
		
		final String vcfifx=req.getParameter(VCFIDX_PARAM);
		if(StringUtil.isBlank(vcfifx)) return null;
		final int idx;
		try { idx=Integer.parseInt(vcfifx);} catch(NumberFormatException err) {return null;}
		if(idx<0 || idx>=this.vcfFiles.size()) return null;
		return vcfFiles.get(idx);
		}
	
	
	
	public void handle(String arg0,
			final org.eclipse.jetty.server.Request jetty,
			final HttpServletRequest req,
			final javax.servlet.http.HttpServletResponse res
			) throws java.io.IOException ,javax.servlet.ServletException
		{
		DelegateHandler delegate=null;
		final File file = this.getVcfFile(req);
		
		if(file==null)
			{
			delegate = new WelcomeHandler(req,res);
			}
		else 
			{
			delegate = new ShowVcfHandler(req,res);
			}
		/*if(delegate==null)
			{
			delegate = new WelcomeHandler(req,res);
			}*/
		try {
			delegate.run();
			}
		catch(final Exception err)
			{
			LOG.error(err);
			}
		finally
			{
			CloserUtil.close(delegate);
			}
		}
	}
@Override
public int doWork(final List<String> args) {
	Server server = null;
	try
		{
		final List<File> vcfFiles = IOUtil.unrollFiles(args.stream().
			map(S->new File(S)).
			collect(Collectors.toList()),
			".vcf",".vcf.gz"
			);
		if(vcfFiles.isEmpty())
			{
			LOG.error("No VCF file defined");
			return -1;
			}
		vcfFiles.forEach(F->IOUtil.assertFileIsReadable(F));
		
		server = new Server(this.port);
		final HandlerList handlers = new HandlerList();
		handlers.addHandler(new ViewVcfHandler(vcfFiles));
		server.setHandler(handlers);
		LOG.info("Starting "+VcfServer.class.getName()+" on http://localhost:"+this.port);
		server.start();
		server.join();
		return 0;
		}
	catch(final Exception err)  {
		LOG.error(err);
		return -1;
		}
	finally
		{
		if(server!=null)
			{
			
			server.destroy();
			}
		}
	}

public static void main(final String args[]) {
	new VcfServer().instanceMainWithExit(args);
	}
}