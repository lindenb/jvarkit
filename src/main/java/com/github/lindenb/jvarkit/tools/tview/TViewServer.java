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
package com.github.lindenb.jvarkit.tools.tview;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Collections;
import java.util.List;
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
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.tview.TView.Formatout;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.JavascriptSamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;


/**

BEGIN_DOC

## Input

Input is a set of indexed BAM files  or a file containing the path to the BAMs.

## Screenshot

https://twitter.com/yokofakun/status/925395272225746944

![twitter](https://pbs.twimg.com/media/DNepyjOW4AA-_rz.jpg "Screenshot")


## Example 

```
$ java -jar dist/tviewserver.jar -R src/test/resources/toy.fa src/test/resources/toy.bam
2017-10-31 16:31:15.281:INFO::main: Logging initialized @1262ms
[INFO][TViewServer]Starting com.github.lindenb.jvarkit.tools.tview.TViewServer on http://localhost:8080
2017-10-31 16:31:15.523:INFO:oejs.Server:main: jetty-9.3.7.v20160115
2017-10-31 16:31:15.690:INFO:oejs.ServerConnector:main: Started ServerConnector@14a2189{HTTP/1.1,[http/1.1]}{0.0.0.0:8080}
2017-10-31 16:31:15.691:INFO:oejs.Server:main: Started @1675ms

```


END_DOC

**/

@Program(name="tviewserver",
	description="Web Server displaying SAM/BAM file. A web interface for jvarkit:tview",
    keywords={"sam","bam","table","visualization","server","web"}
		)
public class TViewServer extends Launcher{
private static final Logger LOG=Logger.build(TViewServer.class).make();
private static final String REGION_PARAM="rgn";
private static final String SAMIDX_PARAM="samfile";
private static final String JAVASCRIPT_PARAM="js";
private static final String SHOWCLIP="clip";
private static final String SHOWINSERT="insert";
private static final String HIDEBASES="bases";
private static final String SHOWNAME="name";
private static final String SHOWALLBAMS="showall";

@Parameter(names={"-R","--reference"},description=Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION)
private File optionalReferenceFile=null;
@Parameter(names={"-P","--port","-port"},description="Server listening port")
private int port=8080;	
@Parameter(names={"-nojs","--no-javascript"},description="Disable Javascript (which is not filesystem-safe).")
private boolean disable_javascript = false;
@Parameter(names={"-m","--max"},description="Max interval Length")
private int max_interval_length = 2000;
@Parameter(names={"--url"},description=Launcher.USER_CUSTOM_INTERVAL_URL_DESC)
private String userCustomUrl=null;
@Parameter(names={"--shutdown-after"},description="Stop the server after 'x' seconds.")
private long shutdownAferSeconds=-1L;


private class SamViewHandler extends AbstractHandler
	{
	private final List<File> samFiles;
	
	
	SamViewHandler(final List<File> samFiles)
		{
		this.samFiles = samFiles;
		}
	
	
	
	private abstract class DelegateHandler implements Closeable
		{
		private String title="";
		protected final HttpServletRequest request;
		protected final HttpServletResponse response;
		protected XMLStreamWriter writer = null;
		SamViewHandler getOwner() { return SamViewHandler.this;}
		DelegateHandler(final HttpServletRequest request,final HttpServletResponse response)
			{ 
			this.request = request ;
			this.response = response ;
			}
		
		protected boolean showAllBamsInOneWindow() {
			return  TViewServer.showAllBamsInOneWindow(this.request);
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
			writeSelectSam();
			
			
			final String flags[]=new String[] {
					SHOWCLIP,"Show Clip",
					HIDEBASES,"Hide Bases",
					SHOWINSERT,"Show Inserts",
					SHOWNAME,"Show Read Names",
					SHOWALLBAMS,"Show All Bams in one Window"
					};


			for(int i=0;i< flags.length;i+=2)
				{
				final String id=flags[i];
				final String lbl=flags[i+1];
				this.writer.writeStartElement("span");
				this.writer.writeEmptyElement("input");
				this.writer.writeAttribute("id",id);
				this.writer.writeAttribute("type","checkbox");
				this.writer.writeAttribute("name",id);
				if("true".equals(request.getParameter(id))) 	this.writer.writeAttribute("checked","true");
				this.writer.writeAttribute("value","true");
				this.writer.writeStartElement("label");
				this.writer.writeAttribute("for",id);
				this.writer.writeCharacters(lbl);
				this.writer.writeCharacters(".");
				this.writer.writeEndElement();//label
				this.writer.writeCharacters(" ");
				this.writer.writeEndElement();//span
				}
		
			
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
			
			//
			if(!TViewServer.this.disable_javascript)
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
				this.writer.writeAttribute("placeholder","see https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/filter/SamRecordFilter.html");
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
		
		
		
		void writeSelectSam()throws XMLStreamException
			{
			if(getOwner().samFiles.size()==1) {
				this.writer.writeEmptyElement("input");
				this.writer.writeAttribute("type","hidden");
				this.writer.writeAttribute("name",SAMIDX_PARAM);
				this.writer.writeAttribute("value","0");
				}
			else
				{
				final File selVcf = getOwner().getSamFile(this.request);
				this.writer.writeStartElement("span");
				this.writer.writeStartElement("label");
				this.writer.writeAttribute("for",SAMIDX_PARAM);
				this.writer.writeCharacters("BAM");
				this.writer.writeEndElement();
				this.writer.writeCharacters(" : ");
				this.writer.writeStartElement("select");
				this.writer.writeAttribute("name",SAMIDX_PARAM);
				this.writer.writeAttribute("id",SAMIDX_PARAM);
				for(int i=0;i< getOwner().samFiles.size();++i)
					{
					this.writer.writeStartElement("option");
					if((selVcf!=null && selVcf.equals(getOwner().samFiles.get(i))) ||i==0)
						{
						this.writer.writeAttribute("selected","1");
						}
					this.writer.writeAttribute("name",SAMIDX_PARAM);
					this.writer.writeAttribute("value",String.valueOf(i));
					this.writer.writeCharacters(getOwner().samFiles.get(i).getPath());
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
			writer.writeCharacters("body {background-color:#F8F8F8;}");
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
				this.writer.writeComment("BEGIN-BAM-SERVER");
				this.writer.writeStartElement("head");
				this.writer.writeEmptyElement("meta");
				this.writer.writeAttribute("charset", encoding);
				this.writer.writeEmptyElement("meta");
					this.writer.writeAttribute("name","Description");
					this.writer.writeAttribute("content","SAM Server");
				this.writer.writeEmptyElement("meta");
					this.writer.writeAttribute("name","keywords");
					this.writer.writeAttribute("content","SAM,BAM,read,server,bioinformatics");
				writeHtmlHead();
				this.writer.writeEndElement();//head
				this.writer.writeStartElement("body");
				this.writer.flush();
				writeHtmlBody();
				writeHtmlFooter();
				this.writer.writeEndElement();//BODY
				this.writer.writeComment("END-BAM-SERVER");
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
		String getTitle() { return TViewServer.class.getSimpleName();}
		@Override
		void writeHtmlBody() throws XMLStreamException {
			writeForm();
			}
		
		}
	
	private class ShowBamHandler extends DelegateHandler
		{
		ShowBamHandler(
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
			int sam_file_index=0;
			writeForm();
			this.writer.writeEmptyElement("hr");
			
			flush();
			do {
				
				final File samFile;
				
				if(this.showAllBamsInOneWindow())
					{
					if(sam_file_index<0 || sam_file_index>=getOwner().samFiles.size())
						{
						writeError("Illegal state sam_file_index");
						return;
						}
					samFile = getOwner().samFiles.get(sam_file_index);
					}
				else
					{
					samFile  = getOwner().getSamFile(this.request);
					}
				
				/* paranoid */
				if(samFile==null)
					{
					writeError("Illegal state sam File is null");
					return;
					}
				
				this.writer.writeStartElement("p");
				this.writer.writeCharacters(samFile.getName());
				this.writer.writeEndElement();
				
				
				
				TView tview = new TView();
				try
					{
					final SAMSequenceDictionary dict= SAMSequenceDictionaryExtractor.extractDictionary(samFile);
					if(dict==null) {
						writeError("no dict in this bam file :"+samFile);
						return;
						}
					final Interval interval;
					
					if(!StringUtil.isBlank(rgn_str)) 
						{
						final IntervalParser parser= new IntervalParser(dict);
						parser.setFixContigName(true);
						parser.setContigNameIsWholeContig(true);
						parser.setRaiseExceptionOnError(false);
						interval = TViewServer.this.trimInterval(parser.parse(rgn_str));	
						}
					else
						{
						final SAMSequenceRecord rec = dict.getSequence(0);
						interval = TViewServer.this.trimInterval(new Interval(rec.getSequenceName(),1,Math.min(100, rec.getSequenceLength())));
						}
					
					tview.setInterval(interval);
					if(optionalReferenceFile!=null) tview.setReferenceFile(optionalReferenceFile);
					tview.setFormatOut(Formatout.html);
					tview.setShowClip("true".equals(this.request.getParameter(SHOWCLIP)));
					tview.setShowReadName("true".equals(this.request.getParameter(SHOWNAME)));
					tview.setShowInsertions("true".equals(this.request.getParameter(SHOWINSERT)));
					tview.setHideBases("true".equals(this.request.getParameter(HIDEBASES)));
					tview.setBamFiles(Collections.singletonList(SamInputResource.of(samFile)));
	
					if(!TViewServer.this.disable_javascript)
						{
						final String js_expr = this.request.getParameter(JAVASCRIPT_PARAM);
						if(!StringUtil.isBlank(js_expr))
							{
							final JavascriptSamRecordFilter filter;
							try 
								{
								StringReader strReader = new StringReader(js_expr);
								filter = new JavascriptSamRecordFilter(
										js_expr,
										SamReaderFactory.makeDefault().getFileHeader(samFile)
										);
								strReader.close();
								}
							catch(Exception err)
								{
								writeException(err);
								return;
								}
							tview.setSamRecordFilter(filter);
							}
						}
					
					if(tview.initialize()!=0)
						{
						writeError("cannot initialize tview");
						return ;
						}
					
					/* Hyperlink to IGV */
					if(!StringUtil.isBlank(TViewServer.this.userCustomUrl)) {
						final String gotostr=Launcher.createUrlFromInterval(
								TViewServer.this.userCustomUrl,
								interval
								);
						if(!StringUtil.isBlank(gotostr)) {
							this.writer.writeStartElement("div");
							this.writer.writeStartElement("a");
							this.writer.writeAttribute("title","URL");
							this.writer.writeAttribute("rel","nofollow");
							this.writer.writeAttribute("href", gotostr );
							this.writer.writeCharacters("[URL]");
							this.writer.writeEndElement();//a
							this.writer.writeEndElement();//div
							this.writer.writeCharacters("");
							}
						}
					
					this.writer.writeStartElement("pre");
					this.writer.writeCharacters("");
					this.writer.flush();
					
					final PrintStream out  = new PrintStream(IOUtils.uncloseableOutputStream(this.response.getOutputStream()));
					tview.paint(out);
					out.flush();
					tview.close();
					tview=null;
					this.writer.flush();
					this.writer.writeCharacters("");
					this.writer.writeEndElement();//pre
					this.writer.writeEmptyElement("hr");
					this.writer.writeCharacters("");
					this.flush();
					}
				catch(final Exception err)
					{
					super.writeException(err);
					}
				finally
					{
					CloserUtil.close(tview);
					}
				
				
				if(!this.showAllBamsInOneWindow()) {
					break;
					}
				sam_file_index++;
			} while(sam_file_index< getOwner().samFiles.size());
			
			}
		}

	
	
	private File getSamFile(final HttpServletRequest req) {
		
		if(this.samFiles.size()==1) return this.samFiles.get(0);
		if(TViewServer.showAllBamsInOneWindow(req)) return null;
		final String samidx=req.getParameter(SAMIDX_PARAM);
		if(StringUtil.isBlank(samidx)) return null;
		final int idx;
		try { idx=Integer.parseInt(samidx);} catch(NumberFormatException err) {return null;}
		if(idx<0 || idx>=this.samFiles.size()) return null;
		return samFiles.get(idx);
		}
	
	
	
	public void handle(String arg0,
			final org.eclipse.jetty.server.Request jetty,
			final HttpServletRequest req,
			final javax.servlet.http.HttpServletResponse res
			) throws java.io.IOException ,javax.servlet.ServletException
		{
		DelegateHandler delegate=null;
		final File samFile = this.getSamFile(req);
		
		if(!TViewServer.showAllBamsInOneWindow(req) && samFile==null)
			{
			delegate = new WelcomeHandler(req,res);
			}
		else 
			{
			delegate = new ShowBamHandler(req,res);
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

private Interval trimInterval(final Interval interval) {
	if(interval.length()<=this.max_interval_length) return interval;
	final Interval interval2 = new Interval(
			interval.getContig(),
			interval.getStart(),
			interval.getStart()+this.max_interval_length-1
			);
	
	LOG.debug("trim interval "+interval+" to "+interval2);
	return interval2;
	}

private static boolean showAllBamsInOneWindow(final HttpServletRequest request) {
	return "true".equals(request.getParameter(SHOWALLBAMS));
	}

@Override
public int doWork(final List<String> args) {
	Server server = null;
	try
		{
		final List<File> samFiles = IOUtil.unrollFiles(args.stream().
			map(S->new File(S)).
			collect(Collectors.toList()),
			".bam"
			);
		if(samFiles.isEmpty())
			{
			LOG.error("No BAM file defined");
			return -1;
			}
		
		
		
		if(this.optionalReferenceFile==null)
			{
			LOG.warn("No reference file defined");
			}
		else
			{
			final SAMSequenceDictionary dict =  SAMSequenceDictionaryExtractor.extractDictionary(this.optionalReferenceFile);
			if(dict==null || dict.isEmpty())
				{
				LOG.error("Empty/No dict in "+this.optionalReferenceFile);
				return -1;
				}
			samFiles.forEach(F->{
				final SAMSequenceDictionary samdict =  SAMSequenceDictionaryExtractor.extractDictionary(F);
				if(!SequenceUtil.areSequenceDictionariesEqual(samdict, dict)) {
					throw new JvarkitException.DictionariesAreNotTheSame(samdict, dict);
					}
				});
			}
		samFiles.forEach(F->{
			IOUtil.assertFileIsReadable(F);
			});
		
		server = new Server(this.port);
		final HandlerList handlers = new HandlerList();
		handlers.addHandler(new SamViewHandler(samFiles));
		server.setHandler(handlers);
		LOG.info("Starting "+TViewServer.class.getName()+" on http://localhost:"+this.port);
		server.start();
		if(this.shutdownAferSeconds>0)
			{
			final Server theServer = server;
			new java.util.Timer().schedule( 
			        new java.util.TimerTask() {
			            @Override
			            public void run() {
			                LOG.info("automatic shutdown after "+shutdownAferSeconds);
			                try {
			                	theServer.stop();
			                	}
			                catch(final Throwable err2) {
			                	LOG.error(err2);
			                	}
			            }
			        }, 
			        1000 * this.shutdownAferSeconds 
					);
			}
		
		server.join();
		return 0;
		}
	catch(final Throwable err)  {
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
	new TViewServer().instanceMainWithExit(args);
	}
}
