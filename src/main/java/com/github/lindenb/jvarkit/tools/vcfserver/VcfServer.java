package com.github.lindenb.jvarkit.tools.vcfserver;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.eclipse.jetty.server.Handler;
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
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

@Program(name="vcfserver",
description="Web Server for vcf2table",
keywords={"vcf","table","visualization","server","web"})
public class VcfServer extends Launcher{
private static final Logger LOG=Logger.build(VcfServer.class).make();
private static final int DEFAULT_LIMIT=100;
private static final String REGION_PARAM="rgn";
private static final String VCFIDX_PARAM="vcf";
private static final String LIMIT_PARAM="limit";
private static final String JEXL_PARAM="jexl";
private static final String SHOW_HEADER_PARAM="header";

@Parameter(names={"-p","--ped","--pedigree"},description="Optional Pedigree file:"+Pedigree.OPT_DESCRIPTION+" If undefined, this tool will try to get the pedigree from the header.")
private File pedigreeFile=null;

@Parameter(names={"-P","--port","-port"},description="Server listening port")
private int port=8080;	


private class ViewVcfHandler extends AbstractHandler
	{
	private final List<File> vcfFiles;
	
	
	ViewVcfHandler(List<File> vcfFiles)
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
		DelegateHandler(HttpServletRequest request,HttpServletResponse response)
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
			this.writer.writeAttribute("value",StringUtil.isBlank(rgn_str)?"":rgn_str);
			this.writer.writeEndElement();//span
			//
			//limit
			final String limit_str = request.getParameter(LIMIT_PARAM);
			this.writer.writeStartElement("span");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for","limitctx");
			this.writer.writeCharacters("Limit");
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
			//jexl
			final String jexl = request.getParameter(JEXL_PARAM);
			this.writer.writeStartElement("span");
			this.writer.writeStartElement("label");
			this.writer.writeAttribute("for","jexl");
			this.writer.writeCharacters("JEXL Expression");
			this.writer.writeEndElement();
			this.writer.writeCharacters(" : ");
			
			this.writer.writeEmptyElement("input");
			this.writer.writeAttribute("id","jexl");
			this.writer.writeAttribute("type","text");
			this.writer.writeAttribute("name",JEXL_PARAM);
			this.writer.writeAttribute("value",StringUtil.isBlank(jexl)?"":jexl);
			this.writer.writeEndElement();//span
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
			writer.writeCharacters(".error {background-color:yellow;color:red;");
			writer.writeEndElement();//style

			}
		void writeHtmlBody() throws XMLStreamException
			{
		
			}
		
		void writeException(final Throwable err) throws XMLStreamException
			{
			LOG.error(err);
			StringWriter sw=new StringWriter();
			PrintWriter pw=new PrintWriter(sw);
			err.printStackTrace(pw);
			pw.close();
			this.writer.writeStartElement("pre");
			this.writer.writeAttribute("class","error");
			this.writer.writeCharacters(sw.toString());
			this.writer.writeEndElement();
			}
		void writeError(final String err) throws XMLStreamException
			{
			LOG.error(err);
			this.writer.writeStartElement("span");
			this.writer.writeAttribute("class","error");
			
			this.writer.writeCharacters(err);
			this.writer.writeEndElement();
			}
		
		void run() {
			VCFFileReader vcfReader=null;
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			try
				{
				String encoding = this.request.getCharacterEncoding();
				if(StringUtil.isBlank(encoding)) encoding="UTF-8";
				this.request.setCharacterEncoding(encoding);
				this.writer = xof.createXMLStreamWriter(this.response.getOutputStream(),encoding);
				vcfReader =new VCFFileReader(getVcfFile(this.request),true);
				this.writer.writeStartElement("html");
				this.writer.writeComment("BEGIN-VCF-SERVER");
				this.writer.writeStartElement("head");
				writeHtmlHead();
				this.writer.writeEndElement();//head
				this.writer.writeStartElement("body");
				this.writer.flush();
				writeHtmlBody();
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
				CloserUtil.close(vcfReader);
				CloserUtil.close(this.writer);
				this.writer=null;
				}
			}
		@Override
		public void close() throws IOException {
			
			}
		}
	
	private class DelegateErrorHandler extends DelegateHandler
		{
		private final Throwable error;
		private final String error_msg;
		DelegateErrorHandler(
				HttpServletRequest request,
				HttpServletResponse response,
				Throwable error
				)
			{ 
			super(request,response);
			this.error = error;
			this.error_msg = error.getMessage();
			}
		DelegateErrorHandler(
				HttpServletRequest request,
				HttpServletResponse response,
				String error
				)
			{ 
			super(request,response);
			this.error = null;
			this.error_msg = error;
			}
		@Override
		String getTitle() {
			return this.error_msg;
			}
		@Override
		void writeHtmlBody() throws XMLStreamException {
			this.writer.writeStartElement("h1");
			this.writer.writeCharacters(getTitle());
			this.writer.writeEndElement();
			
			if(this.error!=null)
				{
				writeException(this.error);
				}
			else
				{
				writeError(this.error_msg);
				}
			}
		}
	
	private class WelcomeHandler extends DelegateHandler
		{
		WelcomeHandler(
				HttpServletRequest request,
				HttpServletResponse response
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
				HttpServletRequest request,
				HttpServletResponse response
				)
			{ 
			super(request,response);
			}
		@Override
		String getTitle() {
			final String rgn=request.getParameter(REGION_PARAM);
			return StringUtil.isBlank(rgn)?"??":rgn;
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
				
				final IntervalParser parser= new IntervalParser(dict);
				parser.setFixContigName(true);
				parser.setRaiseExceptionOnError(false);
				final Interval interval = parser.parse(rgn_str);	
				if(interval==null)
					{
					writeError("Cannot parse interval "+rgn_str);
					return;
					}
				final String jex_str= this.request.getParameter(JEXL_PARAM);
				final Predicate<VariantContext> variantPredicate;
				
				if(!StringUtil.isBlank(jex_str))
					{
					try
						{
						final List<JexlVCMatchExp> exps= VariantContextUtils.initializeMatchExps(
								Collections.singletonList("CUSTOM"),
								Collections.singletonList(jex_str)
								);
						variantPredicate = (V)-> VariantContextUtils.match(V,exps.get(0));
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
				this.flush();
				
				
				
				final VcfToTable.VcfToTableViewer vcfToTable=new VcfToTable.VcfToTableViewer();
				vcfToTable.setOutputFormat(VcfToTable.OutputFormat.html);
				final PrintStream newOut= new PrintStream(IOUtils.uncloseableOutputStream(this.response.getOutputStream()));
				vcfToTable.setOutputStream(newOut);
				vcfToTable.setHideHtmlHeader(!("true".equals(this.request.getParameter(SHOW_HEADER_PARAM))));
				
				vcfToTable.writeHeader(header);
				if(VcfServer.this.pedigreeFile!=null)
					{
					vcfToTable.setPedigreeFile(VcfServer.this.pedigreeFile);
					}
				
				iter = reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
				
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
				
				while(iter!=null && iter.hasNext() && limit>0)
					{
					final VariantContext ctx = iter.next();
					if(!variantPredicate.test(ctx)) continue;
					vcfToTable.add(ctx);
					--limit;
					}
				if(iter!=null && iter.hasNext())
					{
					this.writer.writeStartElement("p");
					this.writer.writeCharacters("WARNING: there are more variants in "+rgn_str);
					this.writer.writeEndElement();
					}
				vcfToTable.close();
				
				newOut.flush();
				newOut.close();
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
		
		final String rgn_str = req.getParameter(REGION_PARAM);
		final File file = this.getVcfFile(req);
		
		if(StringUtil.isBlank(rgn_str) ||file==null)
			{
			delegate = new WelcomeHandler(req,res);
			}
		else 
			{
			delegate = new ShowVcfHandler(req,res);
			}
		if(delegate==null)
			{
			delegate = new WelcomeHandler(req,res);
			}
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
		final List<File> vcfFiles = args.stream().
			map(S->new File(S)).
			collect(Collectors.toList())
			;
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