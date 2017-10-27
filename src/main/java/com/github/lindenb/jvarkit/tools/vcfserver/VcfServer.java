package com.github.lindenb.jvarkit.tools.vcfserver;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.misc.VcfToTable;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class VcfServer extends Launcher{
private static final Logger LOG=Logger.build(VcfServer.class).make();
private static final String REGION_PARAM="rgn";
private static final String VCFIDX_PARAM="vcf";

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
		protected XMLStreamWriter writer;
		ViewVcfHandler getOwner() { return ViewVcfHandler.this;}
		DelegateHandler(HttpServletRequest request,HttpServletResponse response)
			{ 
			this.request = request ;
			this.response = response ;
			}
		
		String getTitle() {
			return this.title;
			}
		
		void writeForm()throws XMLStreamException
			{
			this.writer.writeStartElement("div");
			this.writer.writeStartElement("form");
			writeSelectVcf();
			
			//region
			String rgn_str = request.getParameter(REGION_PARAM);
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
			
			
			this.writer.writeEndElement();//form
			this.writer.writeEndElement();//div
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
			}
		void writeHtmlBody() throws XMLStreamException
			{
		
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
				this.writer.writeStartElement("head");
				writeHtmlHead();
				this.writer.writeEndElement();//head
				this.writer.writeStartElement("body");
				writeHtmlBody();
				this.writer.writeEndElement();
				this.writer.writeEndElement();
				this.writer.flush();
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
			this.writer.writeStartElement("pre");
			if(this.error!=null)
				{
				StringWriter sw=new StringWriter();
				PrintWriter pw=new PrintWriter(sw);
				this.error.printStackTrace(pw);
				pw.close();
				this.writer.writeCharacters(sw.toString());
				}
			else
				{
				this.writer.writeCharacters(this.error_msg);
				}
			this.writer.writeEndElement();
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
			
			
			
			VCFFileReader reader=null;
			CloseableIterator<VariantContext> iter=null;
			try
				{
				reader = new VCFFileReader(vcfFile, true);
				final VCFHeader header = reader.getFileHeader();
				if(header==null)
					{
					writer.writeCharacters("Cannot get header of "+header);
					return;
					}
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				if(dict==null || dict.isEmpty())
					{
					writer.writeCharacters("Empty or null Dictionary in "+vcfFile);
					return;
					}
				
				final IntervalParser parser= new IntervalParser(dict);
				parser.setFixContigName(true);
				parser.setRaiseExceptionOnError(false);
				final Interval interval = parser.parse(rgn_str);	
				if(interval==null)
					{
					writer.writeCharacters("Cannot parse interval "+rgn_str);
					return;
					}
				this.writer.flush();
				
				
				
				final VcfToTable.VcfToTableViewer vcfToTable=new VcfToTable.VcfToTableViewer();
				vcfToTable.setOutputFormat(VcfToTable.OutputFormat.html);
				vcfToTable.setOutputStream(new PrintStream(IOUtils.uncloseableOutputStream(this.response.getOutputStream())));
				vcfToTable.setHideHtmlHeader(true);
				vcfToTable.writeHeader(header);
				
				iter = reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
				
				while(iter!=null && iter.hasNext())
					{
					final VariantContext ctx = iter.next();
					vcfToTable.add(ctx);
					}
				
				vcfToTable.close();
				}
			catch(final Exception err)
				{
				writer.writeCharacters("Error "+err.getMessage());
				LOG.error(err);
				return;
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
			LOG.info("START delegate class "+delegate.getClass());
			delegate.run();
			LOG.info("END");
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