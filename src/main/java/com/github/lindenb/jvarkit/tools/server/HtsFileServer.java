/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.server;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.eclipse.jetty.http.HttpStatus;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## input

input is a set of indexed Vcf file or a file with the suffix `.list` containing the path to the vcfs.
 

END_DOC
 
 */

@Program(name="htsfileserver",
	description="Jetty Based http server serving Vcf and Bam files.",
	creationDate="20200405",
	modificationDate="20200405",
	keywords={"vcf","server"},
	generate_doc=false
	)
public  class HtsFileServer extends Launcher {
	private static final Logger LOG = Logger.build(HtsFileServer.class).make();
	
	@Parameter(names="--port",description="server port.")
	private int serverPort = 8080;
	@Parameter(names= {"--gtf"},description="Optional GTF file. Will be used to retrieve an interval by gene name")
	private Path gtfFile = null;
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;

	
	private abstract class AbstractInput {
		private final Path path;
		private final String md5;
		AbstractInput(final Path path) throws IOException{
			this.path = path;
			IOUtil.assertFileIsReadable(path);
			this.md5 = StringUtils.md5(path.getFileName().toString()).substring(0,10);
			}
		public String getMd5()
			{
			return md5;
			}
		public Path getPath()
			{
			return path;
			}
		abstract String getOutputName();

		abstract void dump(final Locatable loc,final OutputStream os) throws IOException;
		}
	
	private class VcfInput extends AbstractInput {
		VcfInput(final Path path) throws IOException{
			super(path);
			// try open with index
			try(VCFFileReader ignore= new VCFFileReader(path, true)) {
				//nothing
				}
			}
		@Override
		String getOutputName()
			{
			String s= getPath().getFileName().toString();
			if(!s.endsWith(".gz")) s+=".gz";
			return s;
			}
		
		@Override
		void dump(Locatable loc, OutputStream os) throws IOException
			{
			VariantContextWriterBuilder vcb = new VariantContextWriterBuilder();
			vcb.setOutputStream(os);
			vcb.setCreateMD5(false);
			try(VariantContextWriter w=vcb.build()) {
				try(VCFFileReader r = new VCFFileReader(getPath(), true)) {
					final VCFHeader header = r.getFileHeader();
					final SAMSequenceDictionary dict = header.getSequenceDictionary();
					w.writeHeader(header);
					
					CloseableIterator<VariantContext> iter;
					if(loc==null) {
						iter = r.iterator();
						}
					else
						{
						iter = r.query(loc);
						}
					while(iter.hasNext()) {
						w.add(iter.next());
						}
					}
				
				}
			}
		}
	
	private class BamInput extends AbstractInput {
		final CRAMReferenceSource ref;
		BamInput(final Path path,final CRAMReferenceSource ref) throws IOException{
			super(path);
			this.ref = ref;
			}
		@Override
		String getOutputName()
			{
			String s= getPath().getFileName().toString();
			if(!s.endsWith(FileExtensions.BAM)) s+=FileExtensions.BAM;
			return s;
			}
		@Override
		void dump(Locatable loc, OutputStream os) throws IOException
			{
			SAMFileWriterFactory vcb = new SAMFileWriterFactory();
		
			
			try(SamReader r = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).
						referenceSource(this.ref).
						open(getPath())) {
				final SAMFileHeader header = r.getFileHeader();
				try(SAMFileWriter w=vcb.makeBAMWriter(header, true, os)) {
					CloseableIterator<SAMRecord> iter;
					if(loc==null) {
						iter = r.iterator();
						}
					else
						{
						iter = r.query(loc.getContig(), loc.getStart(), loc.getEnd(), false);
						}
					while(iter.hasNext()) {
						w.addAlignment(iter.next());
						}
					}
				
				}
			}
		}

		
		
	private final Map<String,AbstractInput> htsMap = new HashMap<>();
	private SAMSequenceDictionary dictionary = null;
	
	@SuppressWarnings("serial")
	private  class HtsFileServerServlet extends HttpServlet {
		@Override
		protected void doGet(final HttpServletRequest req,final  HttpServletResponse resp) throws ServletException, IOException {
			doPost(req, resp);
			}
		@Override
		protected void doPost(final HttpServletRequest request, final HttpServletResponse response) throws ServletException, IOException {
			final String action = request.getParameter("action");
			if("dump".equals(action))
				{
				dumpData(request,response);
				}
			else
				{
				printPage(request,response);
				}
			}
		}
		
	
	private SimpleInterval parseInterval(final String s) {
		if(StringUtils.isBlank(s)) return null;
		
		/* search in GTF file */
		if(this.gtfFile!=null && !s.contains(":")) {
			final String geneName  = s.trim();
			TabixReader tbr =null;
			
			}
		
		
		return IntervalParserFactory.newInstance(this.dictionary).
			enableSinglePoint().
			make().
			apply(s.trim()).
			orElse(null);
	}
	

			
	
	
	/** dumpData */
	private void dumpData(final HttpServletRequest request,final HttpServletResponse response)	throws IOException, ServletException
		{
		final String fileids[]  = request.getParameterValues("htsfileid");
		final List<AbstractInput> selected = new ArrayList<>();
		if(fileids!=null) {
			for(final String fid: fileids) {
				final AbstractInput input = this.htsMap.get(fid);
				if(input!=null) selected.add(input);
			}
		}
		
		if(selected.isEmpty()) selected.addAll(this.htsMap.values());
		
		Locatable loc = null;
		String intervalstr = request.getParameter("interval");
		if(!StringUtils.isBlank(intervalstr)) {
			 loc = IntervalParserFactory.newInstance().enableWholeContig().dictionary(this.dictionary).make().apply(intervalstr).orElse(null);
		}
		if(!StringUtils.isBlank(intervalstr) && this.gtfFile!=null && loc==null && !intervalstr.contains(':')) {
			TabixReader tbr;
			try {
			String geneName = intervalstr;
			final ContigNameConverter cvt = ContigNameConverter.fromOneDictionary(this.dictionary);
			tbr= new TabixReader(this.gtfFile.toString());
			final GTFCodec codec = new GTFCodec();
			String line;
			while((line=tbr.readLine())!=null) {
				if(StringUtils.isBlank(line) ||  line.startsWith("#")) continue;
				final String tokens[]= CharSplitter.TAB.split(line);
				if(tokens.length<9 ) continue;
				if(!(tokens[2].equals("gene") || tokens[2].equals("transcript"))) continue;
				if(StringUtils.indexOfIgnoreCase(tokens[8],geneName)==-1) continue;
				final GTFLine gtfLine = codec.decode(line);
				if(gtfLine==null) continue;
				
				if(tokens[2].equals("gene") ) {
					if(!(geneName.equals(gtfLine.getAttribute("gene_id")) || geneName.equals(gtfLine.getAttribute("gene_name")))) continue;
				}
				else if(tokens[2].equals("transcript") ) {
					if(!(geneName.equals(gtfLine.getAttribute("transcript_id")))) continue;
				}
				
				final String ctg = cvt.apply(gtfLine.getContig());
				if(StringUtils.isBlank(ctg)) continue;
				tbr.close();
				tbr = null;
				loc = new SimpleInterval(ctg,gtfLine.getStart(),gtfLine.getEnd());
				break;
				}
			}
	catch(final Throwable err) {
		
		}
	finally
		{
		if(tbr!=null) tbr.close();
		}
	}
		
		if(!StringUtils.isBlank(intervalstr) && this.gtfFile!=null && loc==null ) {
			loc = new SimpleInterval("undefined_interval",1,1);
		}
		
		
		String prefix= StringUtils.now()+".";
		if(loc!=null) {
			prefix += new SimpleInterval(loc).toNiceString()+".";
			}
		
		try(PrintStream out = new PrintStream(response.getOutputStream())) {
			final String charset = StringUtils.ifBlank(request.getCharacterEncoding(), "UTF-8");
			response.setCharacterEncoding(charset);

			if(selected.size()==1) {
				final String fname = prefix+selected.get(0).getOutputName();
				response.addHeader("Content-Disposition","form-data; name=\""+ fname +"\"; filename=\""+ fname +"\"");
				response.setContentType("text/html; charset="+charset.toLowerCase());

				selected.get(0).dump(loc, out);
				}
			else
				{
				final String fname = prefix+HtsFileServer.class.getSimpleName();
				 response.addHeader("Content-Disposition","form-data; name=\""+ fname +"\"; filename=\""+ fname +"\"");

				response.setContentType("text/html; charset="+charset.toLowerCase());
				final ZipOutputStream zout = new ZipOutputStream(out, Charset.forName(charset));
				for(AbstractInput input: selected) {
					final ZipEntry zipEntry= new ZipEntry(prefix+input.getOutputName());
					zout.putNextEntry(zipEntry);
					input.dump(loc, zout);
					zout.closeEntry();
					if(out.checkError()) break;
					}
				
				zout.finish();
				}
			out.flush();
			}
		
		
		}

	/** print HTML page */
	private void printPage(final HttpServletRequest request,final HttpServletResponse response)	throws IOException, ServletException
		{
		 
		 final String charset = StringUtils.ifBlank(request.getCharacterEncoding(), "UTF-8");
		 response.setContentType("text/html; charset="+charset.toLowerCase());
		 response.setCharacterEncoding(charset);
		 PrintWriter pw = response.getWriter();
		
		
		 
		 final String title = "HtsFileServer";
				
		 
		 try {
			final XMLStreamWriter w=XMLOutputFactory.newFactory().createXMLStreamWriter(pw);
			w.writeStartElement("html");
			
			w.writeStartElement("head");
			
			w.writeEmptyElement("meta");
			w.writeAttribute("charset", charset);
			
			w.writeStartElement("title");
			w.writeCharacters(title);
			w.writeEndElement();//title
			
			w.writeStartElement("style");
			w.writeCharacters(
					"body {background-color:#f0f0f0;color:#070707;font-size:18px;}"
					+ "h1 {text-align:center;color:#070707;}"
					+ ".span1 {border:1px dotted blue;}"
					+ ".lbl {font-weight: bold;}"
					+ ".bampath {font-family:monospace;font-size:12px; font-style: normal;color:gray;}"
					+ ".headerform {background-color:lavender;text-align:center;font-size:14px;}"
					+ ".comment {background-color:khaki;text-align:center;font-size:14px;}"
					+ ".highlight {background-color:#DB7093;}"
					+ ".parents {font-size:75%;color:gray;}"
					+ ".children {font-size:75%;color:gray;}"
					+ ".allsamples {font-size:125%;}"
					+ ".message {color:red;}"
					+ ".affected {background-color:#e6cccc;}"
					+ ".gtf {background-color:moccasin;text-align:center;}"
					+ ".known {background-color:wheat;text-align:center;}"
					);
			w.writeEndElement();//title
			
			w.writeStartElement("script");
			w.writeEndElement();//script
			
			w.writeEndElement();//head
			
			w.writeStartElement("body");
			
			
			
			w.writeStartElement("h1");
			
			w.writeCharacters(title);
			
			
			
			w.writeEmptyElement("a");
			w.writeAttribute("name", "top");

			
			w.writeComment("BEGIN FORM");
			w.writeStartElement("div");
			w.writeAttribute("class", "headerform");
			w.writeStartElement("form");
			w.writeAttribute("method", "GET");
			w.writeAttribute("action", "/page");
			
			w.writeEmptyElement("input");
			w.writeAttribute("name", "action");
			w.writeAttribute("type", "hidden");
			w.writeAttribute("value","dump");

				
			w.writeEmptyElement("input");
			w.writeAttribute("name", "interval");
			w.writeAttribute("id", "interval");
			w.writeAttribute("type", "text");
			w.writeAttribute("value","");
			
			w.writeStartElement("button");
			w.writeAttribute("class", "btn");
			w.writeAttribute("name", "go");
			w.writeCharacters("GO");
			w.writeEndElement();//button
				
			w.writeEndElement();//form
			w.writeEndElement();//div
			
			w.writeComment("END FORM");
			
			
			
			w.writeStartElement("ul");
			w.writeAttribute("class", "grid-container");
			
			for( AbstractInput htsFile: this.htsMap.values()) {
				w.writeStartElement("li");
			
				w.writeEmptyElement("input");
				w.writeAttribute("type", "checkbox");
				w.writeAttribute("name","htsfileid");
				w.writeAttribute("id",htsFile.getMd5());
				w.writeAttribute("value",htsFile.getMd5());
				
				w.writeStartElement("label");
				w.writeAttribute("for",htsFile.getMd5());
				w.writeCharacters(htsFile.getPath().toString());
				w.writeEndElement();//label
				w.writeEndElement(); //lli
			}
			
			w.writeEndElement();//ul
			
			w.writeEmptyElement("hr");
			w.writeStartElement("div");
			w.writeCharacters("Author: Pierre Lindenbaum. ");
			w.writeCharacters(JVarkitVersion.getInstance().getLabel());
			w.writeEndElement();
			
			w.writeEndElement();//body
			w.writeEndElement();//html
			w.flush();
			w.close();
		 	}
		 catch(XMLStreamException err) { LOG.warn(err);throw new IOException(err);}
		 finally { pw.close(); }
		}
		
	
	
	
	@Override
	public int doWork(final List<String> args) {
	
	
		try {
			final CRAMReferenceSource ref;
			if(this.faidxRef!=null) {
				ref = new ReferenceSource(this.faidxRef);
				this.dictionary  = SequenceDictionaryUtils.extractRequired(this.faidxRef);
				}
			else
				{
				ref = null;
				}
		
			for(final Path path: IOUtils.unrollPaths(args)) {
				final String fn = path.getFileName().toString();
				final AbstractInput input;
				if(fn.endsWith(FileExtensions.BAM) || fn.endsWith(FileExtensions.CRAM)) {
					input = new BamInput(path,ref);
					}
				else if(fn.endsWith(FileExtensions.VCF) || fn.endsWith(FileExtensions.COMPRESSED_VCF)) {
					input = new VcfInput(path);
					}
				else
					{
					LOG.error("unsupported format "+path);
					return -1;
					}
				if(this.htsMap.containsKey(input.getMd5())) {
					LOG.error("duplicate key "+input.getMd5()+" "+input.getPath());
					return -1;
					}
				this.htsMap.put(input.getMd5(), input);
				}
			if(this.htsMap.isEmpty()) {
				LOG.error("No VCF/BAM defined.");
				return -1;
				}
			
			
			final Server server = new Server(this.serverPort);
			
			final ServletContextHandler context = new ServletContextHandler();
	        context.addServlet(new ServletHolder(new HtsFileServerServlet()),"/*");
	        context.setContextPath("/");
	        context.setResourceBase(".");
	        server.setHandler(context);
			
		    
		    
		    LOG.info("Starting server "+getProgramName()+" on port "+this.serverPort);
		    server.start();
		    LOG.info("Server started. Press Ctrl-C to stop. Check your proxy settings ."
		    		+ " Open a web browser at http://localhost:"+this.serverPort+"/htsserver .");
		    server.join();
		    return 0;
			}
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}	


public static void main(final String[] args) throws Exception{
    new HtsFileServer().instanceMainWithExit(args);
	}

}
