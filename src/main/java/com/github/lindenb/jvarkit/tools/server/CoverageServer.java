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
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## Example

```
java -jar dist/coverageserver.jar -r 'RF01:100-200' -R src/test/resources/rotirus_rf.fa src/test/resources/S*.bam

```
 
END_DOC
 
 */

@Program(name="coverageserver",
	description="Jetty Based http server serving VCF, BAM files",
	generate_doc=false,
	creationDate="20200212",
	modificationDate="20200212"
	)
public  class CoverageServer extends Launcher {
	private static final Logger LOG = Logger.build(CoverageServer.class).make();
	
	@Parameter(names="--port",description="server port")
	private int serverPort = 8080;
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;
	@Parameter(names= {"-r","--region","--vcf","--bed","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class)
	private IntervalListProvider intervalsource = IntervalListProvider.empty();

	private SAMSequenceDictionary dictionary;
	private final List<SimpleInterval> intervals = new Vector<>();
	private final List<Path> bamPaths = new Vector<>();
	
	
	private final int image_width= 800;
	private final int image_height= 400;

	
	
	@SuppressWarnings("serial")
	private  class CoverageServlet extends HttpServlet {
		@Override
		protected void doGet(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
			doPost(req, resp);
			}
		@Override
		protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
			  String pathInfo = request.getPathInfo();
			  if(pathInfo==null) pathInfo="/";
				LOG.info(pathInfo);
			if(pathInfo.equals("/getimage"))
				{
				printImage(request,response);
				}
			else
				{
				printPage(request,response);
				}
			}
		}
		
	
	private SimpleInterval parseInterval(final String s) {
		if(StringUtils.isBlank(s)) return null;
		return IntervalParserFactory.newInstance(this.dictionary).
			enableSinglePoint().
			make().
			apply(s).
			orElse(null);
	}

	
	private void printImage( HttpServletRequest request, HttpServletResponse response)	throws IOException, ServletException
	{
		int bam_id;
		try {
			bam_id = Integer.parseInt(StringUtils.ifBlank(request.getParameter("id"),"-1"));
		} catch(Exception err) {
			bam_id=-1;
		}
		final SimpleInterval region = parseInterval(request.getParameter("interval"));
		if(region==null || bam_id<0 || bam_id>=this.bamPaths.size()) {
			response.reset();
			response.sendError(HttpStatus.BAD_REQUEST_400,"id:"+bam_id);
			response.flushBuffer();
			return;
		}
		
		final Path bam = this.bamPaths.get(bam_id);
		final SamReaderFactory srf = SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).referenceSequence(this.faidxRef);
		try(SamReader sr=srf.open(bam)) {
			final SAMFileHeader header = sr.getFileHeader();
			if(region.getLengthOnReference()>10_000_000)  {
				response.reset();
				response.sendError(HttpStatus.BAD_REQUEST_400,"contig:"+region);
				response.flushBuffer();
				return;
			}
			LOG.info("drawing "+bam+" "+region);
			 final int coverage[]=new int[region.getLengthOnReference()];
			 Arrays.fill(coverage, 0);
			 try(CloseableIterator<SAMRecord> iter=sr.query(region.getContig(), region.getStart(), region.getEnd(),false)) {
				 while(iter.hasNext()) {
					 final SAMRecord rec=iter.next();
					 if(rec.getReadUnmappedFlag()) continue;
					 if(rec.getDuplicateReadFlag()) continue;
					 if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					 if(rec.isSecondaryOrSupplementary()) continue;
					 final Cigar cigar = rec.getCigar();
					 if(cigar==null || cigar.isEmpty()) continue;
					 int ref=rec.getAlignmentStart();
					 for(final CigarElement ce:cigar) {
						 final CigarOperator op=ce.getOperator();
						 if(op.consumesReferenceBases()) {
							 if(op.consumesReadBases()) {
								 for(int x=0;x< ce.getLength();++x) {
									 int pos=ref+x;
									 if(pos< region.getStart()) continue;
									 if(pos> region.getEnd()) break;
									 coverage[pos-region.getStart()]++;
								 }
							 }
							 ref+=ce.getLength();
						 }
					 }
				 }
			 }
			 if(coverage.length>image_width) {
				 final int copy[]=Arrays.copyOf(coverage, coverage.length);
				 int len = Math.max(1,coverage.length/1000);
				 
				 for(int i=0;i< coverage.length;i++) {
					 int j=Math.max(0, i-len);
					 double sum=0;
					 int count=0;
					 while(j< i+len && j< copy.length) {
						 sum +=copy[j];
						 j++;
						 count++;
					 }
					 coverage[i]=(int)(sum/count);
				 }
			 }
			 
			 double max_cov= Math.max(10, IntStream.of(coverage).max().orElse(1));
			 final double pixelperbase = image_width/(double)coverage.length;
			 final BufferedImage img = new BufferedImage(image_width, image_height, BufferedImage.TYPE_INT_RGB);
			 Graphics2D g=img.createGraphics();
			 g.setColor(Color.WHITE);
			 g.fillRect(0, 0, image_width+1, image_height+1);
			 
			 final Composite oldComposite = g.getComposite();
			 g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)Math.min(1.0, pixelperbase)));
			 
			 for(int x=0;x< coverage.length;++x) {
				 final double height = image_height*(coverage[x]/max_cov);
				 if(coverage[x]<10) g.setColor(Color.DARK_GRAY);
				 else if(coverage[x]<20) g.setColor(Color.BLUE);
				 else g.setColor(Color.GREEN);
				 
				 g.fill(new Rectangle2D.Double(
						 x*pixelperbase,
						 image_height-height,
						 pixelperbase,
						 height));
			 }
			 g.setComposite(oldComposite);
			 
			 g.drawString("max-cov:"+(int)max_cov+" sample:"+header.getReadGroups().
					 stream().
					 map(R->R.getSample()).filter(S->!StringUtils.isBlank(S)).
					 collect(Collectors.toSet()).
					 stream().
					 collect(Collectors.joining(" "))+" "+
					 region.toNiceString()
					 , 10, 10);

			 g.setColor(Color.BLUE);
			 for(int i=1;i<10;i++) {
				 double cov=Math.ceil(max_cov/10.0*i);
				 if(i>1 && cov==Math.ceil(max_cov/10.0*(i-1))) continue;
				 double y = image_height - image_height/10.0*i;
				 g.drawLine(0, (int)y, 5, (int)y);
				 g.drawString(String.valueOf((int)cov),7,(int)y);
			 }
			 
			 g.setColor(Color.GRAY);
			 g.drawRect(0, 0, img.getWidth(),  img.getHeight());

			 response.setContentType("image/png");
			 ImageIO.write(img, "PNG", response.getOutputStream());
			 response.flushBuffer();
		}
		
		
		
	}
		
	private void printPage( HttpServletRequest request, HttpServletResponse response)	throws IOException, ServletException
		{
		 String message="";
		 
		 final String charset = StringUtils.ifBlank(request.getCharacterEncoding(), "UTF-8");
		 response.setContentType("text/html; charset="+charset.toLowerCase());
		 response.setCharacterEncoding(charset);
		 PrintWriter pw = response.getWriter();
		 SimpleInterval interval = null;
		 
		 if(!StringUtils.isBlank(request.getParameter("custominterval"))) {
			 interval = parseInterval(request.getParameter("custominterval"));
			 if(interval==null) {
				 message+=" Cannot parse user interval '"+request.getParameter("custominterval")+"'.";
			 }
		 }
		 
		 if(interval==null && !StringUtils.isBlank(request.getParameter("interval"))) {
			interval = parseInterval(request.getParameter("interval"));
			if(interval==null) {
				 message+=" Cannot parse interval. using default: "+interval+".";
				}
		 	}
		 
		 if(interval==null && 	!this.intervals.isEmpty()) {
			 interval = this.intervals.get(0);
		 }
		 
		 if(interval==null) {
			 final SAMSequenceRecord ssr=this.dictionary.getSequence(0);
			 interval= new SimpleInterval(ssr.getSequenceName(), 1,Math.min(ssr.getSequenceLength(),100));
		 	}
		
		 if(!StringUtils.isBlank(request.getParameter("move"))) {
			 final String value = request.getParameter("move");
			 double factor;
			 switch(value.length()) {
			 	 case 0: factor=0;break;
				 case 1: factor=0.1;break;
				 case 2: factor=0.475;break;
				 default: factor = 0.95;break;
				 }
			 int length0 = interval.getLengthOnReference();
			 int length1 = (int)(length0*factor);
			 int shift= (value.startsWith(">")?1:-1) * length1;
			 int start = interval.getStart() + shift;
			 int end = interval.getEnd() + shift;
			 if(start<1) start=1;
			 final SAMSequenceRecord ssr=this.dictionary.getSequence(interval.getContig());
			 if(ssr!=null && end>=ssr.getSequenceLength()) {
				end= ssr.getSequenceLength();
			 	}
			 if(ssr!=null && start>=ssr.getSequenceLength()) {
				 	start= ssr.getSequenceLength();
				 	}
			interval = new SimpleInterval(interval.getContig(), start, end);
		 	}
		 
		 for(int side=0;side<2;++side) {
			 final String param = "zoom"+(side==0?"in":"out");
			 String value = request.getParameter(param);
			 if(StringUtils.isBlank(value)) continue;
			 double factor = Double.parseDouble(value);
			 if(side==0) factor=1.0/factor;
			 int length0 = interval.getLengthOnReference();
			 int length1 = (int)(length0*factor);
			 if(length1< 1) length1=1;
			 if(length1> 10_000_000) length1=1;
			 int mid = interval.getStart()+length0/2;
			 int start =mid-length1/2;
			 if(start<1) start=1;
			 int end =mid+length1/2;
			 final SAMSequenceRecord ssr=this.dictionary.getSequence(interval.getContig());
			 if(ssr!=null && end>=ssr.getSequenceLength()) {
				end= ssr.getSequenceLength();
			 	}
			 interval = new SimpleInterval(interval.getContig(), start, end);
			 break;
			 }
		 
		 
		 
		 final String title = interval.toNiceString()+" ("+StringUtils.niceInt(interval.getLengthOnReference())+" bp.)";
				
		 
		 try {
			final XMLStreamWriter w=XMLOutputFactory.newFactory().createXMLStreamWriter(pw);
			w.writeStartElement("html");
			
			w.writeStartElement("head");
			
			w.writeStartElement("title");
			w.writeCharacters(title);
			w.writeEndElement();//title
			
			w.writeStartElement("style");
			w.writeCharacters("");
			w.writeEndElement();//title
			
			w.writeStartElement("script");
			{
				w.writeCharacters(""); 
				w.flush();
				/* messy with < and > characters + xmlstream */
				pw.write(
				"function loadImage(idx) {"+
				"if(idx>="+this.bamPaths.size()+") return;"+
				"var img = document.getElementById(\"bamid\"+idx);"+
				"img.addEventListener('load',(event) => {img.width="+image_width+";img.height="+image_height+";loadImage(idx+1);});"+
				"img.setAttribute(\"src\",\"/getimage?id=\"+idx+\"&interval="+ StringUtils.escapeHttp(interval.toString()) +"\");"+
				"img.setAttribute(\"alt\",\"bam idx\"+idx);"+
				"}"+
				"function init() {"+
				"var span=document.getElementById(\"spaninterval\");"+
				"span.addEventListener('click',(evt)=>{document.getElementById(\"custominterval\").value = span.textContent; });"+
				"var sel=document.getElementById(\"selinterval\");"+
				"sel.addEventListener('change',(evt)=>{document.getElementById(\"custominterval\").value = evt.target.value; });"+
				"loadImage(0);"+
				"}"+
				"window.addEventListener('load', (event) => {init();});"
				);
				pw.flush();
			}
			w.writeEndElement();//script
			
			w.writeEndElement();//head
			
			w.writeStartElement("body");
			
			w.writeStartElement("h1");
			w.writeCharacters(title);
			w.writeEndElement();//h1
			
			if(!StringUtils.isBlank(message)) {
				w.writeStartElement("h2");
				w.writeCharacters(message);
				w.writeEndElement();//h1
			}
			
			w.writeStartElement("div");
			w.writeStartElement("form");
			w.writeAttribute("method", "GET");
			w.writeAttribute("action", "/page");
			
			w.writeEmptyElement("input");
			w.writeAttribute("name", "interval");
			w.writeAttribute("value", interval.toString());
			w.writeAttribute("type", "hidden");

			
			if(!this.intervals.isEmpty()) {
				w.writeStartElement("select");
				
				w.writeEmptyElement("option");
				
				for(final SimpleInterval r:this.intervals) {
					w.writeStartElement("option");
					w.writeAttribute("value",r.toString());
					if(r.equals(interval)) w.writeAttribute("selected", "true");
					w.writeCharacters(r.toNiceString());
					w.writeEndElement();
				}
				
				w.writeEndElement();//select
				
				//
				
				/* zoom in */
				w.writeStartElement("label");
				w.writeCharacters("move");
				w.writeEndElement();
				for(final String mv:new String[]{"<<<","<<","<",">",">>",">>>"}) {
					w.writeEmptyElement("input");
					w.writeAttribute("type", "submit");
					w.writeAttribute("name", "move");
					w.writeAttribute("value", mv);
				}
				
				/* zoom in */
				w.writeStartElement("label");
				w.writeCharacters("zoom in");
				w.writeEndElement();
				for(final String zoom:new String[]{"1.5","3","10"}) {
					w.writeEmptyElement("input");
					w.writeAttribute("type", "submit");
					w.writeAttribute("name", "zoomin");
					w.writeAttribute("value", zoom);
				}
				
				w.writeStartElement("span");
				w.writeAttribute("id", "spaninterval");
				w.writeCharacters(interval.toNiceString());
				w.writeEndElement();//span
				
				w.writeEmptyElement("input");
				w.writeAttribute("name", "custominterval");
				w.writeAttribute("id", "custominterval");
				w.writeAttribute("type", "text");
				w.writeAttribute("value","");
				
				w.writeStartElement("button");
				w.writeAttribute("name", "parseinterval");
				w.writeCharacters("GO");
				w.writeEndElement();//button
				
				
				/* zoom in */
				w.writeStartElement("label");
				w.writeCharacters("zoom out");
				w.writeEndElement();
				for(final String zoom:new String[]{"1.5","3","10","100"}) {
					w.writeEmptyElement("input");
					w.writeAttribute("type", "submit");
					w.writeAttribute("name", "zoomout");
					w.writeAttribute("value", zoom);
					}
				}
			w.writeEndElement();//form
			w.writeEndElement();//div
			
			w.writeStartElement("div");
			for(int i=0;i< this.bamPaths.size();i++) {
				final Path bam = this.bamPaths.get(i);
				if(!Files.exists(bam)) continue;
				w.writeStartElement("div");
				
				w.writeStartElement("h3");
				w.writeCharacters(IOUtils.getFilenameWithoutCommonSuffixes(bam));
				w.writeEndElement();
				w.writeStartElement("div");
				w.writeCharacters(bam.toString());
				w.writeEndElement();

				
				w.writeEmptyElement("img");
				w.writeAttribute("id","bamid"+i);
				w.writeAttribute("src", "https://upload.wikimedia.org/wikipedia/commons/d/de/Ajax-loader.gif");
				w.writeAttribute("width","32");
				w.writeAttribute("height","32");

				
				w.writeEndElement();
				}
			w.writeEndElement();
			
			w.writeEmptyElement("hr");
			w.writeStartElement("div");
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
			this.bamPaths.addAll(IOUtils.unrollPaths(args));
			if(this.bamPaths.isEmpty()) {
				LOG.error("No BAM defined.");
				return -1;
				}
			this.dictionary = SequenceDictionaryUtils.extractRequired(this.faidxRef);
			
			this.intervals.addAll(this.intervalsource.
					dictionary(this.dictionary).
					stream().
					map(T->new SimpleInterval(T)).
					collect(Collectors.toList()));
			final Server server = new Server(this.serverPort);
			
			final ServletContextHandler context = new ServletContextHandler();
	        context.addServlet(new ServletHolder(new CoverageServlet()),"/*");
	        context.setContextPath("/");
	        context.setResourceBase(".");
	        server.setHandler(context);
			
		    
		    
		    LOG.info("Starting server "+getProgramName()+" on port "+this.serverPort);
		    server.start();
		    LOG.info("Server started. Press Ctrl-C to stop");
		    server.join();
		    return 0;
			}
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}	


public static void main(final String[] args) throws Exception{
    new CoverageServer().instanceMainWithExit(args);
	}

}
