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
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
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
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
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
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## input

input is a set of indexed BAM file or a file with the suffix `.list` containing the path to the bams.
 
## Example

```
java -jar dist/coverageserver.jar \
	--pedigree fam.ped \
	--bed roi.bed \
	-o comments.bed \
	-R fasta src/test/resources/S*.bam

```
## Hidden parameters

 * `columns=5` change the number of columns at runtime.

## Screenshot

![https://twitter.com/yokofakun/status/1227932501747871745](https://pbs.twimg.com/media/EQp-Ga4XsAAxNYn?format=png&name=small)

![https://twitter.com/yokofakun/status/1228260742157209601](https://pbs.twimg.com/media/EQuooeGX0AAAHeu?format=jpg&name=medium)

END_DOC
 
 */

@Program(name="coverageserver",
	description="Jetty Based http server serving Bam coverage.",
	creationDate="20200212",
	modificationDate="20200214",
	keywords={"cnb","bam","coverage","server"}
	)
public  class CoverageServer extends Launcher {
	private static final Logger LOG = Logger.build(CoverageServer.class).make();
	
	@Parameter(names="--port",description="server port.")
	private int serverPort = 8080;
	@Parameter(names= {"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxRef = null;
	@Parameter(names= {"--bed","--B"},description="Optional bed file containing user's intervals. 4th column is used as the name of the interval")
	private Path intervalsource = null;
	@Parameter(names= {"--pedigree","-p"},description=PedigreeParser.OPT_DESC)
	private Path pedigreePath = null;
	@Parameter(names= {"--max_-size"},description="Security. Max interval size. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_window_size = 10_000_000;
	@Parameter(names= {"--link","--url","--hyperlink"},description=Hyperlink.OPT_DESC,converter=Hyperlink.StringConverter.class,splitter=NoSplitter.class)
	private Hyperlink hyperlink =Hyperlink.empty();
	@Parameter(names= {"--width"},description="Image width")
	private int image_width= 500;
	@Parameter(names= {"--height"},description="Image height")
	private int image_height= 300;
	@Parameter(names= {"--images-per-row","--ipr"},description="Number of images per row.")
	private int images_per_row= 2;
	@Parameter(names= {"--extend"},description="Extend interval by this factor. e.g: if x='0.5' chr1:100-200 -> chr1:50-250")
	private double extend_factor=1.0;
	@Parameter(names= {"-o","--output","--comment"},description="Output file for writing comments as a BED file.")
	private Path commentPath= null; 
	@Parameter(names= {"--mapq"},description="Min. Read Mapping Quality.")
	private int min_mapq = 0; 
	@Parameter(names= {"--small-length"},description="show reads if the region has a length <= 'x'.")
	private int small_region_size = 1_000; 
	@Parameter(names= {"--vcf","--region","--regions","--intervals"},description="Same as --bed but intervals won't be annotated. "+IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class)
	private IntervalListProvider intervalListProvider = IntervalListProvider.empty();

	
	
	private SAMSequenceDictionary dictionary;
	private final List<ReviewedInterval> named_intervals = new Vector<>();
	private final List<BamInput> bamInput = new Vector<>();
	private Pedigree pedigree = null;
	

	private static class ReviewedInterval extends SimpleInterval {
		final String name;
		boolean reviewed = false;
		ReviewedInterval(final Locatable loc,String name) {
			super(loc);
			this.name=StringUtils.ifBlank(name, "");
			}
		public String getName()  {
			return this.name;
			}
		}
	
	private static class BamInput {
		final Path bamPath;
		String sample;
		BamInput(final Path path) {
			this.bamPath = path;
		}
	}
	
	private static class Coverage {
		private final float array[];
		int count=0;
		Coverage(final int len) {
			this.array=new float[len];
			Arrays.fill(this.array, 0.0f);
			count=0;
			}
		double median() {
			if(count==0) {
				LOG.warn("count=0");
				return 1.0;
				}
			Arrays.sort(this.array,0,count);
			int mid_x = count/2;
			if(count%2==0) {
				return (array[mid_x-1]+array[mid_x])/2.0;
			} else {
				return array[mid_x];
			}
			
		
		}
		void add(final int v) {
			if(count>=this.array.length) {
				LOG.warn("array out of range");
				return;
			}
			this.array[this.count++]=(float)v;
		}
		
	}

	
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
			if(pathInfo.equals("/getimage"))
				{
				printImage(request,response);
				}
			else if(pathInfo.equals("/comment")) {
				saveComment(request,response);
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
	

	private void saveComment( HttpServletRequest request, HttpServletResponse response) 	throws IOException, ServletException{
		final SimpleInterval region = parseInterval(request.getParameter("interval"));
		final String comment =request.getParameter("comment");
		if(this.commentPath==null || region==null || StringUtils.isBlank(comment)) {
			LOG.info("region "+region+"."+request.getParameter("interval"));
			LOG.info("comment "+comment+".");
			response.reset();
			response.sendError(HttpStatus.BAD_REQUEST_400,"bad request");
			response.flushBuffer();
			return;
		}
		try(BufferedWriter pw= Files.newBufferedWriter(this.commentPath,StandardOpenOption.CREATE,StandardOpenOption.APPEND)) {
			pw.append(region.getContig());
			pw.append("\t");
			pw.append(String.valueOf(region.getStart()-1));
			pw.append("\t");
			pw.append(String.valueOf(region.getEnd()));
			pw.append("\t");
			pw.append(StringUtils.normalizeSpaces(comment));
			pw.append("\t");
			pw.append(StringUtils.ifBlank(this.named_intervals.
					stream().
					filter(I->new SimpleInterval(I).equals(region)).
					map(R->R.getName()).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.joining("; ")),"."));
			pw.append("\t");
			pw.append(StringUtils.now());
			pw.append('\n');
			pw.flush();
			
			response.setContentType("text/plain");
			PrintWriter w=response.getWriter();
			w.append("OK. Comment saved.");
			w.flush();
			return;
		} catch(final IOException err) {
			 LOG.error(err);
			 response.reset();
			 response.setContentType("text/plain");
			 PrintWriter w=response.getWriter();
			 w.append("Error. Cannot save comment ("+err.getMessage()+")");
			 w.flush();
			 return;
		}
		
	}
	
	private boolean acceptRead(final SAMRecord rec) {
		 if(
			rec.getReadUnmappedFlag() ||
			rec.getDuplicateReadFlag() ||
	 		rec.getReadFailsVendorQualityCheckFlag() ||
	 		rec.isSecondaryOrSupplementary() ||
	 		rec.getMappingQuality()<this.min_mapq) return false;
		 return true;
		}
	
	private void writeImage(
			final BufferedImage img,
			final BamInput bam,Locatable region,
			final HttpServletResponse response
			) throws IOException{

		 final String basename = bam.sample+"_"+region.getContig()+"_"+region.getStart()+"_"+region.getEnd();
		 response.setContentType("image/png");
		 response.addHeader("Content-Disposition","form-data; name=\""+basename+"\"; filename=\""+basename +".png\"");
		 try {
			 ImageIO.write(img, "PNG", response.getOutputStream());
			 response.flushBuffer();
		 	 }
		 catch(Throwable err) {
		 	}
		}
	
	/** print BAM for small interval, displaying reads */
	private void printRaster(final BamInput bam,final SimpleInterval midRegion,final SimpleInterval region,final HttpServletRequest request,final HttpServletResponse response) throws IOException, ServletException {
		final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*(double)image_width;
		final SamReaderFactory srf = SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).referenceSequence(this.faidxRef);
		final List<List<SAMRecord>> rows = new ArrayList<>();
		try(SamReader sr=srf.open(bam.bamPath)) {
			 try(CloseableIterator<SAMRecord> iter=sr.query(
					 region.getContig(),
					 Math.max(0,region.getStart()-this.small_region_size), //extend to get clip
					 region.getEnd()+this.small_region_size,false)) {
				 while(iter.hasNext()) {
					 final SAMRecord rec=iter.next();
					 if(!acceptRead(rec)) continue;
					 if(rec.getUnclippedEnd() < region.getStart()) continue;
					 if(rec.getUnclippedStart() > region.getEnd()) continue;
					 final Cigar cigar = rec.getCigar();
					 if(cigar==null || cigar.isEmpty()) continue;

					 int y=0;
					 for(y=0;y< rows.size();++y) {
						 final List<SAMRecord> row = rows.get(y);
						 final SAMRecord last = row.get(row.size()-1);
						 if(position2pixel.applyAsDouble(last.getUnclippedEnd()+1) +1  < position2pixel.applyAsDouble(rec.getUnclippedStart())) {
							 row.add(rec);
							 break;
						 	}
					 	}
				     if(y==rows.size()) {
						 final List<SAMRecord> row = new ArrayList<>();
						 row.add(rec);
						 rows.add(row);
				     	} 
				 	}
				 }//end iterator
			}//end samreder
		ReferenceSequence refInInterval=null;
		 try (ReferenceSequenceFile refseq=ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidxRef)) {
			 
			 final SAMSequenceRecord ssr = this.dictionary.getSequence(region.getContig());
			 if(region.getStart()<=ssr.getSequenceLength()) {
				 refInInterval = refseq.getSubsequenceAt(
						 region.getContig(),
						 region.getStart(),
						 Math.min(region.getEnd(),ssr.getSequenceLength())
						 );
			 	}
		 	}
		
		 final BufferedImage img = new BufferedImage(image_width, image_height, BufferedImage.TYPE_INT_RGB);
		 final Graphics2D g=img.createGraphics();
		 g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		 g.setColor(new Color(240,240,240));
		 g.fillRect(0, 0, image_width+1, image_height+1);
 
		 
	     
	     final Stroke oldStroke = g.getStroke();	 
	     g.setStroke(new BasicStroke(0.5f));
	     
	     g.setColor(Color.WHITE);
	     final double mid_start =  position2pixel.applyAsDouble(midRegion.getStart());
	     final double mid_end =  position2pixel.applyAsDouble(midRegion.getEnd()+1);
	     g.fill(new Rectangle2D.Double(mid_start, 0, (mid_end-mid_start),this.image_height));
	     

	     
	     final int int_coverage[]=new int[region.getLengthOnReference()];
	     final int margin_top=12;
	     
	     final double featureHeight= Math.min(20,(this.image_height-margin_top)/Math.max(1.0,(double)rows.size()));
	     
	     double y=image_height-featureHeight;
	     
	     for(final List<SAMRecord> row:rows) {
	    	final double h2= Math.min(featureHeight*0.9,featureHeight-2);

	    	for(final SAMRecord rec: row) {
	    		final Cigar cigar=rec.getCigar();
	    		if(cigar==null || cigar.isEmpty()) continue;
	    		
	    		/* fill coverage array  */
	    		int ref1 = rec.getAlignmentStart();
	    		for(CigarElement ce: cigar.getCigarElements()) {
	    			 if(ref1> region.getEnd()) break;
					 final CigarOperator op=ce.getOperator();
					 if(op.consumesReferenceBases()) {
						 if(op.consumesReadBases()) {
							 for(int x=0;x< ce.getLength();++x) {
								 int pos=ref1+x;
								 if(pos< region.getStart()) continue;
								 if(pos> region.getEnd()) break;
								 int_coverage[pos-region.getStart()]++;
							 }
						 }
						 ref1+=ce.getLength();
					 }
				 }
	    		
	    		/* draw rec itself */
	    		final double midy=y+h2/2.0;
	    		g.setColor(Color.DARK_GRAY);
	    		g.draw(new Line2D.Double(
	    				position2pixel.applyAsDouble(rec.getUnclippedStart()),
	    				midy,
	    				position2pixel.applyAsDouble(rec.getUnclippedEnd()),
	    				midy));
	    		ref1 = rec.getUnclippedStart();
	    		final List<Double> insertions = new ArrayList<>();
	    		for(final CigarElement ce: cigar.getCigarElements()) {
	    			if(ref1> region.getEnd()) break;
	    			final CigarOperator op=ce.getOperator();
	    			Shape shape = null;
	    			Color fill=null;
	    			Color stroke=Color.DARK_GRAY;
	    			switch(op) {
	    				case P: break;
	    				case M://through
	    				case X://through
	    				case EQ://through
	    				case S: //through
	    				case H: 
	    						final double x1=position2pixel.applyAsDouble(ref1);
	    						
	    						shape = new Rectangle2D.Double(
	    						x1, y,
	    						position2pixel.applyAsDouble(ref1+ce.getLength())-x1,h2
	    						);
	    						
	    						ref1+=ce.getLength();
	    						switch(op) {
	    							case H: case S: fill=Color.YELLOW;break;
	    							case X: fill=Color.RED;break;
	    							case EQ: case M: fill=Color.LIGHT_GRAY;break;
	    							default:break;
	    							}
	    						break;
	    				case N://through
	    				case D: shape=null;fill=null;stroke=null;ref1+=ce.getLength();break;
	    				case I: shape=null;fill=null;stroke=null;insertions.add(position2pixel.applyAsDouble(ref1));break;
	    				default: throw new IllegalStateException(""+op);
	    				}
	    			if(ref1 < region.getStart()) continue;
	    			
	    			if(shape!=null) {
	    				if(fill!=null) {g.setColor(fill);g.fill(shape);}
	    				if(stroke!=null && h2>4)  {g.setColor(stroke);g.draw(shape);}
	    				}
	    			} // end loop cigar
	    		
	    		
	    		 /* draw mismatched bases */
	   	     	if(refInInterval!=null && rec.getReadBases()!=null && rec.getReadBases()!=SAMRecord.NULL_SEQUENCE) {
	   	     		final byte bases[]=rec.getReadBases();
	   	     		final IntFunction<Character> baseRead= IDX-> IDX<0 || IDX>=bases.length || bases==SAMRecord.NULL_SEQUENCE?'N':(char)Character.toUpperCase(bases[IDX]);
	   	     		int read0=0;
	   	     		ref1 = rec.getAlignmentStart();
		   	     	for(CigarElement ce: cigar.getCigarElements()) {
		    			if(ref1> region.getEnd()) break;
		    			final CigarOperator op=ce.getOperator();
		    			switch(op) {
			    			case P:break;
			    			case H:break;
			    			case D: case N: ref1+=ce.getLength(); break;
			    			case S: case I: read0+=ce.getLength(); break;
			    			case EQ:case M: case X:
			    				{
			    				for(int j=0;j< ce.getLength();j++) {
			    					if(ref1+j< region.getStart()) continue;
			    					if(ref1+j>=region.getStart()+refInInterval.length()) break;
			    					final int ref_base_idx = ref1-region.getStart()+j;
			    					char ctgBase =(char)(ref_base_idx<0 || ref_base_idx>=refInInterval.length()?'N':Character.toUpperCase(refInInterval.getBases()[ref_base_idx]));
			    					if(ctgBase=='N') continue;
			    					char readBase = baseRead.apply(read0+j);
			    					if(readBase=='N') continue;
			    					if(readBase==ctgBase) continue;
			    					g.setColor(Color.ORANGE);
			    					final double x1 = position2pixel.applyAsDouble(ref1+j);
			    					final double x2 = position2pixel.applyAsDouble(ref1+j+1);
			    					g.fill( new Rectangle2D.Double( x1, y,x2-x1,h2));
			    					}
			    				read0+=ce.getLength();
			    				ref1+=ce.getLength();
			    				break;
			    				}
			    			default:break;
		    				}
		   	     		}
	   	     		
	   	     		}
	    		
	    		
	    		for(double px:insertions) {
	    			g.setColor(Color.RED);
	    			g.draw(new Line2D.Double(px,y-0.5,px,y+h2+0.5));
	    			}
	    		}
	    	y-=featureHeight;
	     	}
	    
	     
	     double max_cov = IntStream.of(int_coverage).max().orElse(0);

	     g.setColor(Color.DARK_GRAY);
	     g.drawString("Sample:"+ bam.sample +" max-cov:"+(int)max_cov+" "+region.toNiceString() , 10, 10);

	     
	     /* plot coverage */
	     final GeneralPath gp = new GeneralPath();
	     for(int i=0;max_cov>0 && i< int_coverage.length;++i) {
	    	 final double x1= position2pixel.applyAsDouble(region.getStart()+i);
	    	 final double x2= position2pixel.applyAsDouble(region.getStart()+i+1);
	    	 final double y1= image_height - (int_coverage[i]/max_cov)*(image_height-margin_top);
	    	 if(i==0) gp.moveTo(x1, y1);
	    	 else gp.lineTo(x1, y1);
	    	 gp.lineTo(x2, y1);
	     	}
	     g.setStroke(new BasicStroke(0.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL, 0f, new float[] { 1f, 2f, 1f }, 0f));
	     g.setColor(Color.BLUE);
	     g.draw(gp);
	     
	     g.setStroke(oldStroke);
	     
	     g.setColor(Color.PINK);
	     g.draw(new Line2D.Double(mid_start,0,mid_start,image_height));
	     g.draw(new Line2D.Double(mid_end,0,mid_end,image_height));
	     
	     
	     
	     writeImage(img,bam,region,response);
		}
	
	private void printImage(final HttpServletRequest request,final HttpServletResponse response) throws IOException, ServletException
	{
		int bam_id;
		try {
			bam_id = Integer.parseInt(StringUtils.ifBlank(request.getParameter("id"),"-1"));
		} catch(Exception err) {
			bam_id=-1;
		}
		final SimpleInterval midRegion = parseInterval(request.getParameter("interval"));
		if(midRegion==null || bam_id<0 || bam_id>=this.bamInput.size()) {
			response.reset();
			response.sendError(HttpStatus.BAD_REQUEST_400,"id:"+bam_id);
			response.flushBuffer();
			return;
		}
		final int extend = (int)(midRegion.getLengthOnReference()*this.extend_factor);
		int xstart = Math.max(midRegion.getStart()-extend,0);
		int xend = midRegion.getEnd()+extend;
		final SAMSequenceRecord ssr = this.dictionary.getSequence(midRegion.getContig());
		if(ssr!=null) {
			xend = Math.min(xend, ssr.getSequenceLength());
		}
		final SimpleInterval region = new SimpleInterval(midRegion.getContig(),xstart,xend);
		if(region.getLengthOnReference()>this.max_window_size)  {
			response.reset();
			response.sendError(HttpStatus.BAD_REQUEST_400,"contig:"+midRegion);
			response.flushBuffer();
			return;
		}
		
		final BamInput bam = this.bamInput.get(bam_id);

		if(region.length() <=this.small_region_size) {
			printRaster(bam,midRegion, region, request, response);
			return;
		}
		
		final boolean normalize = request.getParameter("normalize")!=null;
		
		final SamReaderFactory srf = SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).referenceSequence(this.faidxRef);
		try(SamReader sr=srf.open(bam.bamPath)) {
			
			
			 final int int_coverage[]=new int[region.getLengthOnReference()];
			 Arrays.fill(int_coverage, 0);
			 try(CloseableIterator<SAMRecord> iter=sr.query(region.getContig(), region.getStart(), region.getEnd(),false)) {
				 while(iter.hasNext()) {
					 final SAMRecord rec=iter.next();
					 if(!acceptRead(rec)) continue;
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
									 int_coverage[pos-region.getStart()]++;
								 }
							 }
							 ref+=ce.getLength();
						 }
					 }
				 }
			 }
			 /* smooth coverage */
			 if(int_coverage.length>image_width) {
				 final int copy[]=Arrays.copyOf(int_coverage, int_coverage.length);
				 final int len = Math.max(1,int_coverage.length/100);
				 
				 for(int i=0;i< int_coverage.length;i++) {
					 int j=Math.max(0, i-len);
					 double sum=0;
					 int count=0;
					 while(j< i+len && j< copy.length) {
						 sum +=copy[j];
						 j++;
						 count++;
					 }
					 int_coverage[i]=(int)(sum/count);
				 }
			 }
			 
			final double norm_coverage[] = new double[int_coverage.length];
			final double median;
			/* normalize on median */
			if(normalize)
				{
				final Coverage leftrightcov = new Coverage( extend*2 );
				 for(int x=region.getStart();x<midRegion.getStart();x++) {
						final int idx = x-region.getStart();
						leftrightcov.add(int_coverage[idx]);
					}
				 for(int x=midRegion.getEnd()+1;x<=region.getEnd();x++) {
						final int idx = x-region.getStart();
						leftrightcov.add(int_coverage[idx]);
					}
				 
				median = Math.max(1.0,leftrightcov.median());
				//LOG.info("median is "+median+" "+leftrightcov.median());
				for(int x=0;x< int_coverage.length;++x) {
					norm_coverage[x]=int_coverage[x]/median;
					}
				} 
			else /* no normalisation */
			
				{
				/* won't be used */
				median = Double.NaN;
				for(int x=0;x< int_coverage.length;++x) {
						norm_coverage[x]=int_coverage[x];
					}
				}
			
			 final double real_max_cov = DoubleStream.of(norm_coverage).max().orElse(1.0);
			 final double max_cov= Math.max((normalize?2:10),real_max_cov );
			 final double pixelperbase = image_width/(double)norm_coverage.length;
			 final BufferedImage img = new BufferedImage(image_width, image_height, BufferedImage.TYPE_INT_RGB);
			 final Graphics2D g=img.createGraphics();
			 g.setColor(Color.WHITE);
			 g.fillRect(0, 0, image_width+1, image_height+1);
			 
			 final Composite oldComposite = g.getComposite();
			 //g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)Math.min(1.0, pixelperbase)));
			 
			 for(int x=0;x< norm_coverage.length;++x) {
				 final double height = image_height*(norm_coverage[x]/max_cov);
				
				 if(normalize) g.setColor(Color.DARK_GRAY);
				 else if(max_cov<10) g.setColor(Color.RED);
				 else if(max_cov<20) g.setColor(Color.BLUE);
				 else g.setColor(Color.DARK_GRAY);
				
				 g.fill(new Rectangle2D.Double(
						 x*pixelperbase,
						 image_height-height,
						 pixelperbase,
						 height));
			 	}
			 g.setComposite(oldComposite);
			 
			 g.setColor(Color.DARK_GRAY);
			 g.drawString("max-cov:"+IntStream.of(int_coverage).max().orElse(0)+
					 (normalize?" normalized on median ("+median+")":"")+
					 " sample:"+ bam.sample +" "+
					 region.toNiceString()
					 , 10, 10);

			 /* ticks for vertical axis */
			 g.setColor(Color.MAGENTA);
			 for(int i=1;i<10;i++) {
				 double cov=max_cov/10.0*i;
				 if(!normalize) cov= Math.ceil(cov);
				 final double y = image_height - image_height/10.0*i;
				 if(!normalize && i>0 && (int)cov==Math.ceil(max_cov/10.0*(i-1))) continue;
				 g.drawLine(0, (int)y, 5, (int)y);
				 g.drawString(normalize?String.format("%.2f",cov):String.valueOf((int)cov),7,(int)y);
			 }
			 
			 /* vertical line for original view */
			 g.setColor(Color.PINK);
			 double vertical = ((midRegion.getStart()-region.getStart())/(double)region.getLengthOnReference())*image_width;
			 g.draw(new Line2D.Double(vertical, 0, vertical, image_height));
			 vertical = ((midRegion.getEnd()-region.getStart())/(double)region.getLengthOnReference())*image_width;
			 g.draw(new Line2D.Double(vertical, 0, vertical, image_height));

			 if(normalize) {
				 /* horizontal line for median 0.5 / 1 / 1.5 */
				 for(int t=1;t<4;++t) {
					 g.setColor(t==2?Color.ORANGE:Color.PINK);
					 final double mediany= image_height-((0.5*t)/max_cov)*image_height;
					 g.draw(new Line2D.Double(0,mediany,image_width,mediany));
					 }
				 }

			 
			 g.setColor(Color.GRAY);
			 g.drawRect(0, 0, img.getWidth(),  img.getHeight());
			 
			 writeImage(img,bam,region,response);
			}
		}
	
	/** write generic information for a sample */
	private void writeSample(final XMLStreamWriter w,final Sample sample) throws XMLStreamException {
		if(sample==null) return;
		w.writeCharacters(" ");
		switch(sample.getSex()){
			case male: w.writeEntityRef("#9794");break;
			case female:  w.writeEntityRef("#9792");break;
			default: break;
			}
		w.writeCharacters(" ");
		
		switch(sample.getStatus()) {
			case unaffected:  w.writeEntityRef("#128578"); break;
			case affected:  w.writeEntityRef("#128577"); break;
			default: break;
			}
		w.writeCharacters(" ");
		}
	
	/** print HTML page */
	private void printPage(final HttpServletRequest request,final HttpServletResponse response)	throws IOException, ServletException
		{
		 String message="";
		 
		 final String charset = StringUtils.ifBlank(request.getCharacterEncoding(), "UTF-8");
		 response.setContentType("text/html; charset="+charset.toLowerCase());
		 response.setCharacterEncoding(charset);
		 PrintWriter pw = response.getWriter();
		 SimpleInterval interval = null;
		 
		 int columns_count = this.images_per_row;
		 if(!StringUtils.isBlank(request.getParameter("columns"))) {
			 columns_count = Math.max(1,StringUtils.parseInt(request.getParameter("columns")).orElse(this.images_per_row));
		 }
		 
		 final boolean normalize = request.getParameter("normalize")!=null;
		 
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
		 
		 if(interval==null && 	!this.named_intervals.isEmpty()) {
			 /* first non reviewed */
			 interval = this.named_intervals.stream().
					 filter(R->!R.reviewed).
					 findFirst().
					 map(R->new SimpleInterval(R)).
					 orElse(null);
			 /* all reviewed ? */
			 if(interval==null) interval = new SimpleInterval(this.named_intervals.get(0));
		 }
		 
		 /* still no interval ?*/
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
		 
		 /* ZOOM buttons */
		 for(int side=0;side<2;++side) {
			 final String param = "zoom"+(side==0?"in":"out");
			 String value = request.getParameter(param);
			 if(StringUtils.isBlank(value)) continue;
			 double factor = Double.parseDouble(value);
			 if(side==0) factor=1.0/factor;
			 int length0 = interval.getLengthOnReference();
			 int length1 = (int)(length0*factor);
			 if(length1< 1) length1=1;
			 if(length1> this.max_window_size) length1=this.max_window_size;
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
			
			w.writeEmptyElement("meta");
			w.writeAttribute("charset", charset);
			
			w.writeStartElement("title");
			w.writeCharacters(title);
			w.writeEndElement();//title
			
			w.writeStartElement("style");
			w.writeCharacters(
					"body {background-color:#f0f0f0;color:#070707;font-size:18px;}"
					+ "h1 {text-align:center;color:#070707;}"
					+ ".grid-container {display: grid; grid-template-columns: "+IntStream.range(0, Math.max(1,columns_count)).mapToObj(X->"auto").collect(Collectors.joining(" "))+";grid-gap: 10px;  padding: 10px;}"
					+ ".span1 {border:1px dotted blue;}"
					+ ".lbl {font-weight: bold;}"
					+ ".bampath {font-family:monospace;font-size:12px; font-style: normal;color:gray;}"
					+ ".headerform {background-color:lavender;text-align:center;font-size:14px;}"
					+ ".comment {background-color:khaki;text-align:center;font-size:14px;}"
					+ ".highlight {background-color:#DB7093;}"
					+ ".parents {font-size:75%;color:gray;}"
					+ ".allsamples {font-size:125%;}"
					+ ".message {color:red;}"
					+ ".affected {background-color:#e6cccc;}"
					);
			w.writeEndElement();//title
			
			w.writeStartElement("script");
			{
				w.writeCharacters(""); 
				w.flush();
				/* messy with < and > characters + xmlstream */
				pw.write(
				"function highlight(name) {"+
				"var e = document.getElementById(name+\"_div\"); if(e==null) return;"+
				"e.classList.add(\"highlight\");"+
				"setTimeout(function () {e.classList.remove(\"highlight\");},2000);"+
				"}"+
				"function sendComment() {"+
				"var commentstr=document.getElementById(\"comment\").value;" +
				"console.log(\"comment is\"+commentstr);"+
				"if(commentstr.trim().length==0) return;" +
				"var xhttp = new XMLHttpRequest();"+
				"xhttp.onreadystatechange = function() {if (this.readyState == 4) if(this.status == 200) { alert(this.responseText);} else alert(\"Cannot send Comment.\")};"+
				"xhttp.open(\"GET\",\"/comment?interval="+StringUtils.escapeHttp(interval.toString())+"&comment=\"+encodeURI(commentstr), true);"+
				"xhttp.send();"+
				"}"+
				"function loadImage(idx) {"+
				"if(idx>="+this.bamInput.size()+") return;"+
				"var img = document.getElementById(\"bamid\"+idx);"+
				"img.addEventListener('load',(event) => {img.width="+image_width+";img.height="+image_height+";loadImage(idx+1);});"+
				"img.setAttribute(\"src\",\"/getimage?id=\"+idx+\"&interval="+ StringUtils.escapeHttp(interval.toString()) +(normalize?"&normalize=1":"")+"\");"+
				"img.setAttribute(\"alt\",\"bam idx\"+idx);"+
				"}"+
				"function init() {"+
				"var span=document.getElementById(\"spaninterval\");"+
				"if(span!=null) span.addEventListener('click',(evt)=>{document.getElementById(\"custominterval\").value = span.textContent; });"+
				"var sel=document.getElementById(\"selinterval\");"+
				"if(sel!=null) sel.addEventListener('change',(evt)=>{document.getElementById(\"custominterval\").value = evt.target.value; });"+
				"var cbox=document.getElementById(\"norminput\");"+
				"if(cbox!=null) cbox.addEventListener('change',(evt)=>{cbox.form.submit(); });"+
				"var comment=document.getElementById(\"commentbtn\");" +
				"if(comment!=null) comment.addEventListener('click',(evt)=>{ console.log(\"send comment\");sendComment(); });"+
				"var shortcuts=document.getElementById(\"shortcuts\");"+
				"if(shortcuts!=null) shortcuts.addEventListener('change',(evt)=>{document.getElementById(\"comment\").value += evt.target.value; });"+
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
			
			w.writeCharacters(SequenceDictionaryUtils.getBuildName(this.dictionary).orElse("")+" ");
			
			
			final String outlinkurl = hyperlink.apply(interval);
			if(StringUtils.isBlank(outlinkurl)) {
				w.writeCharacters(title);
			} else {
				w.writeStartElement("a");
				w.writeAttribute("href", outlinkurl);
				w.writeAttribute("target", "_blank");
				w.writeAttribute("title", title);
				w.writeCharacters(title);
				w.writeEndElement();
				}
			w.writeEndElement();//h1
			
			if(!StringUtils.isBlank(message)) {
				w.writeStartElement("h2");
				w.writeAttribute("class", "message");
				w.writeCharacters(message);
				w.writeEndElement();//h2
			}
			
			
			w.writeEmptyElement("a");
			w.writeAttribute("name", "top");

			
			w.writeComment("BEGIN FORM");
			w.writeStartElement("div");
			w.writeAttribute("class", "headerform");
			w.writeStartElement("form");
			w.writeAttribute("method", "GET");
			w.writeAttribute("action", "/page");
			
			w.writeEmptyElement("input");
			w.writeAttribute("name", "interval");
			w.writeAttribute("value", interval.toString());
			w.writeAttribute("type", "hidden");

			/* number of images per row */
			if(!StringUtils.isBlank(request.getParameter("columns"))) {
				w.writeEmptyElement("input");
				w.writeAttribute("name", "columns");
				w.writeAttribute("value",String.valueOf(columns_count));
				w.writeAttribute("type", "hidden");
			 }

			
			
			
			/* write select box with predefined interval */
			if(!this.named_intervals.isEmpty()) {
				w.writeStartElement("select");
				w.writeAttribute("class", "btn");
				w.writeAttribute("id", "selinterval");
				
				w.writeEmptyElement("option");
				
				for(int side=0;side<2;++side) {
					for(final ReviewedInterval r:this.named_intervals) {
						if(side==0 && r.reviewed) continue;
						if(side==1 && !r.reviewed) continue;
						final SimpleInterval simple= new SimpleInterval(r);
						w.writeStartElement("option");
						w.writeAttribute("value",simple.toString());
						if(simple.equals(interval)) {
							r.reviewed=true;
							w.writeAttribute("selected", "true");
						}
						if(r.reviewed) {
							w.writeEntityRef("#x2713");
							w.writeCharacters(" ");
							}
						w.writeCharacters(simple.toNiceString()+(StringUtils.isBlank(r.getName())?"":" ["+r.getName()+"]"));
						w.writeEndElement();
						w.writeCharacters("\n");
						}
					}
				
				w.writeEndElement();//select
				}
				w.writeEmptyElement("br");
				
				//
				
				/* move */
				w.writeStartElement("label");
				w.writeAttribute("class", "lbl");
				w.writeCharacters("move");
				w.writeEndElement();
				for(final String mv:new String[]{"<<<","<<","<",">",">>",">>>"}) {
					w.writeEmptyElement("input");
					w.writeAttribute("class", "btn");
					w.writeAttribute("type", "submit");
					w.writeAttribute("name", "move");
					w.writeAttribute("value", mv);
					}
				
				/* zoom in */
				w.writeStartElement("label");
				w.writeAttribute("class", "lbl");
				w.writeCharacters("zoom in");
				w.writeEndElement();
				for(final String zoom:new String[]{"1.5","3","10"}) {
					w.writeEmptyElement("input");
					w.writeAttribute("class", "btn");
					w.writeAttribute("type", "submit");
					w.writeAttribute("name", "zoomin");
					w.writeAttribute("value", zoom);
				}
				
				
				/* zoom in */
				w.writeStartElement("label");
				w.writeAttribute("class", "lbl");
				w.writeCharacters("zoom out");
				w.writeEndElement();
				for(final String zoom:new String[]{"1.5","3","10","100"}) {
					w.writeEmptyElement("input");
					w.writeAttribute("type", "submit");
					w.writeAttribute("class", "btn");
					w.writeAttribute("name", "zoomout");
					w.writeAttribute("value", zoom);
					}
				/* checkbox normalize */
				
				w.writeEmptyElement("input");
				w.writeAttribute("id", "norminput");
				w.writeAttribute("name", "normalize");
				w.writeAttribute("type", "checkbox");
				if(normalize) w.writeAttribute("checked", "true");
				
				w.writeStartElement("label");
				w.writeAttribute("class", "lbl");
				w.writeAttribute("for", "norminput");
				w.writeCharacters("normalize");
				w.writeEndElement();
				
				
				w.writeEmptyElement("br");
				
				w.writeStartElement("span");
				w.writeAttribute("class", "span1");
				w.writeAttribute("id", "spaninterval");
				w.writeCharacters(interval.toNiceString());
				w.writeEndElement();//span
				
				w.writeEmptyElement("input");
				w.writeAttribute("name", "custominterval");
				w.writeAttribute("id", "custominterval");
				w.writeAttribute("type", "text");
				w.writeAttribute("value","");
				
				w.writeStartElement("button");
				w.writeAttribute("class", "btn");
				w.writeAttribute("name", "parseinterval");
				w.writeCharacters("GO");
				w.writeEndElement();//button
				

				
			w.writeEndElement();//form
			w.writeEndElement();//div
			
			if(this.commentPath!=null) {
				w.writeStartElement("div");
				w.writeAttribute("class", "comment");
				
				w.writeStartElement("label");
				w.writeAttribute("for", "comment");
				w.writeCharacters("Comment:");
				w.writeEndElement();
				w.writeEmptyElement("input");
				w.writeAttribute("id", "comment");
				w.writeAttribute("type", "text");
				w.writeAttribute("placeholder", "Comment about "+ interval.toNiceString());
				
				w.writeStartElement("button");
				w.writeAttribute("class", "btn");
				w.writeAttribute("id", "commentbtn");
				w.writeCharacters("Send comment");
				w.writeEndElement();//button
				w.writeCharacters(" ");
				
				w.writeStartElement("select");
				w.writeAttribute("class", "btn");
				w.writeAttribute("id", "shortcuts");
				
				w.writeEmptyElement("option");
				
				for(final String opt : new String[] {"OK","BAD","Noisy","False positive","Deletion","Later"}) {
					w.writeStartElement("option");
					w.writeAttribute("value"," "+opt+".");
					w.writeCharacters(opt);
					w.writeEndElement();
				}
				
				w.writeEndElement();//select
				
				w.writeEndElement();//div
			}
			
			w.writeComment("END FORM");
			
			
			
			
			
			/* write anchors to samples */
			w.writeStartElement("div");
			w.writeAttribute("class", "allsamples");

			w.writeCharacters("Samples: ");
			for(final BamInput bi:this.bamInput.stream().sorted((A,B)->A.sample.compareTo(B.sample)).collect(Collectors.toList())) {
				w.writeStartElement("a");
				w.writeAttribute("title", bi.sample);
				w.writeAttribute("href", "#"+bi.sample);
				w.writeAttribute("onclick", "highlight('"+bi.sample+"');");
				w.writeCharacters("["+bi.sample);
				if(pedigree!=null) {
					final Sample sn = this.pedigree.getSampleById(bi.sample);
					if(sn!=null && sn.isAffected()){
						w.writeCharacters(" ");
						w.writeEntityRef("#128577"); 
					}
				}
				
				w.writeCharacters("] ");
				w.writeEndElement();
			}
			w.writeEndElement();//div
			w.writeCharacters("\n");
			
			w.writeStartElement("div");
			w.writeAttribute("class", "grid-container");
			
			for(int i=0;i< this.bamInput.size();i++) {
				final BamInput bamInput = this.bamInput.get(i);
				final Path bam = bamInput.bamPath;
				final Sample sample = this.pedigree==null?null:this.pedigree.getSampleById(bamInput.sample);


				w.writeStartElement("div");
				w.writeAttribute("id", bamInput.sample+"_div");
				w.writeAttribute("class", "sample"+(sample!=null && sample.isAffected()?" affected":""));
				w.writeEmptyElement("a");
				w.writeAttribute("name", bamInput.sample);

				w.writeStartElement("h3");
				w.writeCharacters(bamInput.sample);
				w.writeCharacters(" ");
				w.writeStartElement("a");
				w.writeAttribute("href","#top");
				w.writeAttribute("title","top");
				w.writeCharacters("[^]");
				w.writeEndElement();
				
				
					{
					if(sample!=null) {
						writeSample(w,sample);
						w.writeStartElement("span");
						w.writeAttribute("class","parents");
						for(int p=0;p<2;++p) {
							if(p==0 && !sample.hasFather()) continue;
							if(p==1 && !sample.hasMother()) continue;
							final Sample parent = (p==0?sample.getFather():sample.getMother());
							
							boolean has_bam= this.bamInput.stream().anyMatch(B->B.sample.equals(parent.getId()));
							w.writeCharacters(" ");
							w.writeCharacters(p==0?"Father ":"Mother ");
							if(has_bam) {
								w.writeStartElement("a");
								w.writeAttribute("href","#"+parent.getId());
								w.writeAttribute("title",parent.getId());
								w.writeAttribute("onclick", "highlight('"+parent.getId()+"');");
								w.writeCharacters("["+parent.getId()+"].");
								w.writeEndElement();
								}
							else
								{
								w.writeCharacters(parent.getId());
								}
							writeSample(w,parent);
							}
						w.writeEndElement();
						
					}
				}
				
				w.writeEndElement();
				
				w.writeStartElement("div");
				w.writeAttribute("class", "bampath");
				w.writeCharacters(bam.toString());
				w.writeEndElement();

				
				w.writeEmptyElement("img");
				w.writeAttribute("id","bamid"+i);
				w.writeAttribute("src", "https://upload.wikimedia.org/wikipedia/commons/d/de/Ajax-loader.gif");
				w.writeAttribute("width","32");
				w.writeAttribute("height","32");

				
				w.writeEndElement();//div
				w.writeCharacters("\n");
				}
			w.writeEndElement();
			
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
		if(this.image_width<10) {
			LOG.error("low image width");
			return -1;
			}
		if(this.image_height<10) {
			LOG.error("low image height");
			return -1;
			}
		if(this.images_per_row<1) {
			LOG.error("low images_per_row");
			return -1;
			}
		if(this.extend_factor <= 0) {
			LOG.error("bad extend_factor "+this.extend_factor);
			return -1;
			}
		try {
			
			this.bamInput.addAll(IOUtils.unrollPaths(args).stream().map(F->new BamInput(F)).collect(Collectors.toList()));
			if(this.bamInput.isEmpty()) {
				LOG.error("No BAM defined.");
				return -1;
				}
			this.dictionary = SequenceDictionaryUtils.extractRequired(this.faidxRef);
			
			for(final BamInput bi:this.bamInput) {
				final SamReaderFactory srf = SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).referenceSequence(this.faidxRef);
				try(SamReader sr=srf.open(bi.bamPath)) {
					final SAMFileHeader header = sr.getFileHeader();
					SequenceUtil.assertSequenceDictionariesEqual(this.dictionary, SequenceDictionaryUtils.extractRequired(header));
					bi.sample = header.getReadGroups().
					 stream().
					 map(R->R.getSample()).filter(S->!StringUtils.isBlank(S)).
					 findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(bi.bamPath));
				}
			}
			
			if(this.pedigreePath!=null) {
				this.pedigree = new PedigreeParser().parse(this.pedigreePath);
			}
			
			if(this.intervalsource!=null) {
				final ContigNameConverter cvt = ContigNameConverter.fromOneDictionary(this.dictionary);
				final BedLineCodec codec = new BedLineCodec();
				try(BufferedReader br=IOUtils.openPathForBufferedReading(this.intervalsource)) {
					br.lines().
						filter(L->!BedLine.isBedHeader(L)).
						map(L->codec.decode(L)).
						filter(B->B!=null).
						map(B->new ReviewedInterval(new SimpleInterval(B.getContig(), B.getStart(), B.getEnd()),B.getOrDefault(3, ""))).
						map(B->{
							final String ctg= cvt.apply(B.getContig());
							if(StringUtils.isBlank(ctg)) return null;
							if(ctg.equals(B.getContig())) return B;
							return new ReviewedInterval(new SimpleInterval(ctg, B.getStart(), B.getEnd()),B.getName());
							}).
						filter(B->B!=null).
						forEach(B->named_intervals.add(B));
					}
				}
			
			this.intervalListProvider.
				dictionary(this.dictionary).
				skipUnknownContigs().
				stream().
				map(L->new Interval(L)).
				forEach(B->named_intervals.add(new ReviewedInterval(B,"")));
			
			final Server server = new Server(this.serverPort);
			
			final ServletContextHandler context = new ServletContextHandler();
	        context.addServlet(new ServletHolder(new CoverageServlet()),"/*");
	        context.setContextPath("/");
	        context.setResourceBase(".");
	        server.setHandler(context);
			
		    
		    
		    LOG.info("Starting server "+getProgramName()+" on port "+this.serverPort);
		    server.start();
		    LOG.info("Server started. Press Ctrl-C to stop. http://localhost:"+this.serverPort+"/coverage");
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
