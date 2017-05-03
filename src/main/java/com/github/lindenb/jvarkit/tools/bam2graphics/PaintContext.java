package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.Dimension;
import java.awt.Insets;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;

@Deprecated
@Program(name="paintcontext",description="")
public class PaintContext extends Launcher
	{
	private static final Logger LOG = Logger.build(PaintContext.class).make();

	
	private static final String XLINK=com.github.lindenb.jvarkit.util.ns.XLINK.NS;
	private static final String REPLACE="__REGION__";
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private List<BamFile> samFileReaders=new ArrayList<BamFile>();
	private TabixVcfFileReader tabixReader=null;
	private TabixFileReader knownGenes=null;
	private Insets insets=new Insets(100, 100, 100, 100);
	private int screenWidth=2000;
	private int gcWindowSize=5;
	@Parameter(names="-o",description="output file")
	private String fileOutPattern=null;
	
	private Hershey hershey=new Hershey();
	private final int TRANSCRIPT_HEIGHT=30;
	private final int GC_PERCENT_HEIGHT=100;
	private final int BAM_COVERAGE_HEIGHT=100;
	private final int SAMPLE_HEIGHT=50;

	
	private class BamFile
		{
		File file;
		SamReader sfr=null;
		public void close()
			{
			CloserUtil.close(sfr);
			}
		}
	

	
	
	private static class Interval
		{
		String chrom;
		int start;
		int end;
		public String getChrom()
			{	
			return chrom;
			}
		public int getStart()
			{
			return start;
			}
		public int getEnd()
			{
			return end;
			}
		private int distance()
			{
			return this.getEnd()-this.getStart();
			}
		@Override
		public String toString() {
			return chrom+":"+start+"-"+end;
			}
		}
	
	
	/* current interval */
	private Interval interval;
	
	
	
	
	private PaintContext()
		{
		
		}
	
	
	
	
	private double baseToPixel(int pos0)
		{
		return  ((pos0 - this.interval.getStart())/(double)this.interval.distance())*(this.screenWidth-(this.insets.right+this.insets.left) )
				;
		}
	
	
	private void writeBamSection(
			XMLStreamWriter w
			) throws XMLStreamException
		{
		final int drawinAreaWidth= screenWidth-(this.insets.right+this.insets.left);
		int y=0;
		for(BamFile bf:this.samFileReaders)
			{
			w.writeStartElement("g");
			w.writeAttribute("transform", "translate(0,"+y+")");
			int maxDepth=1;
			int depths[]=new int[interval.distance()];
			
			SAMRecordIterator iter=null;
			try
				{
				iter=bf.sfr.query(
						genomicSequence.getChrom(),
						interval.start+1,
						interval.end,
						true
						);
				while(iter.hasNext())
					{
					SAMRecord rec=iter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getNotPrimaryAlignmentFlag()) continue;
					if(rec.getDuplicateReadFlag()) continue;
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					for(int i1= Math.max(interval.start+1,rec.getAlignmentStart());
							i1 <= rec.getAlignmentEnd() && i1 <interval.end;
							++i1)
						{
						int index= (i1-interval.start)-1;
						if(index<0 || index > depths.length)
							{
							System.err.println("boum "+index);
							continue;
							}
						
						depths[index]++;
						maxDepth=Math.max(maxDepth, depths[index]);
						}
					}
				}
			catch(Exception err)
				{
				LOG.error(err);
				throw new RuntimeException(err);
				}
			finally
				{
				CloserUtil.close(iter);
				}
			
			
			w.writeAttribute("title","BAM");
			
			double y_array[]=new double[drawinAreaWidth];
			for(int pixelX=0;pixelX <y_array.length;++pixelX)
				{
				int base_0=(int)(((pixelX+0)/(double)drawinAreaWidth)*this.interval.distance());
				int base_1=(int)(((pixelX+1)/(double)drawinAreaWidth)*this.interval.distance());
				if(base_1==base_0) base_1++;
				
				double d=0;
				for(int i=base_0;i< base_1;++i)
					{
					d+=depths[i];
					}
				d=d/(base_1-base_0);
				d=d/(double)maxDepth;
				if(d>1.0) d=1.0;

				double y2=BAM_COVERAGE_HEIGHT-(float)BAM_COVERAGE_HEIGHT*(d);
				y_array[pixelX]=y2;
				}
			
			List<Point2D.Double> points=new ArrayList<Point2D.Double>();
			points.add(new Point2D.Double(0,BAM_COVERAGE_HEIGHT));//left,bottom
			for(int pixelX=0;pixelX <y_array.length;++pixelX)
				{
				if(pixelX>0 && pixelX+1 !=y_array.length)
					{
					if(y_array[pixelX]==y_array[pixelX+1]) continue;
					}
				Point2D.Double pt=new Point2D.Double(pixelX,y_array[pixelX]);
				points.add(pt);
				}
			points.add(new Point2D.Double(drawinAreaWidth,BAM_COVERAGE_HEIGHT));//right,bottom
			points.add(new Point2D.Double(0,BAM_COVERAGE_HEIGHT));//left,bottom

			StringBuilder sw=new StringBuilder();
			for(Point2D.Double pt:points) sw.append((int)pt.getX()).append(",").append(pt.getY()).append(" ");
			w.writeEmptyElement("polygon");
			w.writeAttribute("class","coverage");
			w.writeAttribute("points",sw.toString().trim());
			
			//label
			int font_size=10;
			String label=String.format("%10s",bf.file.getName());
			w.writeEmptyElement("path");
			w.writeAttribute("title",bf.file.getPath());
			w.writeAttribute("style","stroke:black;");
			w.writeAttribute("d",this.hershey.svgPath(
					label,
					 -Math.min(insets.left,label.length()*font_size),
					0,
					Math.min(insets.left,label.length()*font_size),
					font_size)
					);

			
			//axis
			
			for(int step=1;step<=10;++step)
				{	
				double y1=BAM_COVERAGE_HEIGHT-(BAM_COVERAGE_HEIGHT/10.0)*step;
				double x1=drawinAreaWidth+5;
				w.writeEmptyElement("line");
				w.writeAttribute("class", "xaxis");
				w.writeAttribute("x1",String.valueOf(0));
				w.writeAttribute("y1",String.valueOf(y1));
				w.writeAttribute("x2",String.valueOf(x1));
				w.writeAttribute("y2",String.valueOf(y1));
				w.writeEmptyElement("path");
				label=String.valueOf(maxDepth/10.0*step);
				w.writeAttribute("d",this.hershey.svgPath(
						label,
						x1,
						y1-font_size/2,
						label.length()*font_size,
						font_size)
						);
				}
			
			//frame
			w.writeEmptyElement("rect");
			w.writeAttribute("class","frame");

			w.writeAttribute("x",String.valueOf(0));
			w.writeAttribute("y",String.valueOf(0));
			w.writeAttribute("width",String.valueOf(drawinAreaWidth));
			w.writeAttribute("height",String.valueOf(BAM_COVERAGE_HEIGHT));
			
			
			
			w.writeEndElement();//g
			y+=BAM_COVERAGE_HEIGHT;
			}
		}
	
	
	private void writeKownGeneSection(
			XMLStreamWriter w,
			List<KnownGene> operon
			) throws XMLStreamException
		{
		final int drawinAreaWidth= screenWidth-(this.insets.right+this.insets.left);

		if(operon.isEmpty()) return;
		int y=0;
		for(KnownGene g:operon)
			{
			int cdsHeigh= 5;
			int exonHeight=TRANSCRIPT_HEIGHT-5;
			int midY=TRANSCRIPT_HEIGHT/2;
	
			w.writeStartElement("g");
			
			
			
			w.writeAttribute("transform", "translate(0,"+y+")");
			w.writeAttribute("clip-path","url(#kgclip)");
			w.writeAttribute("title", g.getName());
			
			/* transcript line */
			w.writeEmptyElement("line");
			w.writeAttribute("class","kgtr");
			w.writeAttribute("x1",String.valueOf(baseToPixel(trim(g.getTxStart()))));
			w.writeAttribute("y1",String.valueOf(midY));
			w.writeAttribute("x2",String.valueOf(baseToPixel(trim(g.getTxEnd()))));
			w.writeAttribute("y2",String.valueOf(midY));
			
			
			
			
			/* strand symbols */
			for(double pixX=0;
				pixX< drawinAreaWidth;
				pixX+=30)
				{
				double pos0= interval.start+(pixX/(double)drawinAreaWidth)*interval.distance();
				if(pos0< g.getTxStart()) continue;
				if(pos0> g.getTxEnd()) break;
				w.writeEmptyElement("use");
				w.writeAttribute("class","kgstrand");
				w.writeAttribute("xlink", XLINK, "href", "#strand"+(g.isPositiveStrand()?"F":"R"));
				w.writeAttribute("x",String.valueOf(pixX));
				w.writeAttribute("y",String.valueOf(midY));
				}
			
			/* exons */
			for(KnownGene.Exon exon:g.getExons())
					{
					if(exon.getStart()>= this.interval.end) continue;
					if(exon.getEnd()<= this.interval.start) continue;
					w.writeEmptyElement("rect");
					w.writeAttribute("class","kgexon");
					w.writeAttribute("title",exon.getName());
					w.writeAttribute("x",String.valueOf(baseToPixel(trim(exon.getStart()))));
					w.writeAttribute("y",String.valueOf(midY-exonHeight/2));
					w.writeAttribute("width",String.valueOf(baseToPixel(trim(exon.getEnd()))-baseToPixel(trim(exon.getStart()))));
					w.writeAttribute("height",String.valueOf(exonHeight));
					}
			
			/* coding line */
			w.writeEmptyElement("rect");
			w.writeAttribute("class","kgcds");
			w.writeAttribute("x",String.valueOf(baseToPixel(trim(g.getCdsStart()))));
			w.writeAttribute("y",String.valueOf(midY-cdsHeigh/2));
			w.writeAttribute("width",String.valueOf(baseToPixel(trim(g.getCdsEnd()))-baseToPixel(trim(g.getCdsStart()))));
			w.writeAttribute("height",String.valueOf(cdsHeigh));

			
			String label=String.format("%15s", g.getName());
			w.writeEmptyElement("path");
			double fontHeight=Math.min(10,0.8*TRANSCRIPT_HEIGHT);
			w.writeAttribute("d",this.hershey.svgPath(label,-insets.left,midY-fontHeight/2,insets.left*0.9,fontHeight));

			
			w.writeEndElement();
			y+=TRANSCRIPT_HEIGHT;
			}
		}
	
	private void writeGCPercentSection( XMLStreamWriter w ) throws XMLStreamException
		{
		if(genomicSequence==null || !this.interval.getChrom().equals(genomicSequence.getChrom()))
			{
			this.genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, this.interval.getChrom());
			}		
		
		final int drawinAreaWidth= screenWidth-(this.insets.right+this.insets.left);
		
		/* get GC% */
		float gcPercent[]=new float[drawinAreaWidth];
		for(int x=0;x< gcPercent.length;++x )
			{
			int n=0;
			int countGC=0;
			int pos0_a= this.interval.getStart()+(int)(((x+0)/((double)drawinAreaWidth))*this.interval.distance());
			int pos0_b= this.interval.getStart()+(int)(((x+1)/((double)drawinAreaWidth))*this.interval.distance());
			int win=gcWindowSize;
			if(pos0_b-pos0_a > gcWindowSize)win=Math.max(0, (pos0_b-pos0_a)-gcWindowSize);
			for(int pos0=Math.max(0,pos0_a-win);
					pos0<Math.min(genomicSequence.length(),pos0_b+win);
					++pos0)
				{
				switch(Character.toUpperCase(genomicSequence.charAt(pos0)))
					{
					case 'C':case 'G': case 'S': countGC++; break;
					default:break;
					}
				n++;
				}
			gcPercent[x]=((float)countGC)/((float)n);
			}
		List<Point2D.Double> points=new ArrayList<Point2D.Double>(drawinAreaWidth);
		points.add(new Point2D.Double(0,GC_PERCENT_HEIGHT));//left,bottom
		
		
		int pixelX=0;
		while(pixelX< gcPercent.length)
			{
			int x2=pixelX+1;
			while(x2< gcPercent.length && (float)gcPercent[pixelX]==(float)gcPercent[x2])
				{
				++x2;
				}
			double y2=GC_PERCENT_HEIGHT-(float)GC_PERCENT_HEIGHT*gcPercent[pixelX];
			points.add(new Point2D.Double(pixelX,y2));
			points.add(new Point2D.Double(x2,y2));
			pixelX=x2;
			}
		points.add(new Point2D.Double(drawinAreaWidth,GC_PERCENT_HEIGHT));//right,bottom
		points.add(new Point2D.Double(0,GC_PERCENT_HEIGHT));//left,bottom
		StringBuilder sw=new StringBuilder();
		for(Point2D.Double pt:points) sw.append((int)pt.getX()).append(",").append((float)pt.getY()).append(" ");
		w.writeEmptyElement("polygon");
		w.writeAttribute("points",sw.toString().trim());
		w.writeAttribute("class","gcpercent");
		
		w.writeEmptyElement("rect");
		w.writeAttribute("class", "frame");
		w.writeAttribute("x",String.valueOf(0));
		w.writeAttribute("y","0");
		w.writeAttribute("width",String.valueOf(drawinAreaWidth));
		w.writeAttribute("height",String.valueOf(GC_PERCENT_HEIGHT));
		}
	
	private void paint(
			XMLStreamWriter w
			) throws XMLStreamException,IOException
		{
		final int font_size=10;
		List<KnownGene> operon=new ArrayList<KnownGene>();
		
		if(this.knownGenes!=null)
			{
			Pattern tab=Pattern.compile("[\t]");
			for(Iterator<String> iter=this.knownGenes.iterator(
					this.interval.chrom, this.interval.start,this.interval.end);
					iter.hasNext();
					)
				{	
				String line=iter.next();
				operon.add(new KnownGene(tab.split(line)));
				}
			}
		
		
	
		final Dimension svgSize=new Dimension();
		final int drawinAreaWidth= screenWidth-(this.insets.right+this.insets.left);
		svgSize.width=screenWidth;
		svgSize.height=1500;

		
		
		List<VariantContext> variants=new ArrayList<VariantContext>();
		if(this.tabixReader!=null)
			{
			for(Iterator<VariantContext> iter=tabixReader.iterator(interval.chrom, this.interval.getStart(), this.interval.getEnd());
					iter.hasNext();)
				{
				VariantContext ctx=iter.next();
				variants.add(ctx);
				}
			}
		
		svgSize.height=this.insets.top+this.insets.bottom+
					operon.size()*TRANSCRIPT_HEIGHT+
					samFileReaders.size()*BAM_COVERAGE_HEIGHT+
					GC_PERCENT_HEIGHT+
					variants.size()
					;
		svgSize.height=1550;
		
		w.writeStartElement("svg");
		w.writeDefaultNamespace(SVG.NS);
		w.writeNamespace("xlink", XLINK);
		w.writeAttribute("version", "1.1");
		w.writeAttribute("width",String.valueOf(svgSize.width));
		w.writeAttribute("height",String.valueOf(svgSize.height));
		
		
		w.writeStartElement("title");
		w.writeCharacters(this.interval.chrom+":"+this.interval.getStart()+"-"+this.interval.getEnd());
		w.writeEndElement();
		
		w.writeStartElement("defs");

		//genotypes
		w.writeEmptyElement("rect");
		w.writeAttribute("id","g_"+GenotypeType.HOM_REF); //
		w.writeAttribute("style","fill:none;stroke;black;");
		w.writeAttribute("x", "-15" );
		w.writeAttribute("y", "-15" );
		w.writeAttribute("width", "30" );
		w.writeAttribute("height", "30" );
		
		w.writeEmptyElement("rect");
		w.writeAttribute("id","g_"+GenotypeType.HOM_VAR); //
		w.writeAttribute("style","fill:black;stroke;black;");
		w.writeAttribute("x", "-15" );
		w.writeAttribute("y", "-15" );
		w.writeAttribute("width", "30" );
		w.writeAttribute("height", "30" );
		
		w.writeStartElement("g");
		w.writeAttribute("id","g_"+GenotypeType.HET); //
			w.writeEmptyElement("rect");
			w.writeAttribute("style","fill:none;stroke;black;");
			w.writeAttribute("x", "-15" );
			w.writeAttribute("y", "-15" );
			w.writeAttribute("width", "30" );
			w.writeAttribute("height", "30" );
			w.writeEmptyElement("polygon");
			w.writeAttribute("style","fill:black;stroke;black;");
			w.writeAttribute("points","-15,-15 15,-15 15,15 -15,-15");
		w.writeEndElement();

		

		
		//strand
		w.writeEmptyElement("polyline");
			w.writeAttribute("id","strandF");
			w.writeAttribute("points", "-5,-5 0,0 -5,5" );
		
		w.writeEmptyElement("polyline");
			w.writeAttribute("id","strandR");
			w.writeAttribute("points", "5,-5 0,0 5,5" );
		
		//gradients
			w.writeStartElement("linearGradient");
				w.writeAttribute("id","grad01");
				w.writeAttribute("x1","50%");
				w.writeAttribute("x2","50%");
				w.writeAttribute("y1","0%");
				w.writeAttribute("y2","100%");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","0%");
					w.writeAttribute("style","stop-color:black;stop-opacity:1;");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","50%");
					w.writeAttribute("style","stop-color:white;stop-opacity:1;");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","100%");
					w.writeAttribute("style","stop-color:black;stop-opacity:1;");
			w.writeEndElement();

			w.writeStartElement("linearGradient");
				w.writeAttribute("id","grad02");
				w.writeAttribute("x1","50%");
				w.writeAttribute("x2","50%");
				w.writeAttribute("y1","0%");
				w.writeAttribute("y2","100%");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","0%");
					w.writeAttribute("style","stop-color:steelblue;stop-opacity:1;");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","100%");
					w.writeAttribute("style","stop-color:lightblue;stop-opacity:1;");
			w.writeEndElement();

			w.writeStartElement("linearGradient");
				w.writeAttribute("id","grad03");
				w.writeAttribute("x1","50%");
				w.writeAttribute("x2","50%");
				w.writeAttribute("y1","100%");
				w.writeAttribute("y2","0%");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","0%");
					w.writeAttribute("style","stop-color:red;stop-opacity:1;");
				w.writeEmptyElement("stop");
					w.writeAttribute("offset","20%");
					w.writeAttribute("style","stop-color:lightblue;stop-opacity:1;");
			w.writeEndElement();
		
			
		if(indexedFastaSequenceFile!=null)
			{
			for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
				{
				double width=drawinAreaWidth/(double)this.interval.distance();
				w.writeEmptyElement("path");
				w.writeAttribute("id","base"+base);
				w.writeAttribute("title",base);
				w.writeAttribute("class","base"+base.toUpperCase());
				w.writeAttribute("d",this.hershey.svgPath(
						base,
						0,
						0,
						width*0.95,
						font_size
						));
				}
			}
		
		w.writeEndElement();//defs

		w.writeStartElement("style");
		w.writeCharacters(
				"svg {fill:none; stroke:black;}\n"+
				".ruler-label { stroke:red;}\n"+
				".frame { stroke:black;fill:none;}\n"+
				".kgexon {fill:url(#grad01);stroke:black;}\n"+
				".gcpercent {fill:url(#grad02);stroke:black;}"+
				".coverage {fill:url(#grad03);stroke:black;}"+
				".kgcds {fill:mediumpurple;stroke:black;}\n"+
				".variant{stroke:none;fill:red;opacity:0.2;}\n"+
				".xaxis{stroke:gray;fill:none;opacity:0.2;}"
				);
		w.writeEndElement();//style
		
		w.writeStartElement("g");
		
		int y=insets.top;
		
		
	
		
		/* left and right position */
		{
		w.writeStartElement("g");
		w.writeAttribute("title","ruler");
		w.writeAttribute("transform","translate("+insets.left+","+y+")");
		w.writeAttribute("class", "ruler-label");


		
		NumberFormat fmt = DecimalFormat.getInstance();
		String label=fmt.format(this.interval.getStart());
		
		w.writeEmptyElement("path");
		w.writeAttribute("title",label);
		w.writeAttribute("d",this.hershey.svgPath(label, 0, 0, label.length()*font_size, font_size));
		
		label=fmt.format(this.interval.getEnd());
		w.writeEmptyElement("path");
		w.writeAttribute("title",label);
		w.writeAttribute("d",this.hershey.svgPath(label, drawinAreaWidth - label.length()*font_size, 0,label.length()*font_size, font_size));
		w.writeEndElement();
		
		y+= font_size;
		}
		
		/* print bases / sequences */
		if(this.indexedFastaSequenceFile!=null &&
			font_size*this.interval.distance()< drawinAreaWidth)
			{
			LOG.info("paint DNA sequence");
			if(genomicSequence==null || !this.interval.getChrom().equals(genomicSequence.getChrom()))
				{
				this.genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, this.interval.getChrom());
				}
			w.writeStartElement("g");
			w.writeAttribute("title","sequence");
			w.writeAttribute("transform","translate("+insets.left+","+y+")");
			w.writeAttribute("class", "sequence");
			
			for(int i=0;
					i< this.interval.distance() && (i+this.interval.start)< genomicSequence.length() ;
					++i)
				{
				String base=String.valueOf(this.genomicSequence.charAt(i+this.interval.start));
				double width=drawinAreaWidth/(double)this.interval.distance();
				w.writeEmptyElement("use");
				w.writeAttribute("y","0");
				w.writeAttribute("x",String.valueOf(i*width));
				w.writeAttribute("title",base+"("+String.valueOf(i+this.interval.start+1)+")");
				w.writeAttribute("xlink",XLINK,"href","#base"+base);
				
				}
			w.writeEndElement();
			
			y+= font_size;
			}
		else
			{
			LOG.info("won't display bases."+font_size*this.interval.distance()+" "+drawinAreaWidth);
			}
		
		if(this.indexedFastaSequenceFile!=null)
			{
			/* GC % */
			w.writeStartElement("g");
			w.writeAttribute("title","gc%");
			w.writeAttribute("transform","translate("+insets.left+","+y+")");
			writeGCPercentSection(w);
			w.writeEndElement();
			y+=GC_PERCENT_HEIGHT;
			}
		
		
			{
			w.writeStartElement("g");
			w.writeAttribute("title","kg");
			w.writeAttribute("transform","translate("+insets.left+","+y+")");
	
			writeKownGeneSection(w, operon);
			w.writeEndElement();
			y+=operon.size()*TRANSCRIPT_HEIGHT;
			}
			
			{
			w.writeStartElement("g");
			w.writeAttribute("title","bams");
			w.writeAttribute("transform","translate("+insets.left+","+y+")");
	
			writeBamSection(w);
			w.writeEndElement();
			y+=operon.size()*BAM_COVERAGE_HEIGHT;
			}	
		
			
			
			/* variants */
			{
			w.writeStartElement("g");
			w.writeAttribute("transform","translate("+insets.left+",0)");

			w.writeAttribute("title","variants");
			
			
			for(VariantContext ctx:variants)
				{
				if(ctx.getStart()< (this.interval.start+1)) continue;
				if(ctx.getStart()>= (this.interval.end)) continue;
				
				w.writeStartElement("g");
				w.writeAttribute("title", "("+ctx.getStart()+")");
				
				Rectangle2D.Double rect=new Rectangle2D.Double();
				rect.x= ((ctx.getStart()-(interval.start+1))/(double)interval.distance())*drawinAreaWidth;
				rect.width= (((ctx.getEnd()+1)-(interval.start+1))/(double)interval.distance())*drawinAreaWidth - rect.x;
				rect.y=0;
				rect.height=1000;
				
				w.writeEmptyElement("rect");
				w.writeAttribute("class","variant");
				w.writeAttribute("title", "("+ctx.getStart()+")");
				w.writeAttribute("x", String.valueOf(rect.x));
				w.writeAttribute("y", String.valueOf(rect.y));
				w.writeAttribute("width", String.valueOf(rect.width));
				w.writeAttribute("height", String.valueOf(rect.height));
				
				w.writeEndElement();
				
				}
			
			
			w.writeEndElement();
			}
			
		/* samples */
		if(tabixReader!=null)
			{
			List<String> samples=tabixReader.getHeader().getSampleNamesInOrder();
			w.writeStartElement("g");
			w.writeAttribute("transform","translate("+insets.left+","+y+")");

			w.writeAttribute("title","samples");

			for(int i=0;i< samples.size();++i)
				{
				w.writeStartElement("g");
				w.writeAttribute("transform","translate(0,"+(i*SAMPLE_HEIGHT)+")");

				w.writeAttribute("title",samples.get(i));

				
				w.writeEmptyElement("line");
				w.writeAttribute("class","sample");
				w.writeAttribute("style","stroke:black;");

				w.writeAttribute("x1",String.valueOf(0));
				w.writeAttribute("y1",String.valueOf(SAMPLE_HEIGHT/2));
				w.writeAttribute("x2",String.valueOf(drawinAreaWidth));
				w.writeAttribute("y2",String.valueOf(SAMPLE_HEIGHT/2));
				
				String label=String.format("%10s", samples.get(i));
				w.writeEmptyElement("path");
				w.writeAttribute("title",label);
				w.writeAttribute("style","stroke:black;");
				w.writeAttribute("d",this.hershey.svgPath(
						label,
						-label.length()*font_size,
						SAMPLE_HEIGHT/2-font_size/2,
						label.length()*font_size,
						font_size)
						);
				
				
				
				/* variants */
				for(VariantContext ctx: variants)
					{

					for(String sample:tabixReader.getHeader().getSampleNamesInOrder() )
						{
						Genotype genotype=ctx.getGenotype(sample);
						if(genotype==null || !genotype.isCalled() || !genotype.isAvailable())
							{
							continue;
							}
						double pixX= ((((ctx.getStart()+(ctx.getEnd()+1))/2.0)-(interval.start+1))/(double)interval.distance())*drawinAreaWidth;
						String gType=null;
						switch(genotype.getType())
							{
							case HET: case HOM_REF: case HOM_VAR:
								{
								gType="g_"+genotype.getType().name();
								break;
								}
							default: LOG.error("Cannot handle "+genotype.getType());break;
							}
						
						if(gType==null) continue;
						w.writeEmptyElement("use");
						w.writeAttribute("xlink", XLINK, "href", "#"+gType);
						w.writeAttribute("x",String.valueOf(pixX));
						w.writeAttribute("y",String.valueOf(SAMPLE_HEIGHT/2));

						}
					}

				
				}
			y+=SAMPLE_HEIGHT*samples.size();
			}
		
		
		
		
		w.writeEndElement();//g
		w.writeEndElement();//svg
		}
	
	private int trim(int pos0)
		{
		return Math.min(Math.max(pos0, interval.start),interval.end);
		}
	
	private Interval parseIntervalString(String input)
		{
		Interval R=new Interval();
		R.chrom=null;
		R.start=1;
		R.end=Integer.MAX_VALUE;

		input=input.replace(",","");
		int colon=input.indexOf(':');
		if(colon!=-1)
			{
			R.chrom=input.substring(0,colon);
			String rm=input.substring(colon+1);
			int dash=rm.indexOf('-');
			if(dash==-1)
				{
				R.start=Integer.parseInt(rm);
				}
			else
				{
				R.start=Integer.parseInt(rm.substring(0,dash));
				R.end=Integer.parseInt(rm.substring(dash+1));
				}
			}
		return	R.chrom!=null &&
				R.start<= R.end &&
				R.start>=0 ? R: null;
		}
	
	
	private void run(Interval R)throws IOException,XMLStreamException
		{
		if(R==null || R.end<R.start) return;
		this.interval=R;
		LOG.info("Processing "+this.interval);
		XMLOutputFactory xof=XMLOutputFactory.newFactory();
		xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.TRUE);

	
		OutputStream stream=null;
		XMLStreamWriter w=null;
		if(fileOutPattern!=null)
			{
			File fout=new File(this.fileOutPattern.replaceAll(REPLACE, String.format("%s:%09d-%09d",this.interval.chrom,this.interval.start,this.interval.end)));
			LOG.info("Wrting to "+fout);
			if(fout.getParentFile()!=null)
				{
				fout.getParentFile().mkdirs();
				}
			
			stream=new FileOutputStream(fout);
			if(fout.getName().toLowerCase().endsWith(".gz"))
				{
				stream=new GZIPOutputStream(stream);
				}
			w=xof.createXMLStreamWriter(stream, "UTF-8");
			}
		else
			{
			w=xof.createXMLStreamWriter(System.out, "UTF-8");
			}
		w.writeStartDocument("UTF-8", "1.0");
		paint(w);
		w.writeEndDocument();
		w.flush();
		w.close();
		if(stream!=null)
			{
			stream.flush();
			stream.close();
			}
		}
	
	private void run(LineIterator lr)
		throws IOException,XMLStreamException
		{
		Pattern tab=Pattern.compile("[\t]");
		while(lr.hasNext())
			{
			String line=lr.next();
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line);
			this.interval.chrom=null;
			this.interval.start=1;
			this.interval.end=Integer.MAX_VALUE;

			if(tokens.length>2)
				{
				this.interval.chrom=tokens[0];
				this.interval.start=Integer.parseInt(tokens[1]);
				this.interval.end=Integer.parseInt(tokens[2]);
				}
			else
				{
				this.interval=parseIntervalString(tokens[0]);
				if(this.interval==null)
					{
					LOG.error("info bad interval "+tokens[0]);
					continue;
					}
				}
			run(this.interval);
			}
		}
	
	
	
	@Parameter(names="-r",description="(region chr:start-end")
	private String intervalStr=null;
	@Parameter(names="-k",description="knownGene file")
	private String knownGeneUri=null;
	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File faidx=null;
	@Parameter(names="-B",description="bam files")
	private List<File> bamList=new ArrayList<>();
	@Parameter(names="-V",description="VCF file")
	private String vcfFile=null;
	
	@Override
	public int doWork(List<String> args) {
		Interval userInterval=null;
		String vcfFile=null;
		
		
		
		if(knownGeneUri==null)
			{
			LOG.warning("KnownGene URI undefined.");
			}
		if(faidx==null)
			{
			LOG.warning("Undefined fasta Reference.");
			}
		
		for(File b:bamList)
			{
			BamFile bf=new BamFile();
			bf.file=(b);
			this.samFileReaders.add(bf);
			}
		
		if(intervalStr!=null)
			{
			if((userInterval=parseIntervalString(intervalStr))==null)
				{
				LOG.error("Bad interval. "+intervalStr);
				return -1;
				}
			}
		
		try
			{
			SAMSequenceDictionary dict=null;
			
			if(faidx!=null)
				{
				LOG.info("Opening "+faidx);
				this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
				dict=this.indexedFastaSequenceFile.getSequenceDictionary();
				}
			
			final SamReaderFactory srf= SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT); 
			for(BamFile input: samFileReaders)
				{
				LOG.info("Opening "+input.file);
				input.sfr=srf.open(input.file);
				if(dict!=null && !SequenceUtil.areSequenceDictionariesEqual(
						input.sfr.getFileHeader().getSequenceDictionary(),
						dict))
					{
					input.sfr.close();
					LOG.error("NOT the same sequence dictionaries between "+faidx+" and "+input.file);
					return -1;
					}
				}
				
			if(vcfFile!=null)
				{
				LOG.info("Opening "+vcfFile);
				this.tabixReader=new TabixVcfFileReader(vcfFile);
				
				}
			
			if(knownGeneUri!=null)
				{
				LOG.info("Opening "+knownGeneUri);
				this.knownGenes=new TabixFileReader(knownGeneUri);
				}
			
			if(args.isEmpty())
				{
				if(userInterval!=null)
					{
					run(userInterval);
					return 0;
					}
				
				LOG.info("Reading chr:start-end or BED from stdin");
				LineIterator lr=IOUtils.openStdinForLineIterator();
				run(lr);
				CloserUtil.close(lr);
				}
			else
				{
				for(String filename:args)
					{
					if(userInterval!=null)
						{
						LOG.error("Illegal number of arguments (user interval specified).");
						return -1;
						}
					
					LOG.info("Reading chr:start-end or BED from "+filename);
					LineIterator lr=IOUtils.openURIForLineIterator(filename);
					run(lr);
					CloserUtil.close(lr);
					}
				}
			
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.samFileReaders);
			}
		}
	
	public static void main(String[] args)
		{
		new PaintContext().instanceMainWithExit(args);
		}
	
}
