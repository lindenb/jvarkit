package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.Dimension;
import java.awt.Insets;
import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.broad.tribble.readers.LineIterator;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Hershey;
import com.github.lindenb.jvarkit.util.TabixFileReader;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;

public class PaintContext extends AbstractCommandLineProgram
	{
	private static final String XLINK="http://www.w3.org/1999/xlink";
	private static final String REPLACE="__REGION__";
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;
	private List<SAMFileReader> samFileReaders=new ArrayList<SAMFileReader>();
	private TabixVcfFileReader tabixReader=null;
	private TabixFileReader knownGenes=null;
	private Insets insets=new Insets(100, 100, 100, 100);
	private int screenWidth=2000;
	private int gcWindowSize=5;
	private String fileOutPattern=null;
	private Hershey hershey=new Hershey();
	
	private abstract class ROD
		{
		String name;
		File file;
		}
	
	private class BamROD
		extends ROD
		{
		SAMFileReader sfr;
		}
	
	private class VCFRod
		extends ROD
		{
		TabixVcfFileReader tabixReader;
		}
	
	private class KnownGenesRod
		extends ROD
		{
		private Map<String,List<KnownGene>> chrom2genes=new HashMap<>();
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
		}
	
	
	/* current interval */
	private Interval interval;
	
	
	
	
	private PaintContext()
		{
		
		}
	
	
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	
	
	private double baseToPixel(int pos0)
		{
		return  ((pos0 - this.interval.getStart())/(double)this.interval.distance())*(this.screenWidth-(this.insets.right+this.insets.left) )
				;
		}
	
	private final int TRANSCRIPT_HEIGHT=30;
	private final int GC_PERCENT_HEIGHT=100;
	private final int BAM_COVERAGE_HEIGHT=100;
	private int maxDepth=100;
	
	private void writeBamSection(
			XMLStreamWriter w,
			Interval minMaxBase
			) throws XMLStreamException
		{
		final int drawinAreaWidth= screenWidth-(this.insets.right+this.insets.left);
		int y=0;
		for(SAMFileReader sfr:this.samFileReaders)
			{
			w.writeStartElement("g");
			w.writeAttribute("transform", "translate(0,"+y+")");
			int depth[]=new int[drawinAreaWidth];
			
			SAMRecordIterator iter=null;
			try
				{
				iter=sfr.query(
						genomicSequence.getChrom(),
						minMaxBase.getStart()+1,
						minMaxBase.getEnd(),
						true
						);
				while(iter.hasNext())
					{
					SAMRecord rec=iter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getNotPrimaryAlignmentFlag()) continue;
					if(rec.getDuplicateReadFlag()) continue;
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					for(int i1=Math.max(minMaxBase.getStart()+1, rec.getAlignmentStart());
							(i1-(minMaxBase.getStart()+1)) < depth.length && i1 <= rec.getAlignmentEnd();
							++i1)
						{
						int index = i1-(minMaxBase.getStart()+1);
						depth[index]=Math.max(depth[index]+1,maxDepth);
						}
					}
				}
			catch(Exception err)
				{
				error(err);
				throw new RuntimeException(err);
				}
			finally
				{
				CloserUtil.close(iter);
				}
			
			
			w.writeAttribute("title","BAM");
			List<Point2D.Double> points=new ArrayList<Point2D.Double>();
			points.add(new Point2D.Double(insets.left,BAM_COVERAGE_HEIGHT));//left,bottom
			
			
			int pixelX=0;
			while(pixelX< depth.length)
				{
				int x2=pixelX+1;
				while(x2< depth.length && depth[pixelX]==depth[x2])
					{
					++x2;
					}
				double y2=BAM_COVERAGE_HEIGHT-(float)BAM_COVERAGE_HEIGHT*(depth[pixelX]/(double)maxDepth);
				points.add(new Point2D.Double(pixelX,y2));
				points.add(new Point2D.Double(x2,y2));
				pixelX=x2;
				}
			points.add(new Point2D.Double(drawinAreaWidth,BAM_COVERAGE_HEIGHT));//right,bottom
			StringBuilder sw=new StringBuilder();
			for(Point2D.Double pt:points) sw.append(pt.getX()).append(",").append(pt.getY()).append(" ");
			w.writeEmptyElement("polyline");
			w.writeAttribute("points",sw.toString().trim());
			
			//frame
			w.writeEmptyElement("rec");
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
		w.writeStartElement("defs");
		w.writeStartElement("mask");
		w.writeAttribute("id", "kgclip");
		w.writeEmptyElement("rec");
		w.writeAttribute("x", "0");
		w.writeAttribute("y", "0");
		w.writeAttribute("width",String.valueOf(drawinAreaWidth));
		w.writeAttribute("height",String.valueOf(operon.size()*TRANSCRIPT_HEIGHT));
		w.writeEndElement();//masks
		w.writeEndElement();//defs

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
			w.writeAttribute("x1",String.valueOf(baseToPixel(g.getTxStart())));
			w.writeAttribute("y1",String.valueOf(midY));
			w.writeAttribute("x2",String.valueOf(baseToPixel(g.getTxEnd())));
			w.writeAttribute("y2",String.valueOf(midY));
			
			
			/* coding line */
			w.writeEmptyElement("rect");
			w.writeAttribute("class","kgcds");
			w.writeAttribute("x",String.valueOf(baseToPixel(g.getCdsStart())));
			w.writeAttribute("y",String.valueOf(midY-cdsHeigh/2));
			w.writeAttribute("width",String.valueOf(baseToPixel(g.getCdsEnd())-baseToPixel(g.getCdsStart())));
			w.writeAttribute("height",String.valueOf(cdsHeigh));
			
			
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
					w.writeEmptyElement("rect");
					w.writeAttribute("class","kgexon");
					w.writeAttribute("title",exon.getName());
					w.writeAttribute("x",String.valueOf(baseToPixel(exon.getStart())));
					w.writeAttribute("y",String.valueOf(midY-exonHeight/2));
					w.writeAttribute("width",String.valueOf(baseToPixel(exon.getEnd())-baseToPixel(exon.getStart())));
					w.writeAttribute("height",String.valueOf(exonHeight));
					}
			
			
			
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
			while(x2< gcPercent.length && gcPercent[pixelX]==gcPercent[x2])
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
		for(Point2D.Double pt:points) sw.append(pt.getX()).append(",").append(pt.getY()).append(" ");
		w.writeEmptyElement("polygon");
		w.writeAttribute("points",sw.toString().trim());
		
		w.writeEmptyElement("rec");
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
		svgSize.height=1000;

		
		
		List<VariantContext> variants=new ArrayList<>();
		if(this.tabixReader!=null)
			{
			for(Iterator<VariantContext> iter=tabixReader.iterator(genomicSequence.getChrom(), this.interval.getStart(), this.interval.getEnd());
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
		w.writeEmptyElement("polyline");
			w.writeAttribute("id","strandF");
			w.writeAttribute("points", "-5,-5 0,0 -5,5" );
		
		w.writeEmptyElement("polyline");
			w.writeAttribute("id","strandR");
			w.writeAttribute("points", "5,-5 0,0 5,5" );
		
		
		
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
				".frame { stroke:gray;}"
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
			info("paint DNA sequence");
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
			info("won't display bases."+font_size*this.interval.distance()+" "+drawinAreaWidth);
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
		}
		
		/* samples */
		if(tabixReader!=null)
			{
			for(String sample:tabixReader.getHeader().getSampleNamesInOrder() )
				{
				w.writeEmptyElement("line");
				w.writeAttribute("class","variant");
				w.writeAttribute("x1",String.valueOf(0));
				w.writeAttribute("y1","0");
				w.writeAttribute("x2",String.valueOf(1000));
				w.writeAttribute("y1",String.valueOf(1000));
				}
			}
		/* variants */
		for(VariantContext var: variants)
			{
			double pos0=baseToPixel(var.getStart());
			w.writeEmptyElement("line");
			w.writeAttribute("class","variant");
			w.writeAttribute("x1",String.valueOf(pos0));
			w.writeAttribute("y1","0");
			w.writeAttribute("x2",String.valueOf(pos0));
			w.writeAttribute("y1",String.valueOf(1000));
			for(String sample:tabixReader.getHeader().getSampleNamesInOrder() )
				{
				Genotype genotype=var.getGenotype(sample);
				if(genotype==null || !genotype.isCalled() || !genotype.isAvailable())
					{
					continue;
					}
				
				}
			}
		
		
		w.writeEndElement();//g
		w.writeEndElement();//svg
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
		if(R==null) return;
		this.interval=R;
		info("Processing "+this.interval);
		XMLOutputFactory xof=XMLOutputFactory.newFactory();
		xof.setProperty(XMLOutputFactory.IS_REPAIRING_NAMESPACES, Boolean.TRUE);

	
		OutputStream stream=null;
		XMLStreamWriter w=null;
		if(fileOutPattern!=null)
			{
			File fout=new File(this.fileOutPattern.replaceAll(REPLACE, String.format("%s:%09d-%09d",this.interval.chrom,this.interval.start,this.interval.end)));
			info("Wrting to "+fout);
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
					error("info bad interval "+tokens[0]);
					continue;
					}
				}
			run(this.interval);
			}
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -r (region chr:start-end) . Optional");
		out.println(" -R (file) reference fasta file indexed with faidx . Optional");
		out.println(" -k (file) tabix indexed knownGene list . Optional");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		Interval userInterval=null;
		String vcfFile=null;
		List<File> bamFiles=new ArrayList<File>();
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		File faidx=null;
		String knownGeneUri=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"k:R:B:V:o:r:"))!=-1)
			{
			switch(c)
				{			
				case 'r': 
					{
					if((userInterval=parseIntervalString(opt.getOptArg()))==null)
						{
						error("Bad interval. "+opt.getOptArg());
						return -1;
						}
					break;
					}
				case 'o': fileOutPattern=opt.getOptArg();break;
				case 'B': bamFiles.add(new File(opt.getOptArg()));break;
				case 'R': faidx=new File(opt.getOptArg());break;
				case 'k': knownGeneUri=opt.getOptArg();break;
				case 'V': vcfFile=opt.getOptArg();break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			
			}
		
		if(knownGeneUri==null)
			{
			warning("KnownGene URI undefined.");
			}
		if(faidx==null)
			{
			warning("Undefined fasta Reference.");
			}
		
		
		
		try
			{
			SAMSequenceDictionary dict=null;
			
			if(faidx!=null)
				{
				info("Opening "+faidx);
				this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
				dict=this.indexedFastaSequenceFile.getSequenceDictionary();
				}
			
			
			for(File bamFile:bamFiles)
				{
				info("Opening "+bamFile);
				SAMFileReader samFileReader=new SAMFileReader(faidx);
				if(dict!=null && !SequenceUtil.areSequenceDictionariesEqual(
						samFileReader.getFileHeader().getSequenceDictionary(),
						dict))
					{
					samFileReader.close();
					error("NOT the same sequence dictionaries between "+faidx+" and "+bamFile);
					return -1;
					}
				this.samFileReaders.add(samFileReader);
				}
				
			if(vcfFile!=null)
				{
				info("Opening "+vcfFile);
				this.tabixReader=new TabixVcfFileReader(vcfFile);
				
				}
			
			if(knownGeneUri!=null)
				{
				info("Opening "+knownGeneUri);
				this.knownGenes=new TabixFileReader(knownGeneUri);
				}
			
			if(opt.getOptInd()==args.length)
				{
				if(userInterval!=null)
					{
					run(userInterval);
					return 0;
					}
				
				info("Reading chr:start-end or BED from stdin");
				LineIterator lr=IOUtils.openStdinForLineIterator();
				run(lr);
				CloserUtil.close(lr);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					if(userInterval!=null)
						{
						error("Illegal number of arguments (user interval specified).");
						return -1;
						}
					
					String filename=args[i];
					info("Reading chr:start-end or BED from "+filename);
					LineIterator lr=IOUtils.openURIForLineIterator(filename);
					run(lr);
					CloserUtil.close(lr);
					}
				}
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
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
