package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Dimension;
import java.awt.Insets;
import java.io.BufferedReader;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;

import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC


END_DOC
*/
@Program(name="vcf2postscript",
	description="Print VCF context as Postscript",
		keywords={"vcf","postscript"},
		modificationDate="20191018"
		)
public class VcfToPostscript extends Launcher
	{
	private final static Logger LOG=Logger.of(VcfToPostscript.class);
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;

	private final List<Transcript> genes=new ArrayList<>();
	private final Set<Integer> positions=new HashSet<Integer>();
	private final Map<String,Set<Integer>> sample2positions=new TreeMap<>();
	private int chromStart=Integer.MAX_VALUE;
	private int chromEnd=Integer.MIN_VALUE;
	
	@Parameter(names={"-g","--gtf"},description=GtfReader.OPT_DESC)
	private Path gtfPath = null;


	private final Map<String,List<Transcript>> chrom2transcript=new HashMap<>();
	private Insets margin=new Insets(100, 200, 50, 50);
	private Dimension pageDef=new Dimension(900,700);
	private float fHeight=20;
	private int count_pages_printed=0;

	private VcfToPostscript()
		{
		}

		
	private void clear()
		{
		genes.clear();
		positions.clear();
		sample2positions.clear();
		chromStart=Integer.MAX_VALUE;
		chromEnd=Integer.MIN_VALUE;
		}
	
	private void addGene(final Transcript g)
		{
		chromStart=Math.min(chromStart, g.getStart());
		chromEnd=Math.max(chromEnd, g.getEnd());
		genes.add(g);
		}
	
	private double toPixel(final int pos)
		{
		return (double)margin.left+((pos-(double)this.chromStart)/((double)this.chromEnd-(double)this.chromStart))*(double)pageDef.width;
		}

			
	private void addVariant(final VariantContext ctx)
		{
		if(!ctx.getContig().equals(genes.get(0).getContig())) return;
		if(ctx.getStart()>chromEnd) return;
		if(ctx.getEnd()<chromStart) return;
		positions.add(ctx.getStart());
		for(final String sample: ctx.getSampleNames())
			{
			final Genotype g=ctx.getGenotype(sample);
			if(!g.isAvailable()) continue;
			if(!g.isCalled()) continue;
			if(g.isNoCall()) continue;
			if(g.isNonInformative()) continue;
			Set<Integer> set=sample2positions.get(sample);
			if(set==null)
				{
				set=new HashSet<Integer>();
				sample2positions.put(sample, set);
				}
			set.add(ctx.getStart());
			}
		}
			
	private void print()
		{
		if(this.positions.isEmpty())
			{
			return;
			}
		++count_pages_printed;
		
		final Transcript first=genes.get(0);
		final double fHeight=20;
		
		
		final Dimension localPage=new Dimension(
				(int)(this.pageDef.width),
				(int)(margin.top+margin.bottom + (this.genes.size()+1)*fHeight + this.sample2positions.size()*fHeight)
				);
		
		this.outw.println(  "\n%%Page: "  +  count_pages_printed  +  " " + count_pages_printed);
		this.outw.println("gsave");
		this.outw.println(
				1.0 +" "+ 
				(localPage.height <= this.pageDef.height ? 1.0 : 1.0/(localPage.getHeight()/(float)this.pageDef.getHeight()))
				+ " scale");

		//this.outw.println( "%%BoundingBox: 0 0 "  + page.width +  " "  + page.height  );
		

		float midy=(float)(fHeight/2.0f);

		double cdsHeight=fHeight*0.4;
		double exonHeight=fHeight*0.9;

		this.outw.println(
				"2 "  +  (localPage.height-10 ) +
				" moveto ("  + first.getContig() + ":"+this.chromStart+ "-" +this.chromEnd+") show"
				 )
				;


		this.outw.println("1 0 0 setrgbcolor");
		this.outw.println("0.3 setlinewidth");
		for(Integer r: this.positions )
			{
			this.outw.print(  "newpath " +  (float)toPixel(r)  +  " 0 moveto 0 " +  localPage.height  +  " rlineto stroke\n");
			this.outw.print(  (float)toPixel(r)  +  " "  +  (localPage.height-5)  + " moveto -90 rotate ("  +  (r)  +  ") show 90 rotate\n");
			}



		for(int i=0;i< this.genes.size();++i)
			{
			final Transcript g=this.genes.get(i);
			this.outw.println(  "gsave");
			this.outw.println( "0 " +  (localPage.height - margin.top-(fHeight*i)) + " translate");


			double x1=toPixel(g.getTxStart());
			double x2=toPixel(g.getTxEnd());
			this.outw.print(  "0 0 0 setrgbcolor\n");
			NEWPATH();
			MOVETO(x1,midy);
			LINETO(x2,midy);
			STROKE();
			//draw ticks

			this.outw.print(  "0.2 setlinewidth\n");
			this.outw.print(  "newpath\n");
			this.outw.print(  x1  + " "  +  midy  +  " moveto\n");
			this.outw.print(  x1  +  " "  +  x2  +  (g.isPositiveStrand()?" forticksF":" forticksR") +  "\n");
			this.outw.print(  "closepath stroke\n");


			this.outw.print(  "0.5 setlinewidth\n");
			if(g.hasCodonStartDefined() && g.hasCodonStopDefined()) {
			//draw txStart/txEnd
			this.outw.print(  "0.1 0.1 0.5 setrgbcolor\n"+
					"newpath\n"+
					 +  toPixel(g.getCodonStart().get().getStart())  +  " "+
					 +  (midy-cdsHeight/2.0)  + " "
					 +  (toPixel(g.getCodonStop().get().getEnd())-toPixel(g.getEnd()))  +  " "
					 +  cdsHeight  +  " box closepath fill\n"
					 )
					;
			}
			//draw each exon
			for(int j=0;j< g.getExonCount();++j)
				{
				this.outw.print(  toPixel(g.getExon(j).getStart())  +  " "
					  +  (midy-exonHeight/2.0) +  " "
					  +  (float)(toPixel(g.getExon(j).getEnd())-toPixel(g.getExon(j).getStart()))  +  " "
					  +  exonHeight  +  " gradient\n"
					  )
					 ;
				}
			//draw name
			this.outw.print(  "0 0 0 setrgbcolor\n");
			this.outw.print(  "10 "  +  midy  +  " moveto ("  +  g.getId()  +  ") show\n");
			this.outw.println(  "grestore");
			}
		
		//samples
			{
			double y= localPage.height-margin.top-(fHeight*(this.genes.size()+1));
			for(String sample:this.sample2positions.keySet())
				{
				
				this.outw.print(  "0.2 setlinewidth\n");
				this.outw.print(  "0 0 0 setrgbcolor\n");
				this.outw.print(  "10 "  +  (y-midy+5)  +  " moveto ("  +  sample  +  ") show\n");
				this.outw.print(  "newpath " +  margin.left  +  " " +  y  +  " moveto\n"
						 +  pageDef.width +  " "  +  0  +  " rlineto stroke\n"
						 );
					for(Integer r2 : this.sample2positions.get(sample))
					{
					
					this.outw.print(  "0.8 setlinewidth\n");
					this.outw.print(  "newpath " + toPixel(r2) + " " +  y  +  " circle closepath stroke\n");
					}
				y-=fHeight;
				}
			}
			
		this.outw.println("grestore");
		this.outw.print(  "showpage\n");

		}
		


	
	private void MOVETO(Object x,Object y)
		{
		this.outw.print(x);
		this.outw.print(' ');
		this.outw.print(y);
		this.outw.println(" moveto");
		}
	
	private void LINETO(Object x,Object y)
		{
		this.outw.print(x);
		this.outw.print(' ');
		this.outw.print(y);
		this.outw.println(" lineto");
		}
	private void STROKE()
		{
		this.outw.println("stroke");
		}
	
	private void NEWPATH()
		{
		this.outw.println("newpath");
		}

	


	private void run(final VCFIterator in)
		    {
		    for(;;)
			    {
		    	VariantContext ctx=null;
		    	if(in.hasNext())  ctx=in.next();
		    	
		    	if(ctx==null ||
		    			this.genes.isEmpty() ||
		    			(!this.genes.isEmpty() && !this.genes.get(0).getContig().equals(ctx.getContig())) || 
		    			(!this.genes.isEmpty() && this.chromEnd<=ctx.getStart())
		    			)
		    		{
		    		this.print();
		    		if(this.outw.checkError()) return;
		    		if(ctx==null) return;
		    		this.clear();

		    		
		    		for(Transcript g:chrom2transcript.getOrDefault(ctx.getContig(),Collections.emptyList()))
						{
						if(this.genes.isEmpty())
							{
							if(g.getTxEnd() <=ctx.getStart() || g.getTxStart()> ctx.getEnd() )
								{
								continue;
								}
							this.addGene(g);
							}
						else
							{
							if(!(g.getTxStart()>this.chromEnd || g.getTxEnd()<= this.chromStart))
								{
								this.addGene(g);
								}
							}
		
						}
		    			
		    		}
		    	
		    	if(!genes.isEmpty() &&
		    		ctx.getStart()-1>=this.chromStart &&
		    		ctx.getStart()<= this.chromEnd)
			    	{
			    	this.addVariant(ctx);
			    	}
			    }

		    
		    }

	private PrintStream outw = System.out;

	@Override
	public int doWork(List<String> args) {
		VCFIterator iter=null;
		BufferedReader r=null;
		try
			{
			iter = super.openVCFIterator( oneFileOrNull(args));
			this.outw = super.openPathOrStdoutAsPrintStream(this.outputFile);
			final SAMSequenceDictionary dict=iter.getHeader().getSequenceDictionary();

			if(this.gtfPath!=null) {
				try(GtfReader gtfReader = new GtfReader(this.gtfPath)) {
					if(dict!=null) gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					this.chrom2transcript.putAll(gtfReader.getAllGenes().
							stream().
							flatMap(G->G.getTranscripts().stream()).
							collect(Collectors.groupingBy(X->X.getContig())));
						
					}
				}
			

		    final double ticksH=(fHeight/2.0f)*0.6f;
			final double ticksx=20;

		    this.outw.print( 
		    		"%!PS\n"+
					"%%Creator: Pierre Lindenbaum PhD plindenbaum@yahoo.fr http://plindenbaum.blogspot.com\n"+
					"%%Title: "  +  getClass().getSimpleName()  +  "\n"+
					"%%CreationDate: " +  new Date()   +  "\n"+
					"%%EndComments\n"+
					"%%BoundingBox: 0 0 "  +(this.pageDef.width+margin.left+margin.right) +  " "  + (this.pageDef.height+margin.top+margin.bottom) + "\n"+
					"/Courier findfont 9 scalefont setfont\n"+
					"/circle { 10 0 360 arc} bind def\n"+
					"/ticksF {\n"
					 +  (-ticksH) +  " " +  (ticksH)  +  " rmoveto\n"
					 +  ticksH  +  " " +  (-ticksH)  +  " rlineto\n"
					 +  (-ticksH)  +  " " +  (-ticksH)  +  " rlineto\n"
					 +  ticksH  +  " " +  ticksH  +  " rmoveto\n"
					+ "} bind def\n"
					+ "/ticksR {\n"
					 +  (ticksH) +  " " +  (ticksH)  +  " rmoveto\n"
					 +  (-ticksH)  +  " " +  (-ticksH)  +  " rlineto\n"
					 +  (ticksH)  +  " " +  (-ticksH)  +  " rlineto\n"
					 +  (-ticksH)  +  " " +  (ticksH)  +  " rmoveto\n"
					 + "} bind def\n"+
					 "/forticksF {2 dict begin\n"
					+ "/x2 exch def\n"
					+ "/x1 exch def\n"
					+ "0 1 x2 x1 sub " +  ticksx  + " div {\n"
					+ "ticksF " +  ticksx  +  " 0 rmoveto\n"
					+ "}for\n"
					+ "} bind def\n"

					+ "/forticksR {2 dict begin\n"
					+ "/x2 exch def\n"
					+ "/x1 exch def\n"
					+ "0 1 x2 x1 sub " +  ticksx  + " div {\n"
					+ " ticksR  " +  ticksx  +  " 0 rmoveto\n"
					+ "}for\n"
					+ "} bind def\n"

					+ "/box\n"
					+ "{\n"
					+ "4 dict begin\n"
					+ "/height exch def\n"
					+ "/width exch def\n"
					+ "/y exch def\n"
					+ "/x exch def\n"
					+ "x y moveto\n"
					+ "width 0 rlineto\n"
					+ "0 height rlineto\n"
					+ "width -1 mul 0 rlineto\n"
					+ "0 height -1 mul rlineto\n"
					+ "end\n"
					+ "} bind def\n"

					+ "/gradient\n"
					+ "{\n"
					+ "4 dict begin\n"
					+ "/height exch def\n"
					+ "/width exch def\n"
					+ "/y exch def\n"
					+ "/x exch def\n"
					+ "/i 0 def\n"
					+ "height 2 div /i exch def\n"
					+ "\n"
					+ "0 1 height 2 div {\n"
					+ "	1 i height 2.0 div div sub setgray\n"
					+ "	newpath\n"
					+ "	x  \n"
					+ "	y height 2 div i sub  add\n"
					+ "	width\n"
					+ "	i 2 mul\n"
					+ "	box\n"
					+ "	closepath\n"
					+ 	"	fill\n"
					+ "	i 1 sub /i exch def\n"
					+ "	}for\n"
					+ "newpath\n"
					+ "0 setgray\n"
					+ "0.4 setlinewidth\n"
					+ "x y width height box\n"
					+ "closepath\n"
					+ "stroke\n"
					+ "end\n"
					+ "} bind def\n"
					);
		    
			run(iter);
			iter.close();
			this.outw.print("\n%%Trailer\n%%EOF\n"); 
			this.outw.flush();
			this.outw.close();
			this.outw=null;
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(iter);
			CloserUtil.close(this.outw);
			}
		
		}
	public static void main(final String[] args) {
		new VcfToPostscript().instanceMainWithExit(args);

	}

}
