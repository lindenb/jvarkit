package com.github.lindenb.jvarkit.tools.misc;

import java.awt.Dimension;
import java.awt.Insets;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfToPostscript extends AbstractCommandLineProgram
	{
	private List<KnownGene> genes=new ArrayList<KnownGene>();
	private Set<Integer> positions=new HashSet<Integer>();
	private Map<String,Set<Integer>> sample2positions=new TreeMap<String,Set<Integer>>();
	private int chromStart=Integer.MAX_VALUE;
	private int chromEnd=Integer.MIN_VALUE;

	private Map<String,List<KnownGene>> chrom2knownGenes=new HashMap<String,List<KnownGene>>();

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
	
	private void addGene(KnownGene g)
		{
		chromStart=Math.min(chromStart, g.getTxStart());
		chromEnd=Math.max(chromEnd, g.getTxEnd());
		genes.add(g);
		}
	
	private double toPixel(int pos)
		{
		return (double)margin.left+((pos-(double)this.chromStart)/((double)this.chromEnd-(double)this.chromStart))*(double)pageDef.width;
		}

			
	private void addVariant(final VariantContext ctx)
		{
		if(!ctx.getContig().equals(genes.get(0).getChromosome())) return;
		if(ctx.getStart()>=chromEnd) return;
		if(ctx.getStart()<chromStart) return;
		positions.add(ctx.getStart());
		for(String sample: ctx.getSampleNames())
			{
			Genotype g=ctx.getGenotype(sample);
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
		
		final KnownGene first=genes.get(0);
		final double fHeight=20;
		
		
		final Dimension localPage=new Dimension(
				(int)(this.pageDef.width),
				(int)(margin.top+margin.bottom + (this.genes.size()+1)*fHeight + this.sample2positions.size()*fHeight)
				);
		
		System.out.println(  "\n%%Page: "  +  count_pages_printed  +  " " + count_pages_printed);
		System.out.println("gsave");
		System.out.println(
				1.0 +" "+ 
				(localPage.height <= this.pageDef.height ? 1.0 : 1.0/(localPage.getHeight()/(float)this.pageDef.getHeight()))
				+ " scale");

		//System.out.println( "%%BoundingBox: 0 0 "  + page.width +  " "  + page.height  );
		

		float midy=(float)(fHeight/2.0f);

		double cdsHeight=fHeight*0.4;
		double exonHeight=fHeight*0.9;

		System.out.println(
				"2 "  +  (localPage.height-10 ) +
				" moveto ("  + first.getChromosome() + ":"+this.chromStart+ "-" +this.chromEnd+") show"
				 )
				;


		System.out.println("1 0 0 setrgbcolor");
		System.out.println("0.3 setlinewidth");
		for(Integer r: this.positions )
			{
			System.out.print(  "newpath " +  (float)toPixel(r)  +  " 0 moveto 0 " +  localPage.height  +  " rlineto stroke\n");
			System.out.print(  (float)toPixel(r)  +  " "  +  (localPage.height-5)  + " moveto -90 rotate ("  +  (r)  +  ") show 90 rotate\n");
			}



		for(int i=0;i< this.genes.size();++i)
			{
			KnownGene g=this.genes.get(i);
			System.out.println(  "gsave");
			System.out.println( "0 " +  (localPage.height - margin.top-(fHeight*i)) + " translate");


			double x1=toPixel(g.getTxStart());
			double x2=toPixel(g.getTxEnd());
			System.out.print(  "0 0 0 setrgbcolor\n");
			NEWPATH();
			MOVETO(x1,midy);
			LINETO(x2,midy);
			STROKE();
			//draw ticks

			System.out.print(  "0.2 setlinewidth\n");
			System.out.print(  "newpath\n");
			System.out.print(  x1  + " "  +  midy  +  " moveto\n");
			System.out.print(  x1  +  " "  +  x2  +  (g.isPositiveStrand()?" forticksF":" forticksR") +  "\n");
			System.out.print(  "closepath stroke\n");


			System.out.print(  "0.5 setlinewidth\n");
			//draw txStart/txEnd
			System.out.print(  "0.1 0.1 0.5 setrgbcolor\n"+
					"newpath\n"+
					 +  toPixel(g.getCdsStart())  +  " "+
					 +  (midy-cdsHeight/2.0)  + " "
					 +  (toPixel(g.getCdsEnd())-toPixel(g.getCdsStart()))  +  " "
					 +  cdsHeight  +  " box closepath fill\n"
					 )
					;
			//draw each exon
			for(int j=0;j< g.getExonCount();++j)
				{
				System.out.print(  toPixel(g.getExon(j).getStart())  +  " "
					  +  (midy-exonHeight/2.0) +  " "
					  +  (float)(toPixel(g.getExon(j).getEnd())-toPixel(g.getExon(j).getStart()))  +  " "
					  +  exonHeight  +  " gradient\n"
					  )
					 ;
				}
			//draw name
			System.out.print(  "0 0 0 setrgbcolor\n");
			System.out.print(  "10 "  +  midy  +  " moveto ("  +  g.getName()  +  ") show\n");
			System.out.println(  "grestore");
			}
		
		//samples
			{
			double y= localPage.height-margin.top-(fHeight*(this.genes.size()+1));
			for(String sample:this.sample2positions.keySet())
				{
				
				System.out.print(  "0.2 setlinewidth\n");
				System.out.print(  "0 0 0 setrgbcolor\n");
				System.out.print(  "10 "  +  (y-midy+5)  +  " moveto ("  +  sample  +  ") show\n");
				System.out.print(  "newpath " +  margin.left  +  " " +  y  +  " moveto\n"
						 +  pageDef.width +  " "  +  0  +  " rlineto stroke\n"
						 );
					for(Integer r2 : this.sample2positions.get(sample))
					{
					
					System.out.print(  "0.8 setlinewidth\n");
					System.out.print(  "newpath " + toPixel(r2) + " " +  y  +  " circle closepath stroke\n");
					}
				y-=fHeight;
				}
			}
			
		System.out.println("grestore");
		System.out.print(  "showpage\n");

		}
		


	/*
	private float _pageWidth()
		{
		return (float)pageDef.width;
		}

	private float _pageHeight()
		{
		return (float)pageDef.height;
		}
		

	private float width()
		{
		return pageDef.width-(margin.left+margin.right);
		}*/


	
	private void MOVETO(Object x,Object y)
		{
		System.out.print(x);
		System.out.print(' ');
		System.out.print(y);
		System.out.println(" moveto");
		}
	
	private void LINETO(Object x,Object y)
		{
		System.out.print(x);
		System.out.print(' ');
		System.out.print(y);
		System.out.println(" lineto");
		}
	private void STROKE()
		{
		System.out.println("stroke");
		}
	
	private void NEWPATH()
		{
		System.out.println("newpath");
		}

	


	private void run(VcfIterator in)
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
		    		if(System.out.checkError()) return;
		    		if(ctx==null) return;
		    		this.clear();

		    		if(chrom2knownGenes.containsKey(ctx.getContig()))
		    			{
		    			for(KnownGene g:chrom2knownGenes.get(ctx.getContig()))
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
		    			if(genes.isEmpty())
		    				{
		    				debug("no gene for "+ctx.getContig()+":"+ctx.getStart());
		    				}
		    			}
		    		else
		    			{
		    			debug("not any gene for "+ctx.getContig());
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


	@Override
	public String getProgramDescription() {
		return "Print VCF context as Postscript";
		}
	
	private String ucscKnownGene="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -k (uri/file) knownGene source. Required. Default: "+ucscKnownGene);
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"k:"))!=-1)
			{
			switch(c)
				{	
				case 'k':this.ucscKnownGene=opt.getOptArg();break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		VcfIterator iter=null;
		BufferedReader r=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				iter=new VcfIterator(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String uri=args[opt.getOptInd()];
				info("Reading from "+uri);
				iter=new VcfIterator(IOUtils.openURIForReading(uri));
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			info("Reading "+this.ucscKnownGene);
			SAMSequenceDictionary dict=iter.getHeader().getSequenceDictionary();
			
			int nKG=0;
			r=IOUtils.openURIForBufferedReading(ucscKnownGene);
			Pattern tab=Pattern.compile("[\t]");
			String line;
			while((line=r.readLine())!=null)
				{
				String tokens[]=tab.split(line);
				KnownGene g=new KnownGene(tokens);
				if(dict!=null && dict.getSequence(g.getContig())==null) continue;
				
				List<KnownGene> kg=this.chrom2knownGenes.get(g.getContig());
				if(kg==null)
					{
					kg=new ArrayList<KnownGene>();
					this.chrom2knownGenes.put(g.getContig(),kg);
					}
				kg.add(g);
				++nKG;
				}
			r.close();
			info("Done Reading knownGenes. N="+nKG);

		    final double ticksH=(fHeight/2.0f)*0.6f;
			final double ticksx=20;

		    System.out.print( 
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
			System.out.print("\n%%Trailer\n%%EOF\n"); 
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(iter);
			}
		
		}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new VcfToPostscript().instanceMainWithExit(args);

	}

}
