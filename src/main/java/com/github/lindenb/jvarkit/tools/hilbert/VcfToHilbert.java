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
package com.github.lindenb.jvarkit.tools.hilbert;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Point2D.Double;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

##Example

```bash
$  curl -s "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140123_NA12878_Illumina_Platinum/NA12878.wgs.illumina_platinum.20140404.snps_v2.vcf.gz" | gunzip -c |\
 java -jar dist/vcf2hilbert.jar  -r 1.1 -w 1000 -o hilbert.png
```

END_DOC

 */
@Program(name="vcf2hilbert",
	keywords={"vcf","image","visualization"},
	description="Plot a Hilbert Curve from a VCF file.")
public class VcfToHilbert extends Launcher
	{
	
	private static final Logger LOG=Logger.build(VcfToHilbert.class).make();

	/** graphics context */ 
	private Graphics2D g;
	/** dictionary */
	private SAMSequenceDictionary dict;
    /** level of recursion */
    private int recursionLevel=6;
    /** with/height of the final picture */
    @Parameter(names={"-w","--width"},description="Image width")
    private int imageWidth=1000;
    private double sampleWidth=0;
    private double genomicSizePerCurveUnit=0L;
    /** radius of a point */
    @Parameter(names={"-r","--radius"},description="Radius Size")
    private float radiusSize =3.0f;
    @Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
    private File imgOut =null;
    
    private abstract class HilbertSegmentHandler
    	{
        /** last index of the position in the genome */
    	double prev_base=0;
        /** size (in pb) of an edge */
        private final double d_base=VcfToHilbert.this.genomicSizePerCurveUnit;
    	/** size of an edge */
        private double dist=-1;
        /** last time we plot a segment, point at end */
        private Point2D.Double prevPoint;
        
        
        
        protected HilbertSegmentHandler()
        	{
        	this.dist = VcfToHilbert.this.sampleWidth;
        	for (int i= VcfToHilbert.this.recursionLevel; i>0; i--)
	            {
	            this.dist /= 2.0;
	            }
        	this.prevPoint=new Point2D.Double(
        			this.dist/2.0,
        			this.dist/2.0
        			);
            this.prev_base=0L;
        	}
        
        abstract void segment(
        		Point2D.Double beg,
        		Point2D.Double end,
        		long chromStart,
        		long chromEnd
        		);
        
        void run()
        	{
        	this.HilbertU(recursionLevel);
        	}	
        
        private void lineRel(double deltaX, double deltaY)
        	{
        	Point2D.Double point=new Point2D.Double(
        		this.prevPoint.x + deltaX,
        		this.prevPoint.y + deltaY
        		);
            long chromStart=(long)prev_base;
            long chromEnd=(long)(chromStart+ d_base);
            segment(prevPoint, point, chromStart, chromEnd);
            this.prevPoint= point;
            this.prev_base= chromEnd;;
        	}
 	//make U shaped curve       
        private void HilbertU(int level)
            {
            if (level <= 0) return;
            HilbertD(level-1);    this.lineRel(0, dist);
            HilbertU(level-1);    this.lineRel(dist, 0);
            HilbertU(level-1);    this.lineRel(0, -dist);
            HilbertC(level-1);
            }
     
	//make D shaped rule
        private void HilbertD(int level)
        	{
        	if (level <= 0) return;
            HilbertU(level-1);    this.lineRel(dist, 0);
            HilbertD(level-1);    this.lineRel(0, dist);
            HilbertD(level-1);    this.lineRel(-dist, 0);
            HilbertA(level-1);
    		}
     
	// make C shaped
        private void HilbertC(int level)
        	{
        	if (level <= 0) return;
            HilbertA(level-1);    this.lineRel(-dist, 0);
            HilbertC(level-1);    this.lineRel(0, -dist);
            HilbertC(level-1);    this.lineRel(dist, 0);
            HilbertU(level-1);
        	}
     	//make A shaped
        private void HilbertA(int level) {
        	if (level <= 0) return;
            HilbertC(level-1);    this.lineRel(0, -dist);
            HilbertA(level-1);    this.lineRel(-dist, 0);
            HilbertA(level-1);    this.lineRel(0, dist);
            HilbertD(level-1);
        	}
    	}	
    private class EvalCurve extends HilbertSegmentHandler
    	{
    
    	long count=0;
    	@Override
    	void segment(Double beg, Double end, long chromStart, long chromEnd) {
    		count++;
    		}
    	}
    
    private class PaintSegment extends HilbertSegmentHandler
    	{
    	long genomicStart;
    	long genomicEnd;
    	@Override
    	void segment(
         		final Point2D.Double beg,
         		final Point2D.Double end,
         		final long chromStart,
         		final long chromEnd
         		)
    		{
    		if(this.genomicStart > chromEnd) return;
    		if(this.genomicEnd < chromStart) return;
    		
    		long p0= Math.max(this.genomicStart,chromStart);
    		long p1= Math.min(this.genomicEnd,chromEnd);
    		double segsize=chromEnd-chromStart;
    		if(segsize==0) segsize=1.0;
    		
    		double r0= (p0-(double)chromStart)/segsize;
    		double r1= (p1-(double)chromStart)/segsize;
    		
    		double x0= beg.x + (end.x - beg.x)*r0;
    		double y0= beg.y + (end.y - beg.y)*r0;
    		double x1= beg.x + (end.x - beg.x)*r1;
    		double y1= beg.y + (end.y - beg.y)*r1;
    		    		
    		g.draw(new Line2D.Double(x0,y0,x1,y1));
    		}
    	}
    
    private class PaintVariant extends HilbertSegmentHandler
		{
		VariantContext ctx;
		private long genomic_index;
		private int sample_x;
		private int sample_y;
		PaintVariant(VariantContext ctx,int sample_x,int sample_y)
			{
			this.ctx=ctx;
			this.sample_x=sample_x;
			this.sample_y=sample_y;
			for(SAMSequenceRecord ssr:VcfToHilbert.this.dict.getSequences())
				{
				if(ssr.getSequenceName().equals(ctx.getContig()))
					{
					genomic_index+=ctx.getStart();
					break;
					}
				else
					{
					genomic_index+=ssr.getSequenceLength();
					}
				}
			}
		private void mul(double x,double y)
			{
			g.draw(new Line2D.Double(x-radiusSize, y-radiusSize, x+(radiusSize*2.0), y+(radiusSize*2.0)));
			g.draw(new Line2D.Double(x+radiusSize, y-radiusSize, x-(radiusSize*2.0), y+(radiusSize*2.0)));
			}
		private void cross(double x,double y)
			{
			g.draw(new Line2D.Double(x-radiusSize, y, x+radiusSize, y));
			g.draw(new Line2D.Double(x, y-radiusSize, x, y+radiusSize));
			}
		private void circle(double x,double y)
			{
			g.fill(new Ellipse2D.Double(x-radiusSize, y-radiusSize,(radiusSize*2.0), (radiusSize*2.0)));
			}
		
		@Override
		void segment(
	     		final Point2D.Double beg,
	     		final Point2D.Double end,
	     		final long chromStart,
	     		final long chromEnd
	     		)
			{
			if(genomic_index > chromEnd) return;
			if(genomic_index < chromStart) return;
			
			double segsize=chromEnd-chromStart;
			if(segsize==0) segsize=1.0;
			double r0= (genomic_index-(double)chromStart)/segsize;
			double x = beg.x + (end.x - beg.x)*r0;
			double y = beg.y + (end.y - beg.y)*r0;
			
			Genotype gx=this.ctx.getGenotype(sample_x);
			
	    	Composite composite=g.getComposite();
	    	g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));

			if(this.sample_x == this.sample_y 
				)
				{
				if( gx.isCalled() &&
					!gx.isHomRef())
					{
					boolean uniq_in_context=true;
					for(int i=0;i< ctx.getNSamples();++i)
						{
						if(i==this.sample_x) continue;
						Genotype other= ctx.getGenotype(i);
						if(other.sameGenotype(gx))
							{
							uniq_in_context=false;
							break;
							}
						}
					g.setColor(uniq_in_context?Color.RED:Color.BLACK);
					if(gx.isHet())
						{
						cross(x,y);
						}
					else
						{
						circle(x,y);
						}	
					}
				}
			else
				{
				Genotype gy = this.ctx.getGenotype(sample_y);
				if( gx.isCalled() &&
					gy.isCalled())
					{
					if(!(gx.sameGenotype(gy) && gx.isHomRef())) 
						{
						boolean x_uniq_in_context=true;
						boolean y_uniq_in_context=true;
						for(int i=0; i< ctx.getNSamples();++i)
							{
							if(i==this.sample_x || i==this.sample_y) continue;
							Genotype other= ctx.getGenotype(i);
							if(other.sameGenotype(gx))
								{
								x_uniq_in_context=false;
								}
							if(other.sameGenotype(gy))
								{
								y_uniq_in_context=false;
								}
							}
						g.setColor(Color.GRAY);
						if(y_uniq_in_context && x_uniq_in_context && gx.sameGenotype(gy))
							{
							g.setColor(Color.RED);
							circle(x,y);
							}
						else if(y_uniq_in_context && x_uniq_in_context)
							{
							g.setColor(Color.ORANGE);
							mul(x,y);
							}
						else if(gx.sameGenotype(gy))
							{ 
							g.setColor(Color.GRAY);
							cross(x,y);
							}
						else
							{
							g.setColor(Color.LIGHT_GRAY);
							cross(x,y);
							}
						}
					else
						{
						// both HOM-REF, ignore
						}
					
					}
				else if(gx.isCalled())
					{
					g.setColor(Color.LIGHT_GRAY);
					cross(x,y);
					}
				}
			//g.draw(new Line2D.Double(x0,y0,x1,y1));
	
			g.setComposite(composite);
			}
		}
    
	
	

    
    private void paintReference()
    	{
    	int n_colors=this.dict.getSequences().size();
    	ArrayList<Color> colors=new ArrayList<Color>(n_colors);
    	for(int i=0;i<n_colors;++i)
    		{
    		float gray=(i/(float)n_colors)*0.6f;
    		Color c=new Color( gray, gray, i%2==0?gray:1f-gray);
    		if(i%2==0)
    			{
    			colors.add(c);
    			}
    		else
    			{
    			colors.add(0,c);
    			}
    		}
    	Stroke oldStroke= g.getStroke();
    	Composite composite=g.getComposite();
    	g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.7f));
    	long pos=0L;
    	for(int tid=0;tid<this.dict.getSequences().size();++tid)
    		{
    		g.setStroke(new BasicStroke(0.5f));
    		SAMSequenceRecord ssr= this.dict.getSequence(tid);
    		g.setColor(colors.get(tid%colors.size()));
    		
    		PaintSegment paintSeg=new PaintSegment();
        	paintSeg.genomicStart = pos;
        	paintSeg.genomicEnd   = pos+ssr.getSequenceLength();
        	paintSeg.run();
    		
    		pos+=ssr.getSequenceLength();
    		}
    	g.setComposite(composite);
    	 g.setStroke(oldStroke);
    	}
    
    
    
    @Override
    public int doWork(final List<String> args) {
		if(this.imgOut==null)
			{
			LOG.error("output image file not defined");
			return -1;
			}
		
		if(this.imageWidth<1)
			{
			LOG.error("Bad image size:" +this.imageWidth);
			return -1;
			}
		VCFIterator iter=null;
		try
			{
			iter = this.openVCFIterator(oneFileOrNull(args));
			
			final VCFHeader header=iter.getHeader();
			this.dict = header.getSequenceDictionary();
			if(this.dict == null)
				{
				throw new JvarkitException.FastaDictionaryMissing("no dict in input");
				}
			final List<String> samples=header.getSampleNamesInOrder();
			if(samples.isEmpty())
				{
				throw new JvarkitException.SampleMissing("no.sample.in.vcf");
				}
			LOG.info("N-Samples:"+samples.size());
			double marginWidth = (this.imageWidth-2)*0.05;
			this.sampleWidth= ((this.imageWidth-2)-marginWidth)/samples.size();
			LOG.info("sample Width:"+sampleWidth);
			BufferedImage img = new BufferedImage(this.imageWidth, this.imageWidth, BufferedImage.TYPE_INT_RGB);
			this.g = (Graphics2D)img.getGraphics();
			this.g.setColor(Color.WHITE);
			this.g.fillRect(0, 0, imageWidth, imageWidth);
			
			g.setColor(Color.BLACK);
			final Hershey hershey =new Hershey();

			EvalCurve evalCurve=new EvalCurve();
			evalCurve.run();
			this.genomicSizePerCurveUnit = ((double)dict.getReferenceLength()/(double)(evalCurve.count));
			if(this.genomicSizePerCurveUnit<1) this.genomicSizePerCurveUnit=1;
			LOG.info("genomicSizePerCurveUnit:"+genomicSizePerCurveUnit);
			
			for(int x=0;x< samples.size();++x)
				{
				String samplex =  samples.get(x);
				
				double labelHeight=marginWidth;
				if(labelHeight>50) labelHeight=50;
				
				g.setColor(Color.BLACK);
				hershey.paint(g,
						samplex,
						marginWidth + x*sampleWidth,
						marginWidth - labelHeight,
						sampleWidth*0.9,
						labelHeight*0.9
						);
				
        		AffineTransform old=g.getTransform();
        		AffineTransform tr= AffineTransform.getTranslateInstance(
        				marginWidth ,
        				marginWidth + x*sampleWidth
        				);
        		tr.rotate(Math.PI/2);
        		g.setTransform(tr);
        		hershey.paint(g,
        				samplex,
						0.0,
						0.0,
						sampleWidth*0.9,
						labelHeight*0.9
						);        		//g.drawString(this.tabixFile.getFile().getName(),0,0);
        		g.setTransform(old);
				
				
				double tx = marginWidth+x*sampleWidth;
				for(int y=0;y< samples.size();++y)
					{
					double ty = marginWidth+y*sampleWidth;
					g.translate(tx, ty);
					g.setColor(Color.BLUE);
					g.draw(new Rectangle2D.Double(0, 0, sampleWidth, sampleWidth));
					//paint each chromosome
					paintReference();	
					g.translate(-tx, -ty);
					}
				}
			LOG.info("genomicSizePerCurveUnit:"+(long)genomicSizePerCurveUnit*evalCurve.count+" "+dict.getReferenceLength()+" count="+evalCurve.count);
			LOG.info("Scanning variants");
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.dict);
			while(iter.hasNext())
				{
				final VariantContext var=progress.watch(iter.next());
				for(int x=0;x< samples.size();++x)
					{
					final String samplex =  samples.get(x);
					final Genotype gx = var.getGenotype(samplex);
					if(!gx.isCalled()) continue;
					double tx = marginWidth+x*sampleWidth;

					for(int y=0;y< samples.size();++y)
						{
						final String sampley =  samples.get(y);
						final Genotype gy = var.getGenotype(sampley);
						if(!gy.isCalled()) continue;
						if(gx.isHomRef() && gy.isHomRef()) continue;
						double ty = marginWidth+y*sampleWidth;
						g.translate(tx, ty);
						final PaintVariant paint=new PaintVariant(var, x, y);
						paint.run();
						g.translate(-tx, -ty);
						}
					}
				}
			progress.finish();
			this.g.dispose();
			
			//save file
			LOG.info("saving "+imgOut);
			if(imgOut!=null)
				{
				ImageIO.write(img, imgOut.getName().toLowerCase().endsWith(".png")
						?"PNG":"JPG", imgOut);
				}
			else
				{
				ImageIO.write(img, "PNG", stdout());
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	public static void main(String[] args) {
		new VcfToHilbert().instanceMainWithExit(args);
	}
}
