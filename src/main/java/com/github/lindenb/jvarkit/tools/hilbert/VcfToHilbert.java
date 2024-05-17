/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
import java.awt.image.BufferedImage;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.w3c.dom.Element;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.hilbert.HilbertCurve;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.hershey.Hershey;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

##Example

```bash
$  curl -s "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140123_NA12878_Illumina_Platinum/NA12878.wgs.illumina_platinum.20140404.snps_v2.vcf.gz" | gunzip -c |\
 java -jar dist/jvarkit.jar vcf2hilbert > hilbert.svg
```

END_DOC

 */
@Program(name="vcf2hilbert",
	keywords={"vcf","image","visualization","svg"},
	description="Plot a Hilbert Curve from a VCF file as SVG",
	creationDate="20171201",
	modificationDate="20240517",
	jvarkit_amalgamion = true
	)
public class VcfToHilbert extends Launcher
	{
	
	private static final Logger LOG=Logger.build(VcfToHilbert.class).make();

    /** level of recursion */
    @Parameter(names={"-r","--recursion"},description="Hilbdert Curve level of recursion")
    private int recursionLevel=6;
    /** with/height of the final picture */
    @Parameter(names={"-w","--width"},description="Image width")
    private int imageWidth=1000;
    /** radius of a point */
    @Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
    private Path imgOut =null;
        
    /*
    private class PaintVariant implements HilbertCurve.Handler
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
			for(SAMSequenceRecord ssr:dict.getSequences())
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
		public void segment(
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
   
	
	

    
    private void paintReference(final HilbertCurve hilbertCurve,final List<Interval> intervals,final SVGDocument svgDoc)
    	{
    	final int n_colors= intervals.size();
    	final ArrayList<Color> colors=new ArrayList<Color>(n_colors);
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
    		
    		
    		
    		pos+=ssr.getSequenceLength();
    		}
    	g.setComposite(composite);
    	 g.setStroke(oldStroke);
    	}
     */
    
    
    @Override
    public int doWork(final List<String> args) {
		if(this.imageWidth<1)
			{
			LOG.error("Bad image size:" +this.imageWidth);
			return -1;
			}
		
		
		
		try(VCFIterator iter=this.openVCFIterator(oneFileOrNull(args))) {
			final VCFHeader header=iter.getHeader();
			SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final List<Interval> intervals = new ArrayList<>();
			for(SAMSequenceRecord ssr:dict.getSequences()) {
				intervals.add(new Interval(ssr.getContig(),1,ssr.getSequenceLength(),false,ssr.getSequenceName()));
				}
			
			final long genomeLength=intervals.stream().mapToLong(it->it.getLengthOnReference()).sum();
			LOG.info("genome length "+genomeLength);
			
			final HilbertCurve hilbertCurve = new HilbertCurve(this.imageWidth,genomeLength, this.recursionLevel);
			final SVGDocument svgDoc = new SVGDocument();
			svgDoc.setHeight(this.imageWidth+1);
			svgDoc.setWidth(this.imageWidth+1);
			long pos=1;
			for(int i=0;i< intervals.size();i++) {
				final Interval ssr = intervals.get(i);
				final int color[] = new int[] {0,0,0};
				if(ssr.getContig().matches("(chr)?[X]")) {
					color[2]=255;
					}
				else if(ssr.getContig().matches("(chr)?[Y]")) {
					color[0]=255;
					color[1]=192;
					color[2]=203;
					}
				else if(i%2==0)
					{
					color[0]=205;
					color[1]=133;
					color[2]=63;
					}
				else
					{
					color[0]=100;
					color[1]=100;
					color[2]=100;
					}
				
	        	final List<Point2D.Double> points = hilbertCurve.getPoints(pos, pos+ssr.getLengthOnReference());
	        	if(!points.isEmpty()) {
		        	Element E = svgDoc.polyline(Maps.of("stroke-width",3,"stroke","rgb("+color[0]+","+color[1]+","+color[2]+")"));
		        	svgDoc.rootElement.appendChild(E);
		        	E.setAttribute("points",svgDoc.toString(points));
		        	svgDoc.setTitle(E,ssr.getName());
		        	}
	        	pos+=ssr.getLengthOnReference();
				}
			
			while(iter.hasNext())
				{
				final VariantContext ctx= iter.next();
				if(!ctx.isVariant()) continue;
				long n1=0;
				long n2=0;
				for(SAMSequenceRecord ssr:dict.getSequences()) {
					if(ssr.contigsMatch(ctx)) {
						n2=n1+ctx.getEnd();
						n1+=ctx.getStart();
						break;
						}
					else
						{
						n1+=ssr.getLengthOnReference();
						}
					}
				final List<Point2D.Double> points = hilbertCurve.getPoints(n1,n2);
				if(points.isEmpty()) continue;
				
				double length = 0;
				double x = points.get(0).getX();
				double y = points.get(0).getY();
				for(int i=0;i+1< points.size();i++) {
					length+= points.get(i).distance(points.get(i+1));
					x+= points.get(i+1).getX();
					y+= points.get(i+1).getY();
					}
				x/=points.size();
				y/=points.size();
				Element E;
				if(length<=2) {
					E= svgDoc.circle(x,y,10);
					E.setAttribute("fill","yellow");
					E.setAttribute("stroke","black");
					}
				else
					{
		        	E = svgDoc.polyline(Maps.of("stroke-width",5,"stroke","yellow"));
		        	E.setAttribute("points",svgDoc.toString(points));
					}
				svgDoc.rootElement.appendChild(E);
				svgDoc.setTitle(E,ctx.getContig()+":"+ctx.getStart());
				}
			
			svgDoc.saveToFileOrStdout(this.imgOut);
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		
		}
    
	public static void main(final String[] args) {
		new VcfToHilbert().instanceMainWithExit(args);
		}
	}
	