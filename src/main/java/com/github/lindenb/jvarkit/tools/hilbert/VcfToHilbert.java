/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

import java.awt.geom.Point2D;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.w3c.dom.Element;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.hilbert.HilbertCurve;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
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
	private static final String SIGNATURE="HILBERT_";
	
    /** level of recursion */
    @Parameter(names={"-r","--recursion"},description="Hilbdert Curve level of recursion")
    private int recursionLevel=6;
    /** with/height of the final picture */
    @Parameter(names={"-w","--width"},description="Image width")
    private int imageWidth=1000;
    /** radius of a point */
    @Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
    private Path imgOut =null;
    @Parameter(names={"-B","--bed"},description="reduce input to this BED file. 4th column can be used to name the interval")
    private Path bedIn =null;
	@DynamicParameter(names={"-D"},description="other parameters.")
	private Map<String, String> dynaParams = new HashMap<>();

	private class Region extends SimpleInterval {
		long index=0L;
		String name="";
		String color = "black";
		Region(Locatable delegate) {
			super(delegate);
			}
		
		}
		
    @Override
    public int doWork(final List<String> args) {
		if(this.imageWidth<1)
			{
			LOG.error("Bad image size:" +this.imageWidth);
			return -1;
			}
		
		try(VCFIterator iter=this.openVCFIterator(oneFileOrNull(args))) {
			final VCFHeader header=iter.getHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final Hyperlink hyperlink = Hyperlink.compile(dict);
			final List<Region> intervals = new ArrayList<>();
			if(bedIn==null) {
				for(SAMSequenceRecord ssr:dict.getSequences()) {
					final Region r=new Region(ssr);
					r.name = ssr.getSequenceName();
					intervals.add(r);
					}
				}
			else
				{
				try(BedLineReader br =new BedLineReader(this.bedIn)) {
					while(br.hasNext()) {
						final BedLine bl=br.next();
						if(dict.getSequence(bl.getContig())==null) {
							LOG.warning("contig not in dictionary "+bl.getContig());
							continue;
							}
						final Region r=new Region(bl);
						r.name= bl.getOrDefault(3,bl.toNiceString());
						intervals.add(r);
						}
					}
				Collections.sort(intervals,new ContigDictComparator(dict).createLocatableComparator());
				}
			
			final long genomeLength=intervals.stream().mapToLong(it->it.getLengthOnReference()).sum();
			if(genomeLength<=0L) {
				LOG.error("no interval");
				return -1;
				}
			final IntervalTreeMap<Region> intervalTreeMap=new IntervalTreeMap<>();
			for(Region interval:intervals) {
				intervalTreeMap.put(interval.toInterval(), interval);
				}
			
			long pos=1;
			for(int i=0;i< intervals.size();i++) {
				final Region r=intervals.get(i);
				r.index=pos;
				pos+=r.getLengthOnReference();
				
				if(r.getContig().matches("(chr)?[X]")) {
					r.color="blue";
					}
				else if(r.getContig().matches("(chr)?[Y]")) {
					r.color="pink";
					}
				else if(i%2==0)
					{
					r.color="rgb(205,133,63)";
					}
				else
					{
					r.color="rgb(100,100,100)";
					}
				
				}
			
			final HilbertCurve hilbertCurve = new HilbertCurve(this.imageWidth,genomeLength, this.recursionLevel);
			final SVGDocument svgDoc = new SVGDocument();
			svgDoc.setHeight(this.imageWidth+1);
			svgDoc.setWidth(this.imageWidth+1);
			
			for(Region r:intervals) {
	        	final List<Point2D.Double> points = hilbertCurve.getPoints(r.index, r.index+r.getLengthOnReference());
	        	if(!points.isEmpty()) {
		        	final Element E = svgDoc.polyline(Maps.of(
		        			"stroke-width",dynaParams .getOrDefault("contig.stroke-width","3"),
		        			"stroke",r.color
		        			));
		        	svgDoc.rootElement.appendChild(E);
		        	E.setAttribute("points",svgDoc.toString(points));
		        	svgDoc.setTitle(E,r.name);
		        	}
				}
			
			
			while(iter.hasNext())
				{
				final VariantContext ctx= iter.next();
				if(!ctx.isVariant()) continue;
				for(Region region : intervalTreeMap.getOverlapping(ctx)) {
					final int start = Math.max(ctx.getStart(), region.getStart());
					final int end   = Math.min(ctx.getEnd(), region.getEnd());
					long n1 = region.index + start;
					long n2 = region.index + end;
							
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
					
					final Map<String,String> variantProps = new HashMap<>();
					ctx.getAttributes().entrySet().stream().
						filter(K->K.getKey().startsWith(SIGNATURE)).
						forEach(KV->{
							final String k = KV.getKey().substring(SIGNATURE.length()).toLowerCase().replace('_', '-');
							final String v = String.valueOf(KV.getValue());
							variantProps.put(k, v);
							});
					
					final String fill = variantProps.getOrDefault("fill",this.dynaParams.getOrDefault("fill", "yellow"));
					final String stroke = variantProps.getOrDefault("stroke",this.dynaParams.getOrDefault("stroke", "black"));
					final String opacity = variantProps.getOrDefault("opacity",this.dynaParams.getOrDefault("opacity", "0.6"));
					final String title = variantProps.getOrDefault("title",ctx.getContig()+":"+ctx.getStart());
					final String url = variantProps.getOrDefault("url",hyperlink.apply(ctx).orElse(""));

					Element E;
					if(length<=2) {
						final String radius = variantProps.getOrDefault("radius",this.dynaParams.getOrDefault("radius", "10"));
						E= svgDoc.circle(x,y,Double.parseDouble(radius));
						}
					else
						{
			        	E = svgDoc.polyline(Maps.of( "stroke-width",5));
			        	E.setAttribute("points",svgDoc.toString(points));
						}
					E.setAttribute("fill",fill);
					E.setAttribute("stroke",stroke);
					E.setAttribute("fill-opacity",opacity);

					svgDoc.setTitle(E,title);
					
					E = svgDoc.anchor(E,url);
					svgDoc.rootElement.appendChild(E);
					}
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
	