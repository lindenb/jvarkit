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


History:

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.TabixKnownGeneFileReader;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

## Example

```
# download and gzip refGene from UCSC
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz" |\
   gunzip -c | LC_ALL=C sort -t '	' -k3,3 -k5,5n | bgzip > refGene.txt.gz

# index the refGene file
$ tabix  -0 -s  3 -b 5 -e 6 -f refGene.txt.gz

# run vcf2svg 
$ java -jar dist/vcf2svg.jar \
   --stopAfterFirst \
   -k refGene.txt.gz input.vcf > out.svg
```

## Screenshot

https://twitter.com/yokofakun/status/851875435948462080

![screenshot](https://pbs.twimg.com/media/C9J4LeoXkAEqvIN.jpg)




END_DOC
 */
@Program(name="vcf2svg",
	description="write a vcf to svg , with gene context",
	keywords={"vcf","svg","xlm","visualization"}
		)
public class VcfToSvg extends Launcher {
private static final Logger LOG=Logger.build(VcfToSvg.class).make();
private static final String SEGMENT="__SEGMENT__";

@Parameter(names={"-o","--out"},description="Output SVG file. If defined, MUST Contains the word \'"+SEGMENT+"\'")
private File outputFile=null;

@Parameter(names={"-k","--knownGenes"},description=KnownGene.OPT_KNOWNGENE_DESC,required=true)
private String knownGeneUri=null;

@Parameter(names={"-trim2ctx","--trimToVariant"},
  description="Don't use gene interval for graphics but min/max of variants")
private boolean trimToVariants=false;

@Parameter(names={"-m","--manifest"},description="Manifest file containing the names of the files.")
private File manifestFile=null;
@Parameter(names={"--stopAfterFirst"},
description="when writing multiple SVG docs, stop after the first one. It avoids writing multiple concatenated SVG documents when writing to stdout")
private boolean stop_at_first=false;
@Parameter(names={"-gw","--genotypeWidth"},description="Genotype square width")
private int genotype_width=10;
@Parameter(names={"--nonCoding"},description="Ignore Non-coding genes")
private boolean removeNonCoding = false;
@Parameter(names={"--exon","--exons"},description="Only keep variants in exons")
private boolean variantsInExonOnly = false;
@Parameter(names={"--alphaFILTER"},description="Variant having FILTER!=PASS opacity (0== hide FILTERED variants)")
private double variantFILTEREDOpacity = 1.0;
@Parameter(names={"--alphaINDEL"},description="Variant INDEL opacity (0== hide INDEL variants) ")
private double variantIndelOpacity = 1.0;




private void title(final XMLStreamWriter w,final String title) throws XMLStreamException
	{
	w.writeStartElement("title");
	w.writeCharacters(title);
	w.writeEndElement();
	}

@Override
public int doWork(final List<String> args) {
	if(this.outputFile!=null && !outputFile.getName().contains(SEGMENT))
		{
		LOG.error("output file must contain the word "+SEGMENT+" :"+this.outputFile);
		return -1;
		}
	TabixKnownGeneFileReader tabix = null;
	VCFIterator r=null;
	OutputStream outputStream=null;
	XMLStreamWriter w=null;
	PrintWriter manifestW=null;
	try {
		
		LOG.info("opening knownGene ");
		tabix= new TabixKnownGeneFileReader(knownGeneUri);
		
		if(manifestFile!=null && this.outputFile!=null)
			{
			manifestW= new PrintWriter(manifestFile);
			}
		else
			{
			manifestW= new PrintWriter(new NullOuputStream());
			}

		
		final Set<String> chromosomes = tabix.getChromosomes();
		final XMLOutputFactory xof=XMLOutputFactory.newInstance();
		r= super.openVCFIterator(super.oneFileOrNull(args));
		final VCFHeader header=r.getHeader();
		while(r.hasNext())
			{
			final VariantContext ctx=r.next();
			String tabixContig=ctx.getContig();
			if(!chromosomes.contains(tabixContig))
				{
				if(tabixContig.startsWith("chr"))
					{
					tabixContig=tabixContig.substring(3);
					}
				else if(!tabixContig.startsWith("chr"))
					{
					tabixContig="chr"+tabixContig;
					}
				if(!chromosomes.contains(tabixContig))
					{
					while(r.hasNext())
						{
						final VariantContext ctx2 = r.peek();
						if(!ctx2.getContig().equals(ctx.getContig())) break;
						r.next();
						}
					LOG.error("No chromosome "+ctx.getContig()+" in "+knownGeneUri+". Check the chromosome nomenclature.");
					continue;
					}
				}
			
			final List<VariantContext> variants = new ArrayList<>();
			final List<KnownGene> genes = new ArrayList<>();
			variants.add(ctx);
		
			int chromStart=ctx.getStart()-1;
			int chromEnd = ctx.getEnd();
			/* walk over know gene, loop until there is no overapping transcript 
			 * over that region */
			for(;;) {
				genes.clear();
				/* the max chromEnd, let's see if we can get a bigger */
				int newStart=chromStart;
				int newEnd=chromEnd;
				final Iterator<KnownGene> kgr = tabix.iterator(
						tabixContig,
						chromStart,
						chromEnd
						);
				while(kgr.hasNext())
					{		
					final KnownGene g=kgr.next();
					if(this.removeNonCoding && g.isNonCoding()) continue;
					genes.add(g);
					newStart=Math.min(g.getTxStart(),newStart);
					newEnd=Math.max(g.getTxEnd(),newEnd);
					}
				if( newStart>=chromStart &&  newEnd <= chromEnd)
					{
					break;
					}
				chromStart=newStart;
				chromEnd=newEnd;
				}
			// intergenic, no gene over that variant
			if(genes.isEmpty()) continue;
			
			//fill the variant for that region
			while(r.hasNext())
				{
				final VariantContext ctx2 = r.peek();
				if(!ctx2.getContig().equals(ctx.getContig())) break;
				if(ctx2.getStart()> chromEnd) break;
				
				
				variants.add(r.next());
					
				}
			
			if(this.variantsInExonOnly)
				{
				variants.removeIf(V->{
					for(final KnownGene gene:genes)
						{
						for(final KnownGene.Exon exon: gene.getExons())
							{
							if(V.getEnd()< exon.getStart() || V.getStart()>=exon.getEnd())
								{
								//rien
								}
							else
								{
								
								return false;
								}
							}
						}
					return true;
					});
				}
			if(this.variantFILTEREDOpacity<=0)
				{
				variants.removeIf(V->V.isFiltered());
				}
			if(this.variantIndelOpacity<=0)
				{
				variants.removeIf(V->V.isIndel());
				}

			
			
			if(variants.isEmpty()) continue;
			
			LOG.info("Variants ("+variants.size()+") Transcripts ("+genes.size()+") "+tabixContig+":"+chromStart+"-"+chromEnd);

			if(outputFile!=null)
				{
				File fname=new File(outputFile.getParentFile(),
						outputFile.getName().replaceAll("__SEGMENT__",
						ctx.getContig()+"_"+chromStart+"_"+chromEnd));
				LOG.info("saving as "+fname);
				outputStream=IOUtils.openFileForWriting(fname);
				w=xof.createXMLStreamWriter(outputStream);
				
				manifestW.println(ctx.getContig()+"\t"+chromStart+"\t"+chromEnd+"\t"+
						genes.stream().map(G->G.getName()).collect(Collectors.joining(","))+"\t"+
						genes.size()+"\t"+
						variants.size()+"\t"+
						fname
						);

				
				}
			else
				{
				w=xof.createXMLStreamWriter(stdout());
				}
			
			
			double featureHeight=10;
			double TRANSCRIPT_HEIGHT=featureHeight; 

			final int all_genotypes_width = variants.size()*this.genotype_width;
            if(trimToVariants)
	        	{
	        	chromStart = variants.stream().map(V->V.getStart()-1).
	        			min((A,B)->A.compareTo(B)).get();
	        	chromEnd = variants.stream().map(V->V.getEnd()+1).
	        			max((A,B)->A.compareTo(B)).get();
	        	}	
	        final int drawinAreaWidth=Math.max(all_genotypes_width,1000);
	
	        final Interval interval=new Interval(ctx.getContig(),
	        		chromStart, 
	        		chromEnd
	        		);
			final int interline_weight=6;
			final int margin_top=10;
			final int margin_bottom=10;
			final int margin_right=100;
			final int margin_left=100;
			
			w.writeStartDocument("UTF-8", "1.0");
			
            w.writeStartElement("svg");
            w.writeDefaultNamespace(SVG.NS);
            w.writeNamespace("xlink", XLINK.NS);
            w.writeAttribute("version", "1.1");
            w.writeAttribute("width",String.valueOf(margin_right+margin_right+drawinAreaWidth));
            w.writeAttribute("height",String.valueOf(
            		margin_top+margin_bottom+
            		genes.size()*TRANSCRIPT_HEIGHT+
            		interline_weight*featureHeight+
            		header.getSampleNamesInOrder().size()*this.genotype_width
            		));
            title(w,ctx.getContig()+":"+chromStart+"-"+chromEnd);
            
            w.writeStartElement("desc");
            w.writeCharacters("generated with "+getProgramName()+"\n"+
            		"Author: Pierre Lindenbaum PhD. @yokofakun .");
            w.writeEndElement();
            
            //defs
            w.writeStartElement("defs");
			
    		//genotypes
    		
    		w.writeStartElement("g");
				w.writeAttribute("id","g_"+GenotypeType.HOM_REF); //
					w.writeEmptyElement("rect");
						w.writeAttribute("style","fill:lime;stroke;none;");
						w.writeAttribute("x","0");
						w.writeAttribute("y","0" );
						w.writeAttribute("width",String.valueOf(this.genotype_width));
						w.writeAttribute("height",String.valueOf(this.genotype_width));
    		w.writeEndElement();
    		
    		w.writeStartElement("g");
				w.writeAttribute("id","g_"+GenotypeType.NO_CALL); //
					w.writeEmptyElement("rect");
						w.writeAttribute("style","fill:silver;stroke;gray;");
						w.writeAttribute("x","0");
						w.writeAttribute("y","0" );
						w.writeAttribute("width",String.valueOf(this.genotype_width));
						w.writeAttribute("height",String.valueOf(this.genotype_width));
    		w.writeEndElement();
    		
			w.writeStartElement("g");
				w.writeAttribute("id","g_"+GenotypeType.HOM_VAR); //
					w.writeEmptyElement("rect");
						w.writeAttribute("style","fill:crimson;stroke;none;");
						w.writeAttribute("x","0");
						w.writeAttribute("y","0" );
						w.writeAttribute("width",String.valueOf(this.genotype_width));
						w.writeAttribute("height",String.valueOf(this.genotype_width));
    		w.writeEndElement();	
    		
    		w.writeStartElement("g");
			w.writeAttribute("id","g_"+GenotypeType.MIXED); //
				w.writeEmptyElement("rect");
					w.writeAttribute("style","fill:pink;stroke;none;");
					w.writeAttribute("x","0");
					w.writeAttribute("y","0" );
					w.writeAttribute("width",String.valueOf(this.genotype_width));
					w.writeAttribute("height",String.valueOf(this.genotype_width));
		    w.writeEndElement();	
    		
    		w.writeStartElement("g");
			w.writeAttribute("id","g_"+GenotypeType.UNAVAILABLE); //
				w.writeEmptyElement("rect");
					w.writeAttribute("style","fill:gray;stroke;none;");
					w.writeAttribute("x","0");
					w.writeAttribute("y","0" );
					w.writeAttribute("width",String.valueOf(this.genotype_width));
					w.writeAttribute("height",String.valueOf(this.genotype_width));
		    w.writeEndElement();	

		    
    		w.writeStartElement("g");
    		w.writeAttribute("id","g_"+GenotypeType.HET); //
    			w.writeEmptyElement("rect");
    			w.writeAttribute("style","fill:lime;stroke;black;");
    			w.writeAttribute("x", "0" );
    			w.writeAttribute("y", "0" );
        		w.writeAttribute("width",String.valueOf(genotype_width));
        		w.writeAttribute("height",String.valueOf(genotype_width));
    			w.writeEmptyElement("polygon");
    			w.writeAttribute("style","fill:crimson;stroke;black;");
    			w.writeAttribute("points","0,0 "+genotype_width+",0 0,"+genotype_width+" 0,0");
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
    		
    		w.writeEndElement();//defs

    		w.writeStartElement("style");
    		w.writeCharacters(
    				"svg {fill:none; stroke:black;}\n"+
    				"text {fill:black;stroke:none;font-size:"+ (featureHeight/1.5) +"px;}\n"+
    				".ruler-label { stroke:red;}\n"+
    				".frame { stroke:black;fill:none;}\n"+
    				".kgexon {fill:url(#grad01);stroke:black;}\n"+
    				".gcpercent {fill:url(#grad02);stroke:black;}"+
    				".coverage {fill:url(#grad03);stroke:black;}"+
    				".kgcds {fill:yellow;stroke:black;opacity:0.7;}\n"+
    				".variant{stroke:none;fill:red;opacity:0.2;}\n"+
    				".xaxis{stroke:gray;fill:none;opacity:0.2;}\n"+
    				".postick{font-size:9px;stroke:black;stroke-width:1;}"
    				);
    		w.writeEndElement();//style

            
            	
			final Function<Integer,Integer> trim= new Function<Integer,Integer>() {
				@Override
				public Integer apply(final Integer t) {
					return Math.max(interval.getStart(), Math.min(interval.getEnd(), t));
					}
			};
			final Function<Integer,Double> baseToPixel= new Function<Integer,Double>() {
				@Override
				public Double apply(final Integer t) {
					return  margin_left + drawinAreaWidth*(t-(double)interval.getStart())/((double)interval.length());
					}
			};
			final Function<Integer,Double> variantIndexToPixel= new Function<Integer,Double>() {
				@Override
				public Double apply(final Integer idx) {
					final double variant_width= drawinAreaWidth/(double)variants.size();
					final double midx=variant_width*idx+variant_width/2.0;
					return margin_left+ midx-genotype_width/2.0;
					}
			};
			
			final Function<VariantContext,String> variantTitle= V-> (V.getContig().startsWith("chr")?V.getContig().substring(3):V.getContig())+":"+V.getStart()+" "+V.getReference().getDisplayString();
			
			/** title */
			double y=0;
			w.writeStartElement("text");
			w.writeAttribute("x","0");
			w.writeAttribute("y",String.valueOf(featureHeight));
			w.writeCharacters(interval.toString());
			w.writeEndElement();
			y+= featureHeight;
			
			for(final KnownGene g:genes)
				{
				int cdsHeigh= 5;
				double exonHeight=TRANSCRIPT_HEIGHT-5;
				double midY=TRANSCRIPT_HEIGHT/2;
		
				w.writeStartElement("g");
				
				
				
				w.writeAttribute("transform", "translate(0,"+y+")");
				
				title(w, g.getName());
				
				w.writeStartElement("text");
				w.writeAttribute("x",String.valueOf(margin_left-10));
				w.writeAttribute("y",String.valueOf(featureHeight));
				w.writeAttribute("style","text-anchor:end;");
				w.writeCharacters(g.getName());
				w.writeEndElement();
				
				/* transcript line */
				w.writeEmptyElement("line");
				w.writeAttribute("class","kgtr");
				w.writeAttribute("x1",String.valueOf(baseToPixel.apply(trim.apply(g.getTxStart()))));
				w.writeAttribute("y1",String.valueOf(midY));
				w.writeAttribute("x2",String.valueOf(baseToPixel.apply(trim.apply(g.getTxEnd()))));
				w.writeAttribute("y2",String.valueOf(midY));
				
				
				
				
				/* strand symbols */
				for(double pixX=0;
					pixX< drawinAreaWidth;
					pixX+=30)
					{
					double pos0= interval.getStart()+(pixX/(double)drawinAreaWidth)*interval.length();
						if(pos0+1< g.getTxStart()) continue;
						if(pos0> g.getTxEnd()) break;
						w.writeEmptyElement("use");
						w.writeAttribute("class","kgstrand");
						w.writeAttribute("xlink", XLINK.NS, "href", "#strand"+(g.isPositiveStrand()?"F":"R"));
						w.writeAttribute("x",String.valueOf(margin_left + pixX));
						w.writeAttribute("y",String.valueOf(midY));
					}
			
				/* exons */
				for(KnownGene.Exon exon:g.getExons())
						{
						if(exon.getStart()+1>=  interval.getEnd()) continue;
						if(exon.getEnd()<= interval.getStart()) continue;
						w.writeStartElement("rect");
						w.writeAttribute("class","kgexon");
						
						w.writeAttribute("x",String.valueOf(baseToPixel.apply(trim.apply(exon.getStart()))));
						w.writeAttribute("y",String.valueOf(midY-exonHeight/2));
						w.writeAttribute("width",String.valueOf(baseToPixel.apply(trim.apply(exon.getEnd()))-baseToPixel.apply((trim.apply(exon.getStart())))));
						w.writeAttribute("height",String.valueOf(exonHeight));
						title(w,exon.getName());
						w.writeEndElement();
						}
			
				/* coding line */
				if(!g.isNonCoding())
					{
					w.writeEmptyElement("rect");
					w.writeAttribute("class","kgcds");
					w.writeAttribute("x",String.valueOf(baseToPixel.apply(trim.apply(g.getCdsStart()))));
					w.writeAttribute("y",String.valueOf(midY-cdsHeigh/4.0));
					w.writeAttribute("width",String.valueOf(baseToPixel.apply(trim.apply(g.getCdsEnd()))-baseToPixel.apply((trim.apply((g.getCdsStart()))))));
					w.writeAttribute("height",String.valueOf(cdsHeigh/2.0));
					}
				
				//String label=String.format("%15s", g.getName());
				//w.writeEmptyElement("path");
				//double fontHeight=Math.min(10,0.8*TRANSCRIPT_HEIGHT);
				//w.writeAttribute("d",this.hershey.svgPath(label,-insets.left,midY-fontHeight/2,insets.left*0.9,fontHeight));

				
				w.writeEndElement();
				w.writeCharacters("\n");
				y+=featureHeight;
				}
			
			/* draw lines to variants */
			for(int vidx=0;vidx < variants.size();++vidx)
				{
				final VariantContext vc=variants.get(vidx);
				double x1 =  baseToPixel.apply(vc.getStart());
				double x2 =  baseToPixel.apply(vc.getEnd());
				final double y2= y+featureHeight*interline_weight;
				w.writeStartElement("polygon");
				w.writeAttribute("style","fill:"+(vidx%2==0?"ghostwhite":"lavender")+";stroke:black;opacity:0.6;stroke-width:0.5;");
				w.writeAttribute("points",
						 ""+x1+","+(y-featureHeight/2.0)+
						" "+x2+","+(y-featureHeight/2.0)+
						" "+variantIndexToPixel.apply(vidx)+","+y2+
						" "+(variantIndexToPixel.apply(vidx)+this.genotype_width)+","+y2
						);
				title(w,variantTitle.apply(vc));
				w.writeEndElement();
				}
			for(int vidx=0;vidx < variants.size();++vidx)
				{
				final VariantContext vc=variants.get(vidx);
				final double y2= y+featureHeight*interline_weight;
				w.writeStartElement("text");
				w.writeAttribute("transform",
						"translate("+(String.valueOf(variantIndexToPixel.apply(vidx)+genotype_width/2.0))+","+String.valueOf(y2-5)+") "+
						"rotate(-45)");
				w.writeAttribute("x", "0");
				w.writeAttribute("y", "0");
				w.writeAttribute("class", "postick");
				w.writeCharacters(variantTitle.apply(vc));
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			y+=featureHeight*interline_weight;
			
			w.writeStartElement("g");
			for(final String sample: header.getSampleNamesInOrder())
				{
				for(int vidx=0;vidx < variants.size();++vidx)
					{
					final VariantContext vc=variants.get(vidx);
					final Genotype g=vc.getGenotype(sample);
					double opacity= 1.0;
					
					if(vc.isIndel())  opacity*=this.variantIndelOpacity;
					if(vc.isFiltered())  opacity*=this.variantFILTEREDOpacity;
					if(opacity>1) opacity=1;
					if(opacity<=0) continue;
					if(opacity<1)
						{
						w.writeStartElement("g");
						w.writeAttribute("style","opacity:"+opacity+";");
						}
					
					w.writeEmptyElement("use");
					w.writeAttribute("x",""+variantIndexToPixel.apply(vidx));
					w.writeAttribute("y",String.valueOf(y));
					w.writeAttribute("xlink", XLINK.NS, "href", "#g_"+g.getType());
					
					if(opacity<1)
						{
						w.writeEndElement();
						}
					
					}
				w.writeCharacters("\n");
				w.writeStartElement("text");
				w.writeAttribute("x",String.valueOf(margin_left-10));
				w.writeAttribute("y",String.valueOf(y + this.genotype_width/2.0));
				w.writeAttribute("style","text-anchor:end;");
				w.writeCharacters(sample);
				w.writeEndElement();
				
				y+=this.genotype_width;
				}
			w.writeCharacters("\n");
			
			
			
			w.writeEndDocument();
			w.writeCharacters("\n");
			w.flush();
			w.close();
			if(outputFile!=null)
				{
				outputStream.flush();
				outputStream.close();
				outputStream=null;
				}
			if(stop_at_first) 
				{
				LOG.info("Stop after first SVG document");
				break;
				}
			}
		r.close();
		manifestW.flush();
		manifestW.close();
		manifestW=null;
		return 0;
		}
	catch(Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(r);
		CloserUtil.close(tabix);
		CloserUtil.close(outputStream);
		CloserUtil.close(manifestW);
		}
	}
	
	
public static void main(String[] args) {
	
	new VcfToSvg().instanceMainWithExit(args);
	}
}
