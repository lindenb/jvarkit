package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.hershey.Hershey;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.RuntimeIOException;

public class BasesTrack extends Track {
	
	@Override
	public void paint(SVGContext ctx) {
		if(getReference()==null) return;
		
		double featureWidth=  ctx.image_width/(double)(ctx.loc.getLengthOnReference()); 
		double featureHeight= Math.min(Math.max(5.0,featureWidth),30); 			//base
		if(!isVisible() || featureWidth<5) return;
		final Hershey hershey = new Hershey();
		
		for(char base : new char[] {'A','C','G','T','N'}) {
			ctx.styleNode.appendChild(ctx.text(".b"+base+" {fill:none;stroke:"+AcidNucleics.cssColor(base)+";stroke-width:"+(featureWidth<5?0.5:1.0)+"px}\n"));
			}
		
		for(String base:new String[]{"a","A","t","T","g","G","c","C","n","N"})
			{
			final Element path  = ctx.element("path");
			ctx.defsNode.appendChild(path);
			
			path.setAttribute("id","b"+base);
			
			path.setAttribute("class","b"+base.toUpperCase());
			path.setAttribute("d",hershey.svgPath(
					base,
					0,
					0,
					featureWidth*0.95,
					featureHeight*0.95
					));
			path.appendChild(ctx.title(base));
			}
		final Path faidx = getReference().asPath();
		try(ReferenceSequenceFile refseqfile =  ReferenceSequenceFileFactory.getReferenceSequenceFile(faidx)) {
			final GenomicSequence genomicSequence = new GenomicSequence(refseqfile, ctx.loc.getContig());
			Element g = ctx.element("g");
			for(int pos= ctx.loc.getStart();
					featureWidth>5 && //ignore if too small
					pos<= ctx.loc.getEnd() && pos<=genomicSequence.length();++pos)
				{
				final char c= genomicSequence.charAt(pos-1);
				Element use = ctx.element("use");
				g.appendChild(use);
				use.setAttribute("x",format( ctx.pos2pixel(pos)));
				use.setAttribute("y",format(ctx.y));
				use.setAttribute("href", "#b"+c);
				use.appendChild(ctx.title(String.valueOf(c)+" "+StringUtils.niceInt(pos)));
				}
			ctx.tracksNode.appendChild(g);
			ctx.y+=featureHeight+1;
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}


}
