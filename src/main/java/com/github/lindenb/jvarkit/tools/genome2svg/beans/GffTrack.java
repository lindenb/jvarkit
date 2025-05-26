package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.iterator.LineIterators;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.TabixIteratorLineReader;
import htsjdk.tribble.readers.TabixReader;

public class GffTrack extends Track {

	private String tabix;
	public GffTrack() {
		
		}
	public void setTabix(String tabix) {
		this.tabix = tabix;
		}
	public String getTabix() {
		return tabix;
		}
	
	
	private void recursive(Gff3Feature feat,Pileup<Gff3Feature> pileup) {
		if(feat.getType().equals("gene")) {
			for(Gff3Feature c:feat.getChildren()) {
				pileup.add(c);
				}
			}
		else
			{
			for(Gff3Feature c:feat.getChildren()) {
				recursive(c,pileup);
				}
			}
		}
	
	@Override
	public void paint(SVGContext ctx) {
		ctx.styleNode.appendChild(ctx.text("."+getId()+"strand {fill:none;stroke:gray;stroke-width:0.5px;}\n"));
		ctx.styleNode.appendChild(ctx.text("."+getId()+"kgtr {fill:none;stroke:gray;stroke-width:1px;}\n"));
		ctx.styleNode.appendChild(ctx.text("."+getId()+"kgexon {fill:lightgray;stroke:gray;stroke-width:1px;}\n"));

		

		Element polyline = ctx.element("polyline");
		ctx.defsNode.appendChild(polyline);
		polyline.setAttribute("id", "strandF");
		polyline.setAttribute("class", getId()+"strand");
		polyline.setAttribute("points", "-5,-5 0,0 -5,5");
		
		polyline = ctx.element("polyline");
		ctx.defsNode.appendChild(polyline);
		polyline.setAttribute("id", "strandR");
		polyline.setAttribute("class", getId()+"strand");
		polyline.setAttribute("points", "5,-5 0,0 5,5");
		
		
		final Element g0 = ctx.element("g");
		ctx.tracksNode.appendChild(g0);
		insertTitle(g0,ctx);
		
		int chromStart = ctx.loc.getStart();
		int chromEnd = ctx.loc.getEnd();
		final Pileup<Gff3Feature> pileup = new Pileup<>(ctx.createCollisionPredicate());

		try(TabixReader tbx = new TabixReader(getTabix())) {
			Gff3Codec codec = new Gff3Codec(DecodeDepth.SHALLOW);
			TabixReader.Iterator iter = tbx.query(ctx.loc.getContig(), ctx.loc.getStart(), ctx.loc.getEnd());
			for(;;) {
				final String line = iter.next();
				if(line==null) break;
				final Gff3Feature gff3 = codec.decode(LineIterators.of(Collections.singletonList(line)));
				if(gff3==null) continue;
				if(!(gff3.getType().equals("gene"))) continue;
				chromStart = Math.min(chromStart, gff3.getStart());
				chromEnd = Math.max(chromEnd, gff3.getEnd());
				}
			codec = new Gff3Codec(DecodeDepth.DEEP);
			iter = tbx.query(ctx.loc.getContig(),chromStart, chromEnd);
			try(TabixIteratorLineReader lr = new TabixIteratorLineReader(iter)) {
				final LineIteratorImpl li = new LineIteratorImpl(lr);
				while(!codec.isDone(li)) {
					final Gff3Feature  feat = codec.decode(li);
					if(feat==null) continue;
					recursive(feat,pileup);
					}
				codec.close(li);
				}
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		
		
		final double exonHeight = 10;
		
		for(List<Gff3Feature> row:pileup.getRows()) {
			double midY= ctx.y + exonHeight/2.0;
			for(Gff3Feature tr:row) {
				final Element g = ctx.element("g");
				g0.appendChild(g);
				
				/* transcript line */
				Element line = ctx.element("line");
				g.appendChild(line);
				line.setAttribute("class",getId()+"kgtr");
				line.setAttribute("x1",String.valueOf(ctx.pos2pixel(ctx.trimpos(tr.getStart()))));
				line.setAttribute("y1",String.valueOf(midY));
				line.setAttribute("x2",String.valueOf(ctx.pos2pixel(ctx.trimpos(tr.getEnd()))));
				line.setAttribute("y2",String.valueOf(midY));
				

				
				/* strand symbols */
				for(double pixX=0;
					pixX< ctx.image_width;
					pixX+=30)
					{
					double pos1= ctx.pixel2genomic(pixX);
					if(pos1< tr.getStart()) continue;
					if(pos1> tr.getEnd()) break;
					Element use = ctx.element("use");
					g.appendChild(use);
					use.setAttribute("href", "#strand"+(tr.getStrand()==Strand.POSITIVE?"F":"R"));
					use.setAttribute("x",String.valueOf(pixX));
					use.setAttribute("y",String.valueOf(midY));
					}
				
				/* exons */
				for(final Gff3Feature exon: tr.getChildren())
					{
					if(!exon.getType().equals("exon")) continue;
					if(!exon.overlaps(ctx.loc)) continue;
					Element rect = ctx.element("rect");
					g.appendChild(rect);
					rect.setAttribute("class",getId()+"kgexon");
					
					rect.setAttribute("x",String.valueOf(ctx.pos2pixel(ctx.trimpos(exon.getStart()))));
					rect.setAttribute("y",String.valueOf(midY-exonHeight/2));
					rect.setAttribute("width",String.valueOf(
							ctx.pos2pixel(ctx.trimpos(exon.getEnd()))-ctx.pos2pixel(ctx.trimpos(exon.getStart()))));
					rect.setAttribute("height",String.valueOf(exonHeight));
					rect.appendChild(ctx.title(exon.getName()));
					}
				
				}
			
			// legend
			{
			final Set<String> geneNames = row.stream().flatMap(KG->KG.getAttribute("gene_name").stream()).collect(Collectors.toCollection(LinkedHashSet::new));
			if(geneNames.size()<10) {
				final Element legend = ctx.element("text",String.join(" ", geneNames));
				legend.setAttribute("style", "text-anchor:end;font-size:10px;");
				legend.setAttribute("x", format(-5));
				legend.setAttribute("y",String.valueOf(ctx.y+10));
				g0.appendChild(legend);
				}
			}

			
			ctx.y += exonHeight;
			ctx.y += 2;
			}
		
	}
}