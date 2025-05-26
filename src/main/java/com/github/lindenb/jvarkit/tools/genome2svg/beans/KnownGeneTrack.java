package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.tabix.TabixFileReader;
import com.github.lindenb.jvarkit.tools.genome2svg.SVGContext;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;

import htsjdk.samtools.util.RuntimeIOException;

public class KnownGeneTrack extends Track {

	private PathBean tabix;
	public KnownGeneTrack() {
		
		}
	public void setTabix(PathBean tabix) {
		this.tabix = tabix;
		}
	public PathBean getTabix() {
		return tabix;
		}
	@Override
	public void paint(SVGContext ctx) {
		ctx.styleNode.appendChild(ctx.text("."+getId()+"strand {fill:none;stroke:gray;stroke-width:0.5px;}\n"));
		ctx.styleNode.appendChild(ctx.text("."+getId()+"kgtr {fill:none;stroke:gray;stroke-width:1px;}\n"));
		ctx.styleNode.appendChild(ctx.text("."+getId()+"kgexon {fill:lightgray;stroke:gray;stroke-width:1px;}\n"));

		
		final Element g0 = ctx.element("g");
		ctx.tracksNode.appendChild(g0);
		insertTitle(g0,ctx);
		try(TabixFileReader tbx = new TabixFileReader(getTabix().getPath())) {
			Iterator<String> iter= tbx.iterator(
					ctx.loc.getContig(),
					ctx.loc.getStart(),
					ctx.loc.getEnd()
					);
			final UcscTranscriptCodec codec = new UcscTranscriptCodec();
			final Pileup<UcscTranscript> pileup = new Pileup<>(ctx.createCollisionPredicate());
			
			while(iter.hasNext()) {
				final UcscTranscript tr=codec.decode(iter.next());
				pileup.add(tr);
				}
			if(!pileup.isEmpty()) {
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
				}
			final double exonHeight = 10;
			
			for(List<UcscTranscript> row:pileup.getRows()) {
				double midY= ctx.y + exonHeight/2.0;
				for(UcscTranscript tr:row) {
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
						use.setAttribute("href", "#strand"+(tr.isPositiveStrand()?"F":"R"));
						use.setAttribute("x",String.valueOf(pixX));
						use.setAttribute("y",String.valueOf(midY));
						}
					
					/* exons */
					for(final UcscTranscript.Exon exon: tr.getExons())
						{
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
				final Set<String> geneNames = row.stream().map(KG->StringUtils.ifBlank(KG.getName2(), KG.getTranscriptId())).collect(Collectors.toCollection(LinkedHashSet::new));
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
		catch(IOException err ) {
			throw new RuntimeIOException(err);
			}
		}
	
}
