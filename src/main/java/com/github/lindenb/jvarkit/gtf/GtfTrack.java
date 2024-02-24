package com.github.lindenb.jvarkit.gtf;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.GenomeSVGDocument;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;

public class GtfTrack {
	private Path gtfPath;
	private GenomeSVGDocument svgDoc;
	
	private static class Transcript implements Locatable {
		GTFLine transcript = null;
		final List<GTFLine> exons = new ArrayList<>();
		Transcript() {
			}
		@Override
		public String getContig() {
			return transcript.getContig();
			}
		@Override
		public int getStart() {
			return transcript.getStart();
			}
		@Override
		public int getEnd() {
			return transcript.getEnd();
			}
		}
	
	public GtfTrack() {
		
		}
	public void setGtfPath(Path gtfPath) {
		this.gtfPath = gtfPath;
		}
	public void setDocument(GenomeSVGDocument svgDoc) {
		this.svgDoc = svgDoc;
		}
	
	
	public void paint() {
		final String strandClass = svgDoc.style2class("strand", "fill:none;stroke:gray;stroke-width:0.5px;");
		final String ktrClass = svgDoc.style2class("strand", "fill:none;stroke:gray;stroke-width:1px;");
		final String exonClass = svgDoc.style2class("strand", "fill:lightgray;stroke:gray;stroke-width:1px;");
		
		final AutoMap<Integer,String,Set<String>> y2legend = AutoMap.makeSet();

		

		Element polyline = svgDoc.element("polyline");
		svgDoc.defsElement.appendChild(polyline);
		polyline.setAttribute("id", "strandF");
		polyline.setAttribute("class",strandClass);
		polyline.setAttribute("points", "-5,-5 0,0 -5,5");
		
		polyline = svgDoc.element("polyline");
		svgDoc.defsElement.appendChild(polyline);
		polyline.setAttribute("id", "strandR");
		polyline.setAttribute("class",strandClass);
		polyline.setAttribute("points", "5,-5 0,0 5,5");
		
		
		final Element g0 = svgDoc.group();
		svgDoc.rootElement.appendChild(g0);
		//insertTitle(g0,ctx);
		final GTFCodec codec = new GTFCodec();

		double max_y = svgDoc.lastY;
		for(GenomeSVGDocument.IntervalInfo ii : this.svgDoc.getIntervalInfoList()) {
			double ii_y = svgDoc.lastY;
			final Pileup<Transcript> pileup = new Pileup<>(ii.createCollisionPredicate());
		
			try(TabixReader tbx = new TabixReader(this.gtfPath.toString())) {
				final Map<String,Transcript> id2transcript= new HashMap<>();
				TabixReader.Iterator iter = tbx.query(ii.getContig(), ii.getStart(), ii.getEnd());
				for(;;) {
					final String line = iter.next();
					if(line==null) break;
					final GTFLine gtf = codec.decode(line);
					if(gtf==null) continue;
					if(!gtf.hasAttribute("transcript_id")) continue;
					final String transcript_id = gtf.getAttribute("transcript_id");
					if(StringUtils.isBlank(transcript_id)) continue;
					if(!(gtf.getType().equals("transcript") || gtf.getType().equals("exon"))) {
						continue;
						}
					Transcript tr = id2transcript.get(transcript_id);
					if(tr==null) {
						tr =new Transcript();
						id2transcript.put(transcript_id,tr);
						}
					if(gtf.getType().equals("transcript")) {
						tr.transcript = gtf;
						}
					else if(gtf.getType().equals("exon")) {
						tr.exons.add(gtf);
						}
					}
				id2transcript.values().
					stream().
					filter(T->T.transcript!=null && !T.exons.isEmpty()).
					sorted(LocatableUtils.DEFAULT_COMPARATOR).
					forEach(T->pileup.add(T));
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			
			
		final double featureHeight = 12;
		final double exonHeight = featureHeight;
		
		for(List<Transcript> row:pileup.getRows()) {
			double midY= ii_y + exonHeight/2.0;
			
			for(Transcript tr:row) {
				final Element g = svgDoc.group();
				g0.appendChild(g);
				
				/* transcript line */
				g.appendChild(svgDoc.line(
						ii.pos2pixel(ii.trimPos(tr.getStart())),
						midY,
						ii.pos2pixel(ii.trimPos(tr.getEnd()+1)),
						midY,
						Maps.of("class",ktrClass)
						));
				
				
				/* strand symbols */
				for(double pixX=0;
					pixX<  ii.getWidth();
					pixX+=30)
					{
					double pos1= ii.pixel2genomic(pixX);
					if(pos1< tr.getStart()) continue;
					if(pos1> tr.getEnd()) break;
					g.appendChild(svgDoc.use(pixX, midY,"#strand"+(tr.transcript.isPostiveStrand()?"F":"R")));
					}
				
				/* exons */
				for(final GTFLine exon: tr.exons)
					{
					if(!exon.overlaps(ii)) continue;
					g.appendChild(ii.rect(
						exon,
						midY-exonHeight/2,
						exonHeight,
						Maps.of("class",exonClass)
						));
					}
				
				} // end for transcript
			final Set<String> geneNames = row.stream().
					map(KG->KG.transcript.getAttribute("gene_name")).
					collect(Collectors.toCollection(LinkedHashSet::new));
			y2legend.insert((int)ii_y).addAll(geneNames);
			
			ii_y += exonHeight;
			ii_y += 2;
			max_y = Math.max(ii_y, max_y);
			} // end loop rows
		svgDoc.lastY = max_y;
		} // end interval info
		
	// legend
	for(Integer y: y2legend.keySet())
		{
		if(y2legend.get(y).size()<10) {
			final Element legend = svgDoc.text(
					y+10,
					svgDoc.margin_left-5,
					String.join(" ", y2legend.get(y)),
					Maps.of("style", "text-anchor:end;font-size:10px;")
					);
			g0.appendChild(legend);
			}
		}
}
}