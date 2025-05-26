package com.github.lindenb.jvarkit.gtf;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.GenomeSVGDocument;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Maps;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;

public class GtfTrack {
	private static final Logger LOG = Logger.of(GtfTrack.class);

	private Path gtfPath;
	private GenomeSVGDocument svgDoc;
	private final int max_rows=15;
	
	private static class Gene implements ExtendedLocatable {
		final List<Transcript> transcripts;
		Gene(List<Transcript> transcripts) {
			this.transcripts = transcripts;
			Collections.sort(this.transcripts, LocatableUtils.DEFAULT_COMPARATOR);
			}
		@Override
		public String getContig() {
			return transcripts.get(0).getContig();
			}
		@Override
		public int getStart() {
			return transcripts.get(0).getStart();
			}
		@Override
		public int getEnd() {
			return transcripts.stream().mapToInt(TR->TR.getEnd()).max().getAsInt();
			}
		String getGeneName() {
			return transcripts.get(0).getGeneName();
			}
		}

	
	private static class Transcript implements ExtendedLocatable {
		GTFLine transcript = null;
		int cdsStart=-1;
		int cdsEnd=-1;
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
		String getGeneId() {
			return transcript.getGeneId();
			}
		boolean isCoding() { return cdsStart>0 && cdsEnd>cdsStart;}
		String getGeneName() {
			String s= transcript.getGeneName();
			if(StringUtils.isBlank(s)) s=transcript.getTranscriptId();
			return s;
			}
		}
	
	private static class Params {
		List<Transcript> transcripts;
		String strandClass;
		String ktrClass;
		String exonClass;
		String cdsClass;
		String strandFId;
		String strandRId;
		double max_y;
		}
	
	public GtfTrack() {
		
		}
	public void setGtfPath(Path gtfPath) {
		this.gtfPath = gtfPath;
		}
	public void setDocument(GenomeSVGDocument svgDoc) {
		this.svgDoc = svgDoc;
		}
	
	
	private List<Transcript> fetch(final Locatable loc) {
		final GTFCodec codec = new GTFCodec();
		try(TabixReader tbx = new TabixReader(this.gtfPath.toString())) {
			final AutoMap<String,Transcript,Transcript> id2transcript=AutoMap.make(()->new Transcript());
			TabixReader.Iterator iter = tbx.query(loc.getContig(), loc.getStart(), loc.getEnd());
			for(;;) {
				final String line = iter.next();
				if(line==null) break;
				final GTFLine gtf = codec.decode(line);
				if(gtf==null) continue;
				if(!gtf.hasAttribute("transcript_id")) continue;
				final String transcript_id = gtf.getTranscriptId();
				if(StringUtils.isBlank(transcript_id)) continue;
				if(!(gtf.isExon() || gtf.isTranscript() || gtf.isCDS())) {
					continue;
					}
				final Transcript tr= id2transcript.insert(transcript_id);
				if(gtf.isTranscript()) {
					tr.transcript = gtf;
					}
				else if(gtf.isExon()) {
					tr.exons.add(gtf);
					}
				else if(gtf.isCDS()) {
					tr.cdsStart= tr.cdsStart<0?gtf.getStart():Math.min(tr.cdsStart, gtf.getStart());
					tr.cdsEnd= tr.cdsEnd<0?gtf.getEnd():Math.max(tr.cdsEnd, gtf.getEnd());
					}
				}
			return id2transcript.values().
				stream().
				filter(T->T.transcript!=null && !T.exons.isEmpty()).
				sorted(LocatableUtils.DEFAULT_COMPARATOR).
				collect(Collectors.toList());
			}
		catch(IOException err) {
			LOG.error(err);
			throw new RuntimeIOException(err);
			}
		}
	
	private String toString(GTFLine feat) {
		return feat.toNiceString();
	}
	
	private boolean painDense(final Params params,GenomeSVGDocument.IntervalInfo ii, double ii_y) {
		final AutoMap<String,Transcript,List<Transcript>> gene2tr = AutoMap.makeList();
		for(Transcript tr:params.transcripts) {
			String gn = tr.getGeneId();
			if(StringUtils.isBlank(gn)) continue;
			gene2tr.insert(gn,tr);
			}		
		final List<Locatable> merged_genes = LocatableUtils.mergeIntervals(
				gene2tr.entrySet().stream().
					map(KV->new Gene(KV.getValue())).
					map(LocatableUtils::reduce).
					collect(Collectors.toList()),
				ii.createCollisionPredicate(1).negate()
				);
		
		ii_y += 2;
		final Element g0 = this.svgDoc.group();
		final double featureHeight = 12;
		final double exonHeight = featureHeight;
		final double midY= ii_y + exonHeight/2.0;
		for(Locatable gene:merged_genes) {
			g0.appendChild(svgDoc.anchor( ii.rect(
					gene,
					midY-exonHeight/2,
					exonHeight,
					Maps.of("class",params.exonClass)
					),gene));
			}

		ii_y += exonHeight;
		ii_y += 2;
		params.max_y = Math.max(ii_y, params.max_y);
		svgDoc.rootElement.appendChild(g0);
		return true;
		}

	
	private boolean painPacked(final Params params,GenomeSVGDocument.IntervalInfo ii, double ii_y) {
		final AutoMap<String,Transcript,List<Transcript>> gene2tr = AutoMap.makeList();
		for(Transcript tr:params.transcripts) {
			String gn = tr.getGeneId();
			if(StringUtils.isBlank(gn)) continue;
			gene2tr.insert(gn,tr);
			}
		final Pileup<Gene> pileup = new Pileup<>(ii.createCollisionPredicate());
		pileup.addAll(gene2tr.entrySet().stream().map(KV->new Gene(KV.getValue())).
				sorted(LocatableUtils.DEFAULT_COMPARATOR).
				collect(Collectors.toList()));


		if(pileup.getRowCount() > this.max_rows) {
			LOG.info("too many rows for "+ii+" painPacked.");
			return false;
			}
		final Element g0 = this.svgDoc.group();
		final double featureHeight = 12;
		final double exonHeight = featureHeight;

		for(List<Gene> row:pileup.getRows()) {
			double midY= ii_y + exonHeight/2.0;
			
			for(Gene gene:row) {
				final Element g = svgDoc.group();
				g0.appendChild(svgDoc.anchor(svgDoc.setTitle(g,gene.getGeneName()),gene));
				final List<Locatable> merged_exons = LocatableUtils.mergeIntervals(
						gene.transcripts.stream().flatMap(TR->TR.exons.stream()).collect(Collectors.toList()),
						ii.createCollisionPredicate().negate()
						);

				
				if(merged_exons.size()>1) {
					/* transcript line */
					g.appendChild(ii.line(
							gene,
							midY,
							Maps.of("class",params.ktrClass)
							));
					
					
					/* strand symbols */
					for(double pixX=0;
						pixX<  ii.getWidth();
						pixX+=30)
						{
						final double pos1= ii.pixel2genomic(pixX+ii.getX());
						if(pos1< gene.getStart()) continue;
						if(pos1> gene.getEnd()) break;
						if(merged_exons.stream().anyMatch(EX->EX.getStart()<=pos1 || pos1<=EX.getEnd())) continue;
						g.appendChild(svgDoc.use(pixX, midY,"#"+(gene.transcripts.get(0).transcript.isPostiveStrand()?params.strandFId:params.strandRId)));
						}
					}
				
				/* exons */
				for(final Locatable exon: merged_exons)
					{
					if(!exon.overlaps(ii)) continue;
					g.appendChild(ii.rect(
						exon,
						midY-exonHeight/2,
						exonHeight,
						Maps.of("class",params.exonClass)
						));
					}
				
				} // end for transcript
			
			
			ii_y += exonHeight;
			ii_y += 2;
			params.max_y = Math.max(ii_y, params.max_y);
			} // end of rows
		svgDoc.rootElement.appendChild(g0);
		return true;
		}

	
	private boolean paintFull(final Params params,GenomeSVGDocument.IntervalInfo ii, double ii_y) {
		final Element g0 = this.svgDoc.group();
		final AutoMap<Integer,String,Set<String>> y2legend = AutoMap.makeSet();		
		final Pileup<Transcript> pileup = new Pileup<>(ii.createCollisionPredicate());
		pileup.addAll(params.transcripts);
		
		if(pileup.getRowCount() > this.max_rows) {
			LOG.info("too many rows for "+ii+" full.");
			return false;
			}
		
		final double featureHeight = 12;
		final double cdsHeight = featureHeight-1;
		final double exonHeight = cdsHeight/2;
	
		for(List<Transcript> row:pileup.getRows()) {
			double midY= ii_y + exonHeight/2.0;
			
			for(Transcript tr:row) {
				final Element g = svgDoc.group();
				g0.appendChild(svgDoc.anchor(svgDoc.setTitle(g,tr.getGeneName()+" "+tr.transcript.getTranscriptId()),tr));
				
				
				g.appendChild(svgDoc.line(
						ii.pos2pixel(ii.trimPos(tr.getStart())),
						midY,
						ii.pos2pixel(ii.trimPos(tr.getEnd()+1)),
						midY,
						Maps.of("class",params.ktrClass)
						));
				
				/* strand symbols */
				for(double pixX=0;
					pixX<  ii.getWidth();
					pixX+=30)
					{
					final double pos1= ii.pixel2genomic(pixX+ii.getX());
					if(pos1< tr.getStart()) continue;
					if(pos1> tr.getEnd()) break;
					if(tr.exons.stream().anyMatch(EX->EX.contains((int)pos1))) continue;
					g.appendChild(svgDoc.use(pixX+ii.getX(), midY,"#"+(tr.transcript.isPostiveStrand()?params.strandFId:params.strandRId)));
					}
				
				/* exons */
				for(final GTFLine exon: tr.exons)
					{
					if(!exon.overlaps(ii)) continue;
					// don't plot exon if CDS==exon
					if(tr.isCoding() && CoordMath.encloses(tr.cdsStart,tr.cdsEnd,exon.getStart(),exon.getEnd())) {
						continue;
						}
					g.appendChild(ii.rect(
						exon,
						midY-exonHeight/2,
						exonHeight,
						Maps.of("class",params.exonClass)
						));
					}
				if(tr.cdsStart>0 && tr.cdsEnd>0) {
					for(final GTFLine exon: tr.exons)
						{
						if(!exon.overlaps(ii)) continue;
						if(exon.getEnd() < tr.cdsStart) continue;
						if(exon.getStart() > tr.cdsEnd) continue;
						g.appendChild(ii.rect(
							new SimpleInterval(
									exon.getContig(),
									Math.max(tr.cdsStart, exon.getStart()),
									Math.min(tr.cdsEnd, exon.getEnd())
									),
							midY-cdsHeight/2,
							cdsHeight,
							Maps.of("class",params.cdsClass)
							));
						}
					}

				} // end for transcript
			final Set<String> geneNames = row.stream().
					map(KG->KG.transcript.getGeneName()).
					collect(Collectors.toCollection(LinkedHashSet::new));
			y2legend.insert((int)ii_y).addAll(geneNames);
			
			ii_y += featureHeight;
			ii_y += 2;
			params.max_y = Math.max(ii_y, params.max_y);
			} // end of rows
		
		
		svgDoc.rootElement.appendChild(g0);
		
		// legend
		/*
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
			}*/
		return true;
		}
	
	public void paint() {
		LOG.info("BEGIN GTF");
		if(this.gtfPath==null) return;
		if(this.svgDoc==null) return;
		final Params params = new Params();
		params.max_y = svgDoc.lastY;
		params.strandClass = svgDoc.style2class("fill:none;stroke:gray;stroke-width:0.5px;");
		params.ktrClass = svgDoc.style2class("fill:none;stroke:gray;stroke-width:1px;");
		params.exonClass = svgDoc.style2class("fill:lightgray;stroke:gray;stroke-width:1px;");
		params.cdsClass = svgDoc.style2class("fill:plum;stroke:gray;stroke-width:1px;");
		LOG.info("BEGIN GTF");

		

		params.strandFId= svgDoc.insertDefElement(
				svgDoc.makeNodeBuilder("polyline").
					attribute("class", params.strandClass).
					attribute("points","-5,-5 0,0 -5,5").
					makeElement()
				);
		
		params.strandRId= svgDoc.insertDefElement(
				svgDoc.makeNodeBuilder("polyline").
					attribute("class", params.strandClass).
					attribute("points","5,-5 0,0 5,5").
					makeElement()
				);
		
		svgDoc.addBanner("Genes ("+this.gtfPath.getFileName().toString()+")", 12);
		svgDoc.lastY+=14;
		double y = svgDoc.lastY;

		
		for(GenomeSVGDocument.IntervalInfo ii : this.svgDoc.getIntervalInfoList()) {
			params.transcripts = fetch(ii);
			if(!paintFull(params,ii,y)) {
				if(!painPacked(params,ii,y)) {
					painDense(params,ii,y);
					}
				}
			}
		svgDoc.frame(svgDoc.lastY,params.max_y);
		svgDoc.lastY = params.max_y;
		
		
		
	
		LOG.info("end gtf");
		}
}