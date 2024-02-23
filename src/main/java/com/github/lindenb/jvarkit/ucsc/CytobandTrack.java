package com.github.lindenb.jvarkit.ucsc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.LocatableDelegate;
import com.github.lindenb.jvarkit.svg.GenomeSVGDocument;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.RuntimeIOException;

public class CytobandTrack {
	private Path cytobandPath = null;
	private GenomeSVGDocument svgDoc = null;
	final double featureHeight=20;

	
	private static class CytoInfo extends LocatableDelegate<Cytoband> {
		double x;
		double width;
		CytoInfo(Cytoband c) {
			super(c);
			}
		}
	
	public CytobandTrack() {
		}
	
	public void setCytobandFilePath(final Path cytobandPath) {
		this.cytobandPath = cytobandPath;
		}
	
	public void setDocument(GenomeSVGDocument svgDoc) {
		this.svgDoc = svgDoc;
		}
	
	private double loc2pix(final String contig,final int pos,List<CytoInfo> entries) {
		CytoInfo c = entries.stream().filter(C->C.getContig().equals(contig) && C.contains(pos)).findFirst().orElse(null);
		if(c==null) throw new IllegalArgumentException();
		return c.x + ((pos-c.getStart())/(double)c.getLengthOnReference())* c.width;
		}

	
	public void paint() {
		final AutoMap<String,CytoInfo,List<CytoInfo>> entries= AutoMap.makeList();
		final SAMSequenceDictionary dict = this.svgDoc.getSAMSequenceDictionary();
		final ContigDictComparator ctgCompare = new ContigDictComparator(dict);
		BufferedReader br = null;
		try {
			if(this.cytobandPath==null) {
				final Optional<String> opt = Cytoband.getURLForBuild(dict);
				if(!opt.isPresent()) return;
				br = IOUtils.openURIForBufferedReading(opt.get());
				}
			else
				{
				br = IOUtils.openPathForBufferedReading(this.cytobandPath);
				}
			final UnaryOperator<String> convert= ContigNameConverter.fromOneDictionary(dict);
			for(;;) {
				final String line = br.readLine();
				if(line==null) break;
				final String[] tokens = CharSplitter.TAB.split(line);
				final String ctg = convert.apply(tokens[0]);
				if(StringUtils.isBlank(ctg)) continue;
				tokens[0] = ctg;
				final CytoInfo e = new CytoInfo(Cytoband.of(tokens));
				if(this.svgDoc.getIntervalInfoList().stream().noneMatch(C->C.contigsMatch(e))) continue;
				entries.insert(ctg, e);
				}
			for(List<CytoInfo> L: entries.values()) {
				Collections.sort(L, ctgCompare.createLocatableComparator());
				}
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		finally
			{
			if(br!=null) try {br.close();} catch(Throwable err) {}
			}

		
		
		
		if(entries.isEmpty()) return ;
		double spaceBetweenInterval = 1;
		double y = svgDoc.lastY;
		
		svgDoc.appendStyle(
			".ctgborderp {fill:url(#grad01);stroke:green;}" +
			".ctgborderq {fill:url(#grad01);stroke:green;}" +
			".ctglabel {text-anchor:middle;stroke:none;fill:darkgrey;font: bold 10px Verdana, Helvetica, Arial, sans-serif;}" +
			".cytoband {fill:silver;stroke:none;}" +
			".bedlabel {stroke:red;fill:none;text-anchor:start;font: normal 7px Verdana, Helvetica, Arial, sans-serif;}" +
			".maintitle {stroke:none;fill:darkgrey;text-anchor:middle;font: normal 12px Verdana, Helvetica, Arial, sans-serif;}" +
			".ctgback {fill:gainsboro;stroke:none;filter:url(#filter01);}"
			);
		
		final double sum_length_on_ref = entries.values().stream().
				flatMap(T->T.stream()).
				mapToLong(T->T.getLengthOnReference()).
				sum();
		final List<String> orderedChromosomes = entries.keySet().stream().sorted(ctgCompare).collect(Collectors.toList());
		final double adjust_image_width = this.svgDoc.image_width - Math.max(0,(entries.size()-1)*spaceBetweenInterval);
		
		double x = this.svgDoc.margin_left;
		for(int i=0;i< orderedChromosomes.size();++i) {
			final String ctg =orderedChromosomes.get(i);
			
			final Element gContig = svgDoc.group();
			
			if(i>0) x+=spaceBetweenInterval;
			for(CytoInfo ii: entries.get(ctg)) {
				ii.x=x;
				ii.width = (ii.getLengthOnReference()/(double)sum_length_on_ref)* adjust_image_width;
				x+= ii.width;
				
				final Element rec = svgDoc.rect(
						ii.x,y, ii.width,featureHeight,
						Maps.of("class", svgDoc.style2class("stroke:none;fill:"+ii.getDelegate().getCssColor()))
						);
				svgDoc.setTitle(rec,ii.getDelegate().getName());
				gContig.appendChild(svgDoc.anchor(rec, ii));
				}
			for(GenomeSVGDocument.IntervalInfo info: this.svgDoc.getIntervalInfoList()) {
				if(!info.getContig().equals(ctg)) continue;
				double x1= loc2pix(info.getContig(),info.getStart(),entries.get(ctg));
				double x2= loc2pix(info.getContig(),info.getEnd()+1,entries.get(ctg));
				
				Element polygon = svgDoc.createPolygon().
					lineTo(x1,y+featureHeight).
					lineTo(x2,y+featureHeight).
					lineTo(info.getX()+info.getWidth(),y+featureHeight*2).
					lineTo(info.getX(),y+featureHeight*2).
					make(Maps.of("stroke","red","fill","red"));
				svgDoc.setTitle(polygon,info.toNiceString());
				gContig.appendChild(svgDoc.anchor(polygon, info));	
				}
			
		
			}
		y+= featureHeight*2;
		svgDoc.lastY =y;
		}


}
