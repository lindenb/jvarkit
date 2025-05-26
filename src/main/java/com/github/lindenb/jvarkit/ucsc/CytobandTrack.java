package com.github.lindenb.jvarkit.ucsc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.samtools.util.LocatableDelegate;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.svg.GenomeSVGDocument;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.RuntimeIOException;

public class CytobandTrack {
	private static final Logger LOG = Logger.of(CytobandTrack.class);
	private Path cytobandPath = null;
	private GenomeSVGDocument svgDoc = null;
	final double featureHeight=20;

	
	private static class CytoInfo extends LocatableDelegate<Cytoband> {
		double x;
		double width;
		CytoInfo(Cytoband c) {
			super(c);
			}
		boolean isCentromer() { return getDelegate().getStain().contains("acen");}
		double getMaxX() { return x+width;}
		}
	
	private static class Chromosome implements ExtendedLocatable,Iterable<CytoInfo> {
		List<CytoInfo> bands = new ArrayList<>();
		@Override
		public String getContig() {
			return bands.get(0).getContig();
			}
		@Override
		public int getStart() {
			return bands.get(0).getStart();
			}
		@Override
		public int getEnd() {
			return bands.get(bands.size()-1).getEnd();
			}
		@Override
		public Iterator<CytoInfo> iterator() {
			return bands.iterator();
			}
		void sort() {
			Collections.sort(bands,LocatableUtils.DEFAULT_COMPARATOR);
			}
		OptionalDouble getCentromereX() {
			List<CytoInfo> L = this.bands.stream().filter(C->C.isCentromer()).collect(Collectors.toList());
			if(L.isEmpty()) return OptionalDouble.empty();
			if(L.size()==1) return OptionalDouble.of(L.get(0).x+(L.get(0).width/2));
			if(L.size()==2) return OptionalDouble.of(L.get(0).getMaxX());
			return  OptionalDouble.empty();
			}
		double getX() { return bands.get(0).x;}
		double getMaxX() { return bands.get(bands.size()-1).getMaxX();}
		double getWidth() { return getMaxX()-getX();}
		
		double loc2pix(int pos) {
			return getX()+(pos/(double)getLengthOnReference())*getWidth();
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
	
	


	
	public void paint() {
		if(this.svgDoc==null) return;
		final AutoMap<String,Chromosome,Chromosome> entries= AutoMap.make(()->new Chromosome());
		final SAMSequenceDictionary dict = this.svgDoc.getSAMSequenceDictionary();
		final ContigDictComparator ctgCompare = new ContigDictComparator(dict);
		BufferedReader br = null;
		try {
			if(this.cytobandPath==null) {
				final Optional<String> opt = Cytoband.getURLForBuild(dict);
				if(!opt.isPresent()) return;
				LOG.info("download cytobands from "+opt.get());
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
				entries.insert(ctg).bands.add(e);
				}
			entries.values().stream().forEach(C->C.sort());
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		finally
			{
			if(br!=null) try {br.close();} catch(Throwable err) {}
			}

		
		
		
		if(entries.isEmpty()) {
			LOG.warning("no cytoband was found");
			return ;
			}
		double spaceBetweenInterval = 5;
		
		svgDoc.addBanner("Cytobands", 10);
		
		double y = svgDoc.lastY;
		y+=2;
		
		final double sum_length_on_ref = entries.values().stream().
				mapToLong(T->T.getLengthOnReference()).
				sum();
		final List<Chromosome> orderedChromosomes = entries.
				values().
				stream().
				sorted(ctgCompare.createLocatableComparator()).
				collect(Collectors.toList());
		final double adjust_image_width = this.svgDoc.image_width - Math.max(0,(entries.size()-1)*spaceBetweenInterval);
		
		final Element g = svgDoc.group();
		g.appendChild(svgDoc.comment("BEGIN cytoband"));
		
		final Map<String,String> stain2gradientid=new HashMap<>();
		
		double x = this.svgDoc.margin_left;
		for(Chromosome chromosome : orderedChromosomes) {
			final Element gContig = svgDoc.group();
			g.appendChild(gContig);
			
			
			
			final double contigX1 = x;

			
			for(CytoInfo ii: chromosome.bands) {
				ii.x=x;
				ii.width = (ii.getLengthOnReference()/sum_length_on_ref)* adjust_image_width;
				x+= ii.width;
				}
			
			final String clipId = svgDoc.nextId();
			final Element clipPath = svgDoc.clipPath();
			svgDoc.defsElement.appendChild(clipPath);
			clipPath.setAttribute("id", clipId);

			final double roundRec = featureHeight/2;
			final OptionalDouble centromere = chromosome.getCentromereX();
			List<Element> at_the_end_clone_and_stroke = new ArrayList<>();
			
			if(centromere.isPresent()) {
				Element r= svgDoc.rect(
						chromosome.getX(),
						y+featureHeight,
						centromere.getAsDouble() - chromosome.getX(),
						featureHeight,
						Maps.of("rx",roundRec)
						);
				at_the_end_clone_and_stroke.add(r);
				clipPath.appendChild(r);
				
				r = svgDoc.rect(
						centromere.getAsDouble(),
						y+featureHeight,
						chromosome.getMaxX()-centromere.getAsDouble(),
						featureHeight,
						Maps.of("rx",roundRec)
						);
				clipPath.appendChild(r);
				at_the_end_clone_and_stroke.add(r);
				}
			else
				{
				final Element r= svgDoc.rect(
						chromosome.getX(),
						y+featureHeight,
						chromosome.getWidth(),
						featureHeight,
						Maps.of("rx",roundRec)
						);
				clipPath.appendChild(r);
				at_the_end_clone_and_stroke.add(r);
				}
			
			for(CytoInfo ii: chromosome.bands) {
				String gradientid = stain2gradientid.get(ii.getDelegate().getStain());
				if(gradientid==null) {
					gradientid = svgDoc.createVerticalLinearGradient("black", ii.getDelegate().getCssColor());
					stain2gradientid.put(ii.getDelegate().getStain(), gradientid);
					}
				
				
				final Element rec = svgDoc.rect(
						ii.x,y+featureHeight, ii.width,featureHeight,
						Maps.of("fill", "url(#"+gradientid+")")
						);
				if(clipId!=null) rec.setAttribute("clip-path", "url(#"+clipId+")");
				svgDoc.setTitle(rec,ii.getContig()+ii.getDelegate().getName());
				gContig.appendChild(svgDoc.anchor(rec, ii));
				}
			
			// create frame over clipped elements
			for(Element rec0: at_the_end_clone_and_stroke) {
				final Element rec= Element.class.cast(rec0.cloneNode(false));
				rec.setAttribute("class",svgDoc.style2class("stroke:black;fill:none;"));
				gContig.appendChild(rec);
				}
			
			
			gContig.appendChild(svgDoc.text(
					(contigX1+x)/2.0,
					y+featureHeight,
					chromosome.getContig(),
					Maps.of("text-align", "center","fill","black"))
					);
			
			for(GenomeSVGDocument.IntervalInfo info: this.svgDoc.getIntervalInfoList()) {
				if(!info.getContig().equals(chromosome.getContig())) continue;
				final double x1= chromosome.loc2pix(info.getStart());
				final double x2= chromosome.loc2pix(info.getEnd());
				
				final Element r = this.svgDoc.rect(
					x1,
					y+featureHeight-2,
					Math.max(1,x2-x1),
					featureHeight+4,
					Maps.of("stroke","red","fill","red","opacity",0.6)
					);
				svgDoc.setTitle(r,info.getName());
				gContig.appendChild(svgDoc.anchor(r, info));	
				}
			
			x+=spaceBetweenInterval;
			}
		svgDoc.rootElement.appendChild(g);
		g.appendChild(svgDoc.comment("END cytoband"));
		
		y+= featureHeight*2+3;
		svgDoc.lastY =y;
		}


}
