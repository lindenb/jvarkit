/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bed2svg;

import java.awt.Insets;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiPredicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.w3c.dom.Element;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.tools.bed2vcf.BedToVcf;
import com.github.lindenb.jvarkit.tools.misc.CytobandToSvg;
import com.github.lindenb.jvarkit.ucsc.Cytoband;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.vcf.MultiIntervalVariantIterator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

## Motivation

Convert BED file to SVG

## Example



```

END_DOC
*/
@Program(name="bed2svg",
description="Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele",
keywords={"bed","svg"},
creationDate="20261013",
modificationDate="20261013",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class BedToSVG extends Launcher {
	private static final Logger LOG = Logger.of(BedToSVG.class);
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"--cytoband"},description=Cytoband.OPT_DESC)
	private Path cytoBandFile=null;
	@Parameter(names={"-R"},description=DICTIONARY_SOURCE,required = true)
	private Path faidx=null;
	@Parameter(names={"--width"},description="window width in pixels")
	private int screen_width=1000;
	@Parameter(names={"--regions"},description="regions of interest, multiple comma separated.")
	private List<String> regions_str = new ArrayList<>();
	@Parameter(names={"--regions-file"},description="regions of interest, multiple comma separated.")
	private Path regions_file = null;
	@Parameter(names={"--min-feature-length"},description="discard  features where the length of reference is greater than x. " + DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter = NoSplitter.class)
	private int min_feature_length=0;
	@Parameter(names={"--max-feature-length"},description="discard  features where the length of reference is lower than x. "+ DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter = NoSplitter.class)
	private int max_feature_length=-1;

	private String ssr_regex = "";
	private long min_ssr_length = 0L;
	
	private final List<AreaOfInterest> areas_of_interest = new ArrayList<>();
	private SAMSequenceDictionary dictionary;
	
	
	private interface AreaOfInterest extends HasName {
		public boolean accept(Locatable loc);
		default BiPredicate<Locatable,Locatable> createNCollisionTest() {
			return (left,right)->{
				double x0 = genomic2pixel(left.getContig(), trim(left.getContig(),left.getEnd()));
				double x1 = genomic2pixel(right.getContig(), trim(right.getContig(),right.getStart()));
				return x0 +5 < x1;
				};
			}
		/** trim position so it fit in the area */
		int trim(String contig,int pos1);
		double genomic2pixel(String contig,int pos);
		Set<String> getContigs();
		double getX();
		double getWidth();
		}
	
	private class IntervalOfInterest extends SimpleInterval implements AreaOfInterest {
		double x;
		double width;
		IntervalOfInterest(Locatable loc) {
			super(loc);
			}
		@Override
		public String getName() {
			return toNiceString();
			}
		@Override
		public boolean accept(Locatable loc) {
			return this.overlaps(loc);
			}
		@Override
		public int trim(String contig,int pos1) {
			if(!this.getContig().equals(contig)) throw new IllegalArgumentException();
			return Math.min(Math.max(pos1, this.getStart()),this.getEnd());
			}
		@Override
		public double genomic2pixel(String contig, int pos) {
			return this.x+ ((pos-this.getStart())/(double)this.getLengthOnReference())*this.width;
			}
		@Override
		public BiPredicate<Locatable, Locatable> createNCollisionTest() {
			return null;
			}
		@Override
		public Set<String> getContigs() {
			return Collections.singleton(getContig());
			}
		@Override
		public double getX() {
			return x;
			}
		@Override
		public double getWidth() {
			return this.width;
			}
		}
	private class WholeGenomeArea  implements AreaOfInterest {
		final SAMSequenceDictionary dict;
		final Map<String, Long> contig2idx=new HashMap<>();
		final long genome_length;
		WholeGenomeArea(final SAMSequenceDictionary dict) {
			this.dict=dict;
			this.genome_length = dict.getReferenceLength();
			long n=0L;
			for(SAMSequenceRecord ssr:dict.getSequences()) {
				this.contig2idx.put(ssr.getContig(), n);
				n+=ssr.getLengthOnReference();
				}
			}
		@Override
		public String getName() {
			return "Whole Genome";
			}
		@Override
		public boolean accept(final Locatable loc) {
			return contig2idx.containsKey(loc.getContig());
			}
		@Override
		public int trim(String contig,int pos1) {
			final SAMSequenceRecord ssr = this.dict.getSequence(contig);
			if(ssr==null) throw new IllegalArgumentException();
			return Math.min(Math.max(pos1,1),ssr.getLengthOnReference());
			}
		@Override
		public double genomic2pixel(String contig,int pos) {
			final long ctg_idx = this.contig2idx.get(contig);
			return ((ctg_idx+pos)/(double)this.genome_length)* BedToSVG.this.screen_width;
			}
		@Override
		public Set<String> getContigs() {
			return Collections.unmodifiableSet(this.contig2idx.keySet());
			}
		@Override
		public double getX() {
			return 0;
			}
		@Override
		public double getWidth() {
			return BedToSVG.this.screen_width;
			}
		}
	
	private class Track {
		final List<Locatable> records= new ArrayList<>();
		String name="";
		double y=0;
		double max_y=0;
		}
	
	private final Map<String,Track> id2track = new LinkedHashMap<>();
	
	private Track createTrack(final String trackName) {
		Track track= id2track.get(trackName);
		if(track==null) {
			track = new  Track();
			track.name = trackName;
			this.id2track.put(trackName, track);
			}
		return track;
		}
	
	private Map<String,Object> getAttributesFor(Locatable loc)  {
		if(loc instanceof BedLine) {
			Map<String,Object > attributes =new HashMap<>();
			try {
			 final StreamTokenizer st=new StreamTokenizer(new StringReader(BedLine.class.cast(loc).getOrDefault(4,"")));
             st.wordChars('_', '_');
             String key=null;
             while(st.nextToken() != StreamTokenizer.TT_EOF)
                 {
                 String s=null;
                 switch(st.ttype)
                     {
                     case StreamTokenizer.TT_NUMBER: s=String.valueOf(st.nval);break;
                     case '"': case '\'' : case StreamTokenizer.TT_WORD: s=st.sval;break;
                     case ';':break;
                     default:break;
                     }
                 if(s==null) continue;
                 if(key==null)
                     {
                     key=s;
                     }
                 else
                     {
                     attributes.put(key,s);
                     key=null;
                     }
                 }
			} catch(IOException err) {
				LOG.warn(err);
				}
            if(!attributes.isEmpty()) return attributes;
			}
		else if(loc instanceof VariantContext) {
			String svtype=((VariantContext) loc).getAttributeAsString(VCFConstants.SVTYPE, "");
			String fill="gray";
			if(svtype.equals("DEL")) fill="red";
			else if(svtype.equals("DUP")) fill="green";
			else if(svtype.equals("INV")) fill="pink";
			return Map.of("stroke", "darkgray","fill",fill);
			}
		
		return Map.of("stroke", "pink","fill","green");
		}
	
	private boolean accept(final Locatable loc) {
		if(max_feature_length>=0 && loc.getLengthOnReference()> this.max_feature_length) return false;
		if(loc.getLengthOnReference()< this.min_feature_length) return false;
		if(this.areas_of_interest.stream().noneMatch(R->R.accept(loc))) return false;
		return true;
		}
	
	/** guess the file format from file suffix */
	private void readFileAnyFormat(final Path filePath) throws IOException {
		final String fname = filePath.getFileName().toString();
		if(fname.endsWith(FileExtensions.BED) || fname.endsWith(".bed.gz") || fname.endsWith(".bed.bgz")) {
			try(BufferedReader r = IOUtils.openPathForBufferedReading(filePath)) {
				readBed(r,IOUtils.getFilenameWithoutCommonSuffixes(filePath));
				}
			}
		else if(FileExtensions.VCF_LIST.stream().anyMatch(x->fname.endsWith(x))) {
			readVcf(filePath);
			}
		else
			{
			throw new IllegalArgumentException("cannot read data of type: "+filePath);
			}
	}
	
	/** read BED file */
	private void readBed(BufferedReader r,final String filename) throws IOException {
		try(BedLineReader blr = new BedLineReader(r,filename)) {
			blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(this.dictionary));
			while(blr.hasNext()) {
				final BedLine rec = blr.next();
				if(!accept(rec)) continue;
				final String trackName= (filename==null?"":filename);
				final Track track= createTrack(trackName);
				track.records.add(rec);
				}
			}
		}
	/** open file file, open index if needed/available */
	private void readVcf(final Path vcfFile)  throws IOException {
		boolean require_index= this.areas_of_interest.stream().allMatch(R->R instanceof IntervalOfInterest) && 
				vcfFile.getFileName().toString().endsWith(FileExtensions.VCF);
		
		try(VCFFileReader r= new VCFFileReader(vcfFile,require_index)) {
			if(require_index) {
				try(CloseableIterator<VariantContext> iter= MultiIntervalVariantIterator.query(r, this.areas_of_interest.stream().map(R->(IntervalOfInterest)R).collect(Collectors.toList()))) {
					readVcf(r.getFileHeader(),iter,vcfFile);
					}
				}
			else
				{
				try(CloseableIterator<VariantContext> iter= r.iterator()) {
					readVcf(r.getFileHeader(),iter,vcfFile);
					}
				}
			}
		}
	
	private void readVcf(final VCFHeader hdr,final CloseableIterator<VariantContext> iter,final Path filename)  throws IOException {
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			if(!ctx.isVariant()) continue;
			if(!accept(ctx)) continue;
			
			
			if(ctx.hasGenotypes()) {
				for(Genotype g: ctx.getGenotypes()) {
					if(!g.hasAltAllele()) continue;
					String trackName = g.getSampleName();
					final Track track= createTrack(trackName);
					track.records.add(
						new VariantContextBuilder(ctx)
							.noGenotypes()
							.genotypes(Collections.singletonList(g))
							.make()
						);
					}
				}
			else
				{
				String trackName= IOUtils.getFilenameWithoutCommonSuffixes(filename);
				final Track track= createTrack(trackName);
				track.records.add(ctx);
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final SAMSequenceDictionary original_dict=SequenceDictionaryUtils.extractRequired(this.faidx);
			final Hyperlink hyperlink = Hyperlink.compile(original_dict);
			this.dictionary = new SAMSequenceDictionary(
						original_dict
							.getSequences()
							.stream()
							.filter(SR->SR.getSequenceLength() >= this.min_ssr_length)
							.filter(SR->StringUtils.isBlank(ssr_regex) || SR.getContig().matches(ssr_regex))
							.collect(Collectors.toList())
						);
			if(this.dictionary.isEmpty()) {
				LOG.info("no chromosome to display");
				}
			final boolean splitted_intervals;
			if(!this.regions_str.isEmpty() || this.regions_file!=null)
				{
				splitted_intervals = true;
				/** load regions of interest */
				final List<SimpleInterval> regions=new ArrayList<>();
				final ContigNameConverter ctgConvert = ContigNameConverter.fromOneDictionary(this.dictionary);
				final IntervalParser parser= new IntervalParser(this.dictionary);
				this.regions_str
					.stream()
					.flatMap(S->Arrays.stream(CharSplitter.COMMA.split(S)))
					.filter(S->!StringUtils.isBlank(S))
					.map(S->parser.apply(S).orElse(null))
					.filter(R->R!=null)
					.forEach(R->regions.add(R))
					;
				
				if(this.regions_file!=null) {
					try(BedLineReader br=new BedLineReader(this.regions_file)) {
						br.setContigNameConverter(ctgConvert);
						br.stream().forEach(R->regions.add(new SimpleInterval(R)));
						}
					}
				final List<IntervalOfInterest> L= LocatableUtils.mergeIntervals(regions)
					.stream()
					.map(R->new IntervalOfInterest(R))
					.sorted(new ContigDictComparator(this.dictionary).createLocatableComparator())
					.collect(Collectors.toList());
					
				if(L.isEmpty()) {
					LOG.error("No region of interest contains a valid interval");
					return -1;
					}
				final double sum= L.stream().mapToInt(R->R.getLengthOnReference()).sum();
				double x=0;
				for(IntervalOfInterest r: L) {
					r.x = x;
					r.width = (r.getLengthOnReference()/sum)*this.screen_width;
					x+= r.width;
					}
				L.forEach(R->this.areas_of_interest.add(R));
				}
			else
				{
				splitted_intervals = false;
				this.areas_of_interest.add(new WholeGenomeArea(this.dictionary));
				}
			
			
			final List<Path> bedPaths= IOUtils.unrollPaths(args);
			
			if(bedPaths.isEmpty()) {
				try(BufferedReader r=new BufferedReader(new InputStreamReader(stdin(), "UTF-8"))) {
					readBed(r,null);
					}
				}
			else
				{
				for(Path bedPath : bedPaths) {
					readFileAnyFormat(bedPath);
					}
				}
			
			
			
			final SVGDocument svg = new SVGDocument();
			if(this.areas_of_interest.size()==1) {
				svg.setTitle(this.areas_of_interest.get(0).getName());
				}
			
			svg.styleElement.appendChild(svg.text(
				"rect.contig0 {stroke:dimgray;stroke-width:0.8px;fill:gainsboro;}\n" +
				"rect.contig1 {stroke:dimgray;stroke-width:0.8px;fill:lightgrey;}\n" +
				".tracklabel {text-anchor:end; font-size:12px;} "+
				".sep {stroke-width:0.8px;} " +
				".areaframe {fill:none;stroke:stroke:blue;} " +
				"rect.background0 {fill:white;stroke:none;} " +
				"rect.background1 {fill:gainsboro;stroke:none;}"
				));
			
			double y=0;
			final Element document_g= svg.group();
			svg.rootElement.appendChild(document_g);
			final Insets insets = new Insets(50, 50, 50, 50);
			
			final Element main_g= svg.group();
			document_g.setAttribute("id", "main");
			main_g.setAttribute("transform", "translate("+insets.top+","+insets.left+")");
			document_g.appendChild(main_g);
			
			/** draw chromosomes at TOP */
			final long genome_length = this.dictionary.getReferenceLength();
			Map<String,Long> chromosome2offset=new HashMap<>();
			double featureHeight=10.0;
			long offset = 0L;
			final Element contigs_g = svg.group();
			main_g.appendChild(contigs_g);
			
			for(SAMSequenceRecord ssr:this.dictionary.getSequences()) {
				chromosome2offset.put(ssr.getContig(), offset);
				final double x = (offset/(double)genome_length)*this.screen_width;
				final double w = (ssr.getLengthOnReference()/(double)genome_length)*this.screen_width;
				
				final Element g= svg.group();
				final Element rect = svg.rect(
					x,
					y,
					w,
					featureHeight,
					Maps.of("class", "contig"+(ssr.getSequenceIndex()%2))
					);
				g.appendChild(rect);
				
				final Element text = svg.text(
					x+w/2.0,
					y+featureHeight*2,
					ssr.getContig(),Maps.of(
					"class", "contig"+(ssr.getSequenceIndex()%2),
					"text-anchor","middle",
					"font-size",(int)featureHeight
					));
				g.appendChild(text);
				
				svg.setTitle(g, ssr.getContig());
				contigs_g.appendChild(svg.anchor(g, hyperlink.apply(ssr).orElse(null)));
				
				offset+=ssr.getLengthOnReference();
				}
			y+=featureHeight*2;
			/*
			 * DRAW cytobands, if any
			 */
			if(this.cytoBandFile!=null) {
				try(BufferedReader br=IOUtils.openPathForBufferedReading(this.cytoBandFile)) {
					final Pattern tab=Pattern.compile("[\t]");
					List<Cytoband> cytobands= br.lines()
						.map(L->tab.split(L))
						.map(T->Cytoband.of(T))
						.filter(C->dictionary.getSequence(C.getContig())!=null)
						.sorted(new ContigDictComparator(this.dictionary).createLocatableComparator())
						.collect(Collectors.toList());
					for(Cytoband cytoband: cytobands) {
						Element rect = svg.rect(
								(cytoband.getStart()/(double)genome_length)*this.screen_width,
								y,
								(cytoband.getLengthOnReference()/(double)genome_length)*this.screen_width,
								featureHeight,
								Maps.of("", "")
								);
						svg.setTitle(rect, cytoband.getName());
						main_g.appendChild(svg.anchor(rect, hyperlink.apply(cytoband).orElse(null)));
						}
					}
				}
			
			y+=featureHeight;
			
			/*
			 * Loop over each track
			 */
			for(Track track:this.id2track.values()) {
				y+=2;
				final Element track_g = svg.group();
				track.y = y;
				final List<Element> frames = new ArrayList<Element>();
				for(int area_idx=0;area_idx< this.areas_of_interest.size();++area_idx) {
					final AreaOfInterest area = this.areas_of_interest.get(area_idx);
					final Element area_g = svg.group();
					
					
					
					
					
					// loop over each contig, because pileup doesn't handle multiple contigs
					for(String contig: area.getContigs()) {
						
						
						final Pileup<Locatable> pileup = new Pileup<Locatable>(area.createNCollisionTest());
						for(Locatable loc:track.records) {
							if(!loc.getContig().equals(contig)) continue;
							if(!area.accept(loc)) continue;
							pileup.add(loc);
							}
						double area_y=track.y+2;
						for(List<Locatable> row: pileup) {
							for(Locatable loc: row) {
								final double x1 = area.genomic2pixel(loc.getContig(),area.trim(loc.getContig(), loc.getStart()));
								final double x2 = area.genomic2pixel(loc.getContig(),area.trim(loc.getContig(), loc.getEnd()+1));
								final Element rect= svg.rect(
										x1,
										area_y,
										Math.max(0.5, x2-x1),
										featureHeight*0.9,
										getAttributesFor(loc)
										);
								area_g.appendChild(svg.anchor(rect, hyperlink.apply(loc).orElse(null)));
								}
							area_y+=featureHeight;
							}
						track.max_y = Math.max(track.max_y,area_y);
						}
					
					final Element frame = svg.rect(
							area.getX(),
							track.y,
							area.getWidth(),
							1,
							Maps.of("class", "areaframe")
							);
					frames.add(frame);
					area_g.appendChild(frame);
					
					track_g.appendChild(area_g);
					}
				track.max_y+=2;
				
				//update frame height
				for(Element frame: frames) {
					frame.setAttribute("height", svg.format(track.max_y-track.y));
					}
				
				
				final Element text = svg.text(
						-5,
						(track.y+track.max_y)/2.0,
						track.name,
						Maps.of("class","tracklabel"));
				track_g.appendChild(text);
				
				y = track.max_y+featureHeight;
				main_g.appendChild(track_g);
				}
			
			
			svg.setWidth(this.screen_width+insets.left+insets.right+1);
			svg.setHeight(y+insets.top+insets.bottom+1);
			
			if(this.outputFile==null) {
				svg.saveTo(stdout());
				}
			else
				{
				svg.saveTo(this.outputFile);
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new BedToSVG().instanceMainWithExit(args);
		}
	
}
