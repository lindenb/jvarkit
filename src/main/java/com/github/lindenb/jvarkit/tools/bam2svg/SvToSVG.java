/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.File;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Result;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.DocumentFragment;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.Text;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.lang.StaticCodeExtractor;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.svg.SVG;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;


/**
BEGIN_DOC

## Example

```bash
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA20845/high_coverage_alignment/NA20845.wgs.ILLUMINA.bwa.GIH.high_cov_pcr_free.20140203.bam" "7:8352157-8356597" > jeter.bam && samtools index jeter.bam
$ java -jar dist/svg2svg.jar jeter.bam > jeter.svg
```

### A translocation


Translocation described in [PMC5932280](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5932280/)

> Identification of Balanced Chromosomal Rearrangements Previously Unknown Among Participants in the 1000 Genomes Project: Implications for Interpretation of Structural Variation in Genomes and the Future of Clinical Cytogenetics

```
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02260/alignment/HG02260.mapped.ILLUMINA.bwa.PEL.low_coverage.20130415.bam" "9:137229907-137231907" > jeter1.bam
$ samtools view -b "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02260/alignment/HG02260.mapped.ILLUMINA.bwa.PEL.low_coverage.20130415.bam" "14:79838174-79840174" > jeter2.bam
$ samtools merge jeter3.bam jeter1.bam jeter2.bam
$ samtools index jeter3.bam
$ java -jar dist/sv2svg.jar -r "9:137229907-137231907" -r "14:79838174-79840174"  jeter3.bam > jeter.svg
```


## Gallery

[https://gist.github.com/lindenb/bf48989b8da31eeafdc2caa0694361eb](https://gist.github.com/lindenb/bf48989b8da31eeafdc2caa0694361eb)

[https://twitter.com/yokofakun/status/1063369955406725120](https://twitter.com/yokofakun/status/1063369955406725120)

![https://imgur.com/EVOrXuc](https://i.imgur.com/EVOrXuc.gif)

[https://twitter.com/yokofakun/status/1063511215832539136](https://twitter.com/yokofakun/status/1063511215832539136)

![https://pbs.twimg.com/media/DsJZKRrWsAE_QqA.jpg](https://pbs.twimg.com/media/DsJZKRrWsAE_QqA.jpg)

[https://gist.github.com/lindenb/877d1d00d9f19c618f2d8505a2fe5614](https://gist.github.com/lindenb/877d1d00d9f19c618f2d8505a2fe5614)

[https://twitter.com/yokofakun/status/1064484537059684355](https://twitter.com/yokofakun/status/1064484537059684355)

![https://twitter.com/yokofakun/status/1064484537059684355](https://pbs.twimg.com/media/DsXObwuXoAAepSw.jpg)

[https://twitter.com/yokofakun/status/1064503996285681666](https://twitter.com/yokofakun/status/1064503996285681666)

![https://pbs.twimg.com/media/DsXgnKrWwAAYiBS.jpg](https://pbs.twimg.com/media/DsXgnKrWwAAYiBS.jpg)


[https://gist.github.com/lindenb/88bb702478cafe732d00e2694e77bc09](https://gist.github.com/lindenb/88bb702478cafe732d00e2694e77bc09)

![https://twitter.com/yokofakun/status/1448993550746759169](https://pbs.twimg.com/media/FBva8sTWUAMNApn?format=jpg&name=large)

END_DOC
 */
@Program(name="sv2svg",
description="BAM to SVG. Used to display structural variations.",
keywords={"bam","alignment","graphics","visualization","svg","structural","variant"},
modificationDate="20211015",
creationDate="20181115",
biostars=405059
)
public class SvToSVG extends Launcher
	{
	private static final Logger LOG = Logger.build(SvToSVG.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-r","-i","--interval","--region"},description="interval CHROM:START-END",required=true)
	private List<String> intervalStrList = new ArrayList<>();
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION+". Optional: if defined will be used to display the mismatches.")
	private Path fastaFile = null ;
	@Parameter(names={"-d","--duration"},description="Animation duration, in secs. <=0 disable animation.")
	private int svgDuration=10;
	@Parameter(names= {"--repeat-count"},description="SVG animation repeat count")
	private String svgRepeatCount="indefinite";
	@Parameter(names= {"--variant","-V"},description="optional indexed VCF file.")
	private Path vcfFile=null;
	@Parameter(names= {"--coverage","--depth"},description="Coverage height. Don't print if cov '<=0'.")
	private int coverageHeight=70;
	@Parameter(names= {"--mapq","--Q"},description="min mapping quality")
	private int min_mapq = 1;
	@Parameter(names= {"--mismatch"},description="do not display bases mismatches between read and REF.")
	private boolean hideMismatch = false;
	@Parameter(names= {"--strange"},description="Keep only non-properly-paired , soft clipped and SA:X:* reads.")
	private boolean remove_simple_reads = false;

	
	private final List<Sample> sampleList =new ArrayList<>();
	private double featureHeight = 10;
	private final DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private final DecimalFormat niceIntFormat = new DecimalFormat("###,###");
	private Document document = null;
	private final double arrow_w = 5;
	private ReferenceSequenceFile indexedFastaSequenceFile = null;
	private VCFReader vcfFileReader = null;
	
	private final String DEBUG_READ="___";
	
	private abstract class YHeight {
		public abstract double getY();
		public abstract double getHeight();
	}
	
	private class Sample extends YHeight {
		final String sampleName;
		final List<Region> regions = new ArrayList<>();
		double y = 0;
		double height = 0.0;
		
		/** abstract base class for a short read, wrap a SAMRecord and a SVG element */
		private abstract class ShortRead
			extends YHeight
			implements Locatable
			{
			Element g_translate1 = null;
			double y;
			final Region region;
			
			protected ShortRead(final Region region) {
				this.region = region;
				}
			
			abstract SAMRecord getRecord();
			
			@Override
			public double getY() {
				return this.y;
				}
			@Override
			public double getHeight() {
				return featureHeight;
				}

			
			@Override
			public final String getContig() {
				return getRecord().getContig();
				}
			@Override
			public final int getStart() {
				return getRecord().getUnclippedStart();
				}
			@Override
			public final int getEnd() {
				return getRecord().getUnclippedEnd();
				}
			String getSampleName() {
				final SAMReadGroupRecord rg = getRecord().getReadGroup();
				if(rg==null) return null;
				final String sn = rg.getSample();
				return sn;
				}
			final  boolean isNegativeStrand() {
				return getRecord().getReadNegativeStrandFlag();
				}
			public final double getWidth() {
				return getPixelEnd()-getPixelStart();
				}
			public final double getCenterX() {
				return (getPixelStart()+getPixelEnd())/2.0;
			}
			public final double getCenterY() {
				return y + getHeight()/2.0;
			}

			
			public final double getPixelStart() {
				return getRegion().baseToPixel(this.getStart());
			}
			
			public final double getPixelEnd() {
				return getRegion().baseToPixel(this.getEnd());
			}
			
			public final double getPixelLength() {
				return getPixelEnd() - getPixelStart();
			}
			
			
			abstract boolean isDefaultShortRead();
			abstract boolean isSplitRead();
			/* get owner */
			final Region getRegion() {
				return this.region;
			}
			@Override
			public String toString() {
				return getRecord().getReadName()+"["+getContig()+":"+getStart()+"-"+getEnd()+"] in-region:"+getRegion()+" split:"+isSplitRead()+ " y:"+y;
				}
			}
		
		/* default implementation of ShortRead */
		private class DefaultShortRead
			extends ShortRead
			{
			final SAMRecord record;
			
			DefaultShortRead(final Region region,final SAMRecord record) {
				super(region);
				this.record = record;
				}


			@Override
			SAMRecord getRecord() {
				return this.record;
				}
			
			@Override boolean isDefaultShortRead() {
				return true;
				}
			@Override boolean isSplitRead() {
				return false;
				}
			}
		
		/* implementation of SplitRead , wraps supplementary and it's source */
		private class SplitRead extends ShortRead
			{
			/* bad name for those 'g' element, in fact g_translate2 contains the rotate animation */
			Element g_rotate = null;
			Element g_translate2 = null;

			final SAMRecord record;
			final DefaultShortRead source;
			SplitRead(final Region region,final SAMRecord record,final DefaultShortRead source) {
				super(region);
				this.record = record;
				this.source = source;
				}
			
			@Override
			SAMRecord getRecord() {
				return this.record;
				}
			
			@Override boolean isDefaultShortRead() {
				return false;
				}
			@Override boolean isSplitRead() {
				return true;
				}
			@Override
			public String toString() {
				return super.toString()+ (this.source.getRegion()==this.getRegion()?" => ":" -> ")+this.source.toString();
				}
			}

		
		/** a region */
		private class Region extends YHeight implements Locatable {
			double y;
			double height;
			final SimpleInterval _interval;
			final List<ShortRead> beforePileup=new ArrayList<>();
			final List<List<ShortRead>> lines=new ArrayList<>();			
			
			Region(final SimpleInterval interval) {
				this._interval = interval;
				}
			
			@Override
			public double getY() {
				return y;
				}
			
			@Override
			public double getHeight() {
				return height;
				}
			
			@Override
			public String getContig() {
				return _interval.getContig();
				}
			@Override
			public int getStart() {
				return _interval.getStart();
				}
			@Override
			public int getEnd() {
				return _interval.getEnd();
				}
			
			
			void addSplitReads() {
				final List<SplitRead> addToSelf = new ArrayList<>();
				for(final ShortRead shortRead : this.beforePileup) {
					if(!shortRead.isDefaultShortRead()) continue;
					final List<SAMRecord> others = SAMUtils.getOtherCanonicalAlignments(shortRead.getRecord());
					for(final SAMRecord rec:others) {
						final Cigar cigar = rec.getCigar();
						if(cigar==null || cigar.isEmpty()) continue;
						/* loop over all regions for this samples */
						for(final Region otherRgn : Sample.this.regions) {
							/* check supplementary reads belong to that region */
							if(!rec.getContig().equals(otherRgn.getContig())) continue;
							/* check supplementary reads belong to that interval */
							if(!CoordMath.overlaps(
									rec.getUnclippedStart(), rec.getUnclippedEnd(),
									otherRgn.getStart(), otherRgn.getEnd())) continue;
							if(otherRgn==this) {
								//cannot add to self for now because we're looping
								addToSelf.add(new SplitRead(this,rec,(DefaultShortRead)shortRead));
								}
							else {
								final SplitRead splitRead = new SplitRead(otherRgn,rec,(DefaultShortRead)shortRead);
								otherRgn.beforePileup.add(splitRead);
								}
							}
						}
					}
				this.beforePileup.addAll(addToSelf);
				}
			
			/* pileup reads in that region */
			void pileup() {
				 this.beforePileup.sort((A,B)->Integer.compare(A.getStart(), B.getStart()));
				
				for(final ShortRead shortRead : this.beforePileup) {
					/* pileup */
					int y=0;
					for(y=0;y< this.lines.size();++y)
						{
						final List<ShortRead> line= this.lines.get(y);
						final ShortRead last = line.get(line.size()-1);
						if( baseToPixel(last.getEnd()) + 2*arrow_w < baseToPixel(shortRead.getStart()) )
							{
							line.add(shortRead);
							break;
							}
						}
					
					if( y == this.lines.size())
						{
						final List<Sample.ShortRead> line= new ArrayList<>();
						line.add(shortRead);
						this.lines.add(line);
						}
					}
				// not needed anymore
				this.beforePileup.clear();
				}
			
			double baseToPixel(final int pos)
				{
				return  ((pos - this.getStart())/(double)(this.getLengthOnReference()))*(SvToSVG.this.drawinAreaWidth);
				}
			Stream<ShortRead> shortReadStream() {
				return this.lines.stream().flatMap(L->L.stream());
				}
			
			@Override
			public String toString() {
				return this._interval.toString();
				}
			} /* end region */
		
		Sample(final String sampleName) {
			this.sampleName = sampleName;
			}
		
		Stream<ShortRead> shortReadStream() {
			return this.regions.stream().flatMap(R->R.shortReadStream());
			}
		
		
		@Override
		public double getY() {
			return this.y;
			}
		@Override
		public double getHeight() {
			return this.height;
			}
		
		private List<Sample.ShortRead> getReadsByName(final String readName) {
			return this.shortReadStream().
					filter(S->S.getRecord().getReadName().equals(readName)).
					collect(Collectors.toList());
			}
		}/* end sample */
	
	
	/** convert double to string */
	private String format(final double v)
		{
		return this.decimalFormater.format(v);
		}
	
	
	private Element buildSample(
			final Sample sample,
			final DocumentFragment animationLayer
			)
		{
		final Element sampleRoot = element("g");
		double curr_y = sample.getY();
		final Element sampleLabel= element("text",sample.sampleName);
		sampleLabel.setAttribute("x", "5");
		sampleLabel.setAttribute("y", format(curr_y +14));
		sampleLabel.setAttribute("class", "samplename");
		sampleRoot.appendChild(sampleLabel);
		curr_y += 20;
		
	
		
		// loop over regions 
		for(final Sample.Region region : sample.regions) {
			region.y = curr_y;
			// genomic sequence will be used to test if there is a mismatch. No garantee it exists.
			final GenomicSequence genomicSequence;
			if(this.indexedFastaSequenceFile !=null) {
				final ContigNameConverter ctgConver = ContigNameConverter.fromOneDictionary(this.indexedFastaSequenceFile.getSequenceDictionary());
				final String sname = ctgConver.apply(region.getContig());
				if(StringUtil.isBlank(sname)) {
					genomicSequence = null;
				} else
					{
					genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, sname);
					}
				}
			else
				{
				genomicSequence = null;
				}
			
			
			final Element regionRoot = element("g");
			sampleRoot.appendChild(regionRoot);
			
			final Element rgnLabel= element("text",
					new SimpleInterval(region).toNiceString()+" Length:"+this.niceIntFormat.format((region.getLengthOnReference())));
			rgnLabel.setAttribute("x", "5");
			rgnLabel.setAttribute("y", format(curr_y+12));
			rgnLabel.setAttribute("class", "samplename");
			regionRoot.appendChild(rgnLabel);
			curr_y += 20;
			
			final double y_top_region = curr_y;
			
			if(this.coverageHeight>0)
				{
				/* draw coverage */
				final TreeMap<Integer,Long> pos2cov = region.
					shortReadStream().
					filter(R->R.isDefaultShortRead()).
					flatMapToInt(R->R.getRecord().getAlignmentBlocks().stream().flatMapToInt(RB->java.util.stream.IntStream.rangeClosed(RB.getReferenceStart(),RB.getReferenceStart()+RB.getLength()))).
					filter(P->P>=region.getStart()).
					filter(P->P<=region.getEnd()).
					mapToObj(P->P).
					collect(Collectors.groupingBy( Function.identity(),()->new TreeMap<>(), Collectors.counting()));
					;
				final long max_cov = pos2cov.values().stream().mapToLong(L->L.longValue()).max().orElse(1L);
				final Element covPath = element("path");
				covPath.setAttribute("class", "coverage");
				regionRoot.appendChild(covPath);
				final StringBuilder sb = new StringBuilder();
				sb.append( "M 0 "+format(curr_y+this.coverageHeight));
				for(int k=0;k< this.drawinAreaWidth;k++)
					{
					final int pos1 = region.getStart()+ (int)(((k+0)/(double)this.drawinAreaWidth)*(double)region.getLengthOnReference());
					final int pos2 = region.getStart()+ (int)(((k+1)/(double)this.drawinAreaWidth)*(double)region.getLengthOnReference());
					final double dp = IntStream.range(pos1, pos2).
						filter(P->pos2cov.containsKey(P)).
						mapToLong(P->pos2cov.get(P)).
						max().
						orElseGet(()->0L);
					final double dpy= curr_y + this.coverageHeight - (dp/(double)max_cov)*(double)this.coverageHeight;
					sb.append(" L "+format(k)+" "+format(dpy));
					}
				sb.append(" L "+format(this.drawinAreaWidth)+" "+format(curr_y+this.coverageHeight));
				sb.append(" Z");
				covPath.setAttribute("d", sb.toString());
				covPath.appendChild(element("title","Coverage. Max:"+niceIntFormat.format(max_cov)));
				curr_y += this.coverageHeight;
				curr_y += 2;
				}
				
			
			
			/* print all lines */
			for(int nLine=0;nLine< region.lines.size();++nLine)
				{
				final List<Sample.ShortRead> line= region.lines.get(nLine);
				//loop over records on that line
				for(final Sample.ShortRead shortRead: line)
					{
					boolean trace = shortRead.getRecord().getReadName().equals(DEBUG_READ);
					shortRead.y  = curr_y;

					
					int ref1= shortRead.getStart();
					final double leftX = shortRead.getPixelStart();
					final double midy =  this.featureHeight/2.0;
					final double maxy = this.featureHeight;

					final Element g1 = element("g");
					final Element g3;
					if(shortRead.isSplitRead()) {
							final Sample.SplitRead splitRead =  Sample.SplitRead.class.cast(shortRead);
							//final double dx = -0.5*shortRead.getPixelLength();
							//final double dy = -0.5*shortRead.getHeight();
							
							splitRead.g_rotate = element("g");
							splitRead.g_rotate.setAttribute("transform","rotate(0)");
							g1.appendChild(splitRead.g_rotate);
							g3 = element("g");
							g3.setAttribute("transform","translate(0,0)");
							splitRead.g_rotate.appendChild(g3);
							splitRead.g_translate2 = g3;
					} else {
						g3 = g1;
						}
					g1.setAttribute("transform","translate(" +
						format(shortRead.getPixelStart()) + "," +
						format(shortRead.getY()) +
						")");
					
					// https://developer.mozilla.org/en-US/docs/Web/CSS/transform-origin
					//readElement.setAttribute("style","transform-origin:center");
					if(shortRead.isSplitRead()) { //always in front
						animationLayer.appendChild(g1);
						}
					else
						{
						regionRoot.insertBefore(g1, regionRoot.getFirstChild());
						}
					
					final BiPredicate<Integer, Integer> isMismatch = (readpos0,refpos1)-> {
						if(genomicSequence==null) return false;
						if(refpos1<1 || refpos1>genomicSequence.length()) return false;
						char refC = Character.toUpperCase(genomicSequence.charAt(refpos1-1));
						if(refC=='N') return false;
						final byte bases[] = shortRead.getRecord().getReadBases();
						if(bases==null || bases==SAMRecord.NULL_SEQUENCE ) return false;
						if(readpos0<0 || readpos0>=bases.length) {
							LOG.warn("P="+readpos0+" 0<x<"+bases.length);
							return false;
						}
						char readC = (char)Character.toUpperCase(bases[readpos0]);
						if(readC=='N') return false;
						return readC!=refC;
						};
					
					shortRead.g_translate1 = g1;
					int readpos0=0;
					int readpos0_clipped=0;
					
					
					final Element title = element("title",shortRead.getRecord().getPairedReadName()+" "+shortRead.getRecord().getCigarString()+" "+(shortRead.getRecord().getReadNegativeStrandFlag()?"-":"+"));
					g3.appendChild(title);
					
					final DocumentFragment insertionsFragment = this.document.createDocumentFragment();
					
					final Cigar cigar = shortRead.getRecord().getCigar();
					for(int cigarIdx=0;cigarIdx< cigar.numCigarElements();cigarIdx++) {
						final CigarElement ce = cigar.getCigarElement(cigarIdx);
						final CigarOperator op = ce.getOperator();
						int next_ref = ref1;
						int next_read = readpos0;
						switch(op)
							{
							case I: 
								{
								next_read += ce.getLength();
								readpos0_clipped+= ce.getLength();
								if(!(ref1+ce.getLength() < region.getStart() || ref1> region.getEnd())) {
									final double ce_length = region.baseToPixel(region.getStart()+ce.getLength());
									final Element rectInsert = element("rect");
									rectInsert.setAttribute("class", "insert");
									rectInsert.setAttribute("x", format(region.baseToPixel(ref1)-leftX));
									rectInsert.setAttribute("y", format(0));
									rectInsert.setAttribute("width", format(this.svgDuration>0  && ce_length >1 ?ce_length:1));
									rectInsert.setAttribute("height", format(this.featureHeight));
									rectInsert.appendChild(element("title", this.niceIntFormat.format(ce.getLength())));
									if(this.svgDuration>0 && ce_length >1)
										{
										final Element anim = element("animate");
										rectInsert.appendChild(anim);
										anim.setAttribute("attributeType","XML");
										anim.setAttribute("attributeName","width");
										anim.setAttribute("begin","0s");
										anim.setAttribute("from",format(ce_length));
										anim.setAttribute("to","1");
										anim.setAttribute("dur",String.valueOf(this.svgDuration)+"s");
										anim.setAttribute("repeatCount",this.svgRepeatCount);
										anim.setAttribute("fill","freeze");
										}
									
									insertionsFragment.appendChild(rectInsert);
									}
								break;
								}
							case P: continue;
							case S:case H:
							case M:case X: case EQ:
								{
								next_ref+= ce.getLength();
								
								if(!op.equals(CigarOperator.H)) {
									next_read+= ce.getLength();
									}
								if(!(next_ref < region.getStart() || ref1> region.getEnd())) {
									final double distance_pix = region.baseToPixel(next_ref)-region.baseToPixel(ref1);
									
									
									
									final StringBuilder sb=new StringBuilder();
									
									final Element path = element("path");
									path.setAttribute("class", "op"+op.name()+(shortRead.isDefaultShortRead()?"":"x"));
									if(trace) path.setAttribute("style", "fill:blue;");
									if(!op.isClipping() && shortRead.isDefaultShortRead() && shortRead.getRecord().getReadPairedFlag() &&  !shortRead.getRecord().getProperPairFlag())
										{
										if(!shortRead.getContig().equals(shortRead.getRecord().getMateReferenceName()))
											{
											path.setAttribute("style", "fill:orchid;");
											}
										else
											{
											path.setAttribute("style", "fill:lightblue;");
											}
										}
									
									// arrow <--
									if(cigarIdx==0 && shortRead.isNegativeStrand()) {
										sb.append( "M ").append(format(region.baseToPixel(ref1)-leftX)).append(',').append(0);
										sb.append(" h ").append(format(distance_pix));
										sb.append(" v ").append(format(maxy));
										sb.append(" h ").append(format(-(distance_pix)));
										sb.append(" l ").append(format(-arrow_w)).append(',').append(-featureHeight/2.0);
										sb.append(" Z");
										}
									// arrow -->
									else if(cigarIdx+1==cigar.numCigarElements() && !shortRead.isNegativeStrand()) {
										sb.append( "M ").append(format(region.baseToPixel(ref1)-leftX+distance_pix)).append(',').append(0);
										sb.append(" h ").append(format(-(distance_pix)));
										sb.append(" v ").append(format(maxy));
										sb.append(" h ").append(format(distance_pix));
										sb.append(" l ").append(format(arrow_w)).append(',').append(-featureHeight/2.0);
										sb.append(" Z");
										
										}
									else
										{
										sb.append( "M ").append(format(region.baseToPixel(ref1)-leftX)).append(',').append(0);
										sb.append(" h ").append(format(distance_pix));
										sb.append(" v ").append(format(maxy));
										sb.append(" h ").append(format(-(distance_pix)));
										sb.append(" Z");
										}
									path.setAttribute("d", sb.toString());
									g3.appendChild(path);
									
									if(op.isAlignment() && genomicSequence!=null && !hideMismatch)
										{
										for(int x=0;x< ce.getLength();++x)
											{
											if(!isMismatch.test(readpos0_clipped+x, ref1+x)) continue;
											if(ref1+x < region.getStart()) continue;
											if(ref1+x > region.getEnd()) continue;
											final Element rectMismatch = element("rect");
											rectMismatch.setAttribute("class", "mismatch");
											rectMismatch.setAttribute("x", format(region.baseToPixel(ref1+x)-leftX));
											rectMismatch.setAttribute("y", format(0));
											rectMismatch.setAttribute("width", format((region.baseToPixel(ref1+x+1)-region.baseToPixel(ref1+x))));
											rectMismatch.setAttribute("height", format(this.featureHeight));
											g3.appendChild(rectMismatch);
											}
										}
									}
								if(!op.isClipping()) {
									readpos0_clipped+=ce.getLength();
									}
								break;
								}
							case D: case N:
								{
								next_ref+= ce.getLength();
								if(!(next_ref < region.getStart() || ref1> region.getEnd())) {
									final double distance_pix = region.baseToPixel(next_ref)-region.baseToPixel(ref1);
									final Element lineE = element("line");
									lineE.setAttribute("class", "opD");
									lineE.setAttribute("x1", format(region.baseToPixel(ref1)-leftX));
									lineE.setAttribute("y1", format(midy));
									lineE.setAttribute("x2", format(region.baseToPixel(ref1)-leftX+distance_pix));
									lineE.setAttribute("y2", format(midy));
									g3.insertBefore(lineE, g3.getFirstChild());
									}
								break;
								}
							default: throw new IllegalStateException(op.name());
							}
						ref1 = next_ref;
						readpos0 = next_read; 
						}
					g3.appendChild(insertionsFragment);
					}
				
				curr_y += (this.featureHeight+3);
				}
			
			//add non-properly paired reads
			/* region.shortReadStream().
				filter(R->R.isDefaultShortRead() && 
						R.getRecord().getReadPairedFlag() && 
						R.getRecord().getFirstOfPairFlag() && 
						!R.getRecord().getProperPairFlag()).
				forEach(R->{
					final Sample.Region.ShortRead mate = region.shortReadStream().filter(
							R2->R2.isDefaultShortRead() &&
							R2.getRecord().getReadPairedFlag() &&
							R2.getRecord().getSecondOfPairFlag() &&
							!R2.getRecord().getProperPairFlag() && 
							R2.getRecord().getReadName().equals(R.getRecord().getReadName())).
							findFirst().orElse(null);
					if(mate==null) return;
					final double x1 = R.isNegativeStrand()?R.getPixelStart()-arrow_w:R.getPixelEnd()+arrow_w;
					final double x2 = mate.isNegativeStrand()?mate.getPixelStart()-arrow_w:mate.getPixelEnd()+arrow_w;
					
					final Element line = element("line");
					line.setAttribute("class", "discordant");
					line.setAttribute("x1", format(x1));
					line.setAttribute("y1", format(R.y+this.featureHeight/2.0));
					line.setAttribute("x2", format(x2));
					line.setAttribute("y2", format(mate.y+this.featureHeight/2.0));	
					animationLayer.appendChild(line);
					});*/
			
			
				for(int x=1;x<10;++x)
					{
					final double xx = (this.drawinAreaWidth/10.0)*x;
					final Element line = element("line");
					line.setAttribute("class", "ruler");
					line.setAttribute("x1", format(xx));
					line.setAttribute("y1", format(y_top_region));
					line.setAttribute("x2", format(xx));
					line.setAttribute("y2", format(curr_y));
					line.appendChild(element("title",niceIntFormat.format(region.getStart()+(int)(region.getLengthOnReference()/10.0)*x)));
					regionRoot.insertBefore(line, regionRoot.getFirstChild());
					}
				
				/** print variants */
				if(this.vcfFileReader!=null)
					{
					final VCFHeader header = this.vcfFileReader.getHeader();
					final SAMSequenceDictionary vcfdict = header.getSequenceDictionary();
					final ContigNameConverter ctgConver = (vcfdict==null?null:ContigNameConverter.fromOneDictionary(vcfdict));
					final String sname = ctgConver==null?null:ctgConver.apply(region.getContig());
					if(!StringUtil.isBlank(sname))
						{
						final CloseableIterator<VariantContext> iter = this.vcfFileReader.query(region.getContig(), region.getStart(), region.getEnd());
						while(iter.hasNext())
							{
							final VariantContext ctx = iter.next();
							if(!ctx.isVariant()) continue;
							final double x1 = Math.max(0, region.baseToPixel(ctx.getStart()));
							final double x2 = Math.min(region.baseToPixel(ctx.getEnd()+1),this.drawinAreaWidth);
							final Genotype g = ctx.getGenotype(sample.sampleName);
							
							final Element rect = element("rect");
							rect.setAttribute("class", "variant"+(g==null?"":g.getType().name()));
							rect.setAttribute("x", format(x1));
							rect.setAttribute("y", format(y_top_region));
							rect.setAttribute("width", format(x2-x1));
							rect.setAttribute("height", format(curr_y-y_top_region));
							rect.appendChild(element("title",
									niceIntFormat.format(ctx.getStart())+"-"+
											niceIntFormat.format(ctx.getEnd())
											+" "+ctx.getReference().getDisplayString()));
							regionRoot.insertBefore(rect, regionRoot.getFirstChild());
							
							}
						iter.close();
						}
					}
				region.height =  curr_y - region.y;
				}
		
		
		sample.height = curr_y - sample.y;

		
		final Element frame= element("rect");
		frame.setAttribute("class","frame");
		frame.setAttribute("x",format(0));
		frame.setAttribute("y",format(sample.getY()));
		frame.setAttribute("width",format(drawinAreaWidth));
		frame.setAttribute("height",format(sample.getHeight()));
		sampleRoot.appendChild(frame);
		
		return sampleRoot;
		}
		
	
		private Element element(final String tag) {
			return this.document.createElementNS(SVG.NS, tag);
			}
		private Text text(final Object o) {
			return this.document.createTextNode(o==null?"":String.valueOf(o));
			}
		private Element element(final String tag,final Object content) {
			final Element E = element(tag);
			E.appendChild(text(content));
			return E;
			}
		private void buildDocument() 
			{
			final Element svgRoot = element("svg");
			this.document.appendChild(svgRoot);
			svgRoot.setAttribute("width",format(this.drawinAreaWidth+1));
			
			int doc_height =0;
			
			final Element title = element("title");
			svgRoot.appendChild(title);
			title.appendChild(text(String.join(" ",this.intervalStrList)));
			
			final Element descr = element("desc");
			svgRoot.appendChild(descr);
			descr.appendChild(text("Author: Pierre Lindenbaum"));
			
			final Element style = element("style");
			svgRoot.appendChild(style);
			style.appendChild(text(StaticCodeExtractor.forClass(SvToSVG.class).extract("CSS").get()));
			
			final Element mainG = element("g");
			mainG.setAttribute("class","maing");
			svgRoot.appendChild(mainG);
			
			final DocumentFragment animationLayer = this.document.createDocumentFragment();
			
			//loop over each sample
			for(final Sample sample:this.sampleList)
				{
				sample.y = doc_height;
				final Element div = buildSample(sample,animationLayer);
				doc_height += sample.getHeight();
				mainG.appendChild(div);
				}
			
			//loop over discordant reads and print arrows
			for(final Sample sample:this.sampleList)
				{
				sample.shortReadStream().
					filter(R->R.isDefaultShortRead()).
					filter(R->R.getRecord().getReadPairedFlag()).
					filter(R->!R.getRecord().getMateUnmappedFlag()).
					filter(R->!R.getRecord().getProperPairFlag()).
					//filter(R->!R.getRecord().getContig().equals(R.getRecord().getMateReferenceName())).
					map(SR->SR.getRecord().getReadName()).
					collect(Collectors.toSet()).
					forEach(READNAME->{
					final List<Sample.ShortRead> sr = sample.getReadsByName(READNAME);
					for(int x=0;x+1 < sr.size();++x)
						{
						final Sample.ShortRead srx = sr.get(x);
						if(!srx.isDefaultShortRead()) continue;
						if(!srx.getRecord().getFirstOfPairFlag()) continue;
						for(int y=x+1;y  < sr.size();++y)
							{
							final  Sample.ShortRead sry = sr.get(y);
							if(!sry.isDefaultShortRead()) continue;
							if(!sry.getRecord().getSecondOfPairFlag()) continue;
							//if(srx.getContig().equals(sry.getContig())) continue;
							
							final double x1 = srx.isNegativeStrand()?srx.getPixelStart()-arrow_w:srx.getPixelEnd()+arrow_w;
							final double x2 = sry.isNegativeStrand()?sry.getPixelStart()-arrow_w:sry.getPixelEnd()+arrow_w;
							
							final Element path = element("path");
							path.setAttribute("class", "discordant");
							final StringBuilder sb = new StringBuilder();
							sb.append( "M ").append(format(x1)).append(" ").append(format(srx.y+featureHeight/2.0));
							sb.append( "Q ").
								append(format((x1+x2)/2.0)).append(" ").append(format((srx.y+sry.y)/2.0)).
								append(" ").
								append(format(x2)).append(" ").append(format(sry.y+featureHeight/2.0));
							path.setAttribute("d", sb.toString());
							animationLayer.appendChild(path);
							}
						}
					});
				}
			
			// move split reads
			this.sampleList.stream().
				flatMap(SN->SN.shortReadStream()).
				filter(R->R.isSplitRead()).
				map(R->Sample.SplitRead.class.cast(R)).
				forEach(SR->{
					if(this.svgDuration<=0) return;
					
					boolean rotate = SR.getContig().equals(SR.source.getContig()) && SR.isNegativeStrand()!=SR.source.isNegativeStrand();
					//add a few 0.01 seconds to avoir same shape when rotation
					final String beginSec = format(Math.random()/2.0)+"s";
					
					Element anim = element("animateTransform");
					SR.g_rotate.appendChild(anim);
					anim.setAttribute("attributeType","XML");
					anim.setAttribute("attributeName","transform");
					anim.setAttribute("type","translate");
					anim.setAttribute("begin",beginSec);
					anim.setAttribute("from","0,0");
					anim.setAttribute("to",
							format(SR.source.getPixelStart() - SR.getPixelStart())+" "+
							format(SR.source.getY() - SR.getY())
							);
					anim.setAttribute("dur",String.valueOf(this.svgDuration)+"s");
					anim.setAttribute("repeatCount",this.svgRepeatCount);
					anim.setAttribute("fill","freeze");
					
					if(rotate)
						{
						
						anim = element("animateTransform");
						SR.g_translate2.appendChild(anim);
						
						anim.setAttribute("attributeType","XML");
						anim.setAttribute("attributeName","transform");
						anim.setAttribute("type","rotate");
						anim.setAttribute("begin",beginSec);
						anim.setAttribute("from","0 "+format(SR.getWidth()/2.0)+" "+format(SR.getHeight()/2.0));
						anim.setAttribute("to",(Math.random()<0.5?"-":"")+"180 "+format(SR.getWidth()/2.0)+" "+format(SR.getHeight()/2.0));
						anim.setAttribute("dur",String.valueOf(this.svgDuration)+"s");
						anim.setAttribute("repeatCount",this.svgRepeatCount);
						anim.setAttribute("fill","freeze");
						}

					
				});
			
			mainG.appendChild(animationLayer);
			
			// frame
			final Element frame = element("rect");
			frame.setAttribute("class", "frame");
			frame.setAttribute("x", "0");
			frame.setAttribute("y", "0");
			frame.setAttribute("width",format(this.drawinAreaWidth));
			frame.setAttribute("height",format(doc_height));
			mainG.appendChild(frame);
			
			svgRoot.setAttribute("height",format(doc_height+1));
			}
		
		
		
		
		@Override
		public int doWork(final List<String> args) {
			/* parse interval */
			if(this.intervalStrList.isEmpty() )
				{
				LOG.error("interval(s) undefined");
				return -1;
				}
			
			this.drawinAreaWidth = Math.max(100,this.drawinAreaWidth );
			try
				{
				if(this.fastaFile!=null) {
					this.indexedFastaSequenceFile =ReferenceSequenceFileFactory.getReferenceSequenceFile(this.fastaFile);
					}
				
				if(this.vcfFile!=null)
					{
					this.vcfFileReader = VCFReaderFactory.makeDefault().open(this.vcfFile,true);
					}
				
				final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				dbf.setNamespaceAware(true);
				final DocumentBuilder db = dbf.newDocumentBuilder();
				this.document=db.newDocument();

				
				final List<Path> bamFiles = IOUtils.unrollPaths(args);
				if(bamFiles.isEmpty()) {
					LOG.error("bam(s) undefined");
					return -1;
					}
				
				final SamReaderFactory srf = super.createSamReaderFactory().referenceSequence(this.fastaFile);
				
				/* loop over each bam file */
				for(final Path bamFile : bamFiles) {
					try(final SamReader sr = srf.open(bamFile)) {
						final SAMFileHeader header = sr.getFileHeader();
						final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
						final Function<String,Optional<SimpleInterval>> parser=IntervalParserFactory.newInstance().dictionary(dict).make();
						final Set<String> samples = header.getReadGroups().stream().
								map(G->G.getSample()).
								filter(S->!StringUtil.isBlank(S)).
								collect(Collectors.toSet());
						
						if(samples.size()!=1) {
							LOG.error("expected on sample in "+bamFile+" but got "+samples.size()+" "+samples);
							return -1;
							}
						final Sample sample = new Sample(samples.iterator().next());
						this.sampleList.add(sample);
						/* loop over each region for that sample */
						for(final String intervalStr:this.intervalStrList) {
							LOG.info("scanning "+intervalStr+" for "+bamFile);
							final SimpleInterval interval = parser.apply(intervalStr).orElseThrow(IntervalParserFactory.exception(intervalStr));
							if(interval==null) {
								LOG.error("Cannot parse "+intervalStr+" for "+bamFile);
								return  -1;
								}
							/* create new region */
							final Sample.Region region = sample.new Region(interval);
							if(!sample.regions.isEmpty() && region.getLengthOnReference()!=sample.regions.get(0).getLengthOnReference()) {
								LOG.warning("intervals should all have the same length");
								}
							
							sample.regions.add(region);
							
							try(CloseableIterator<SAMRecord> iter = sr.query(interval.getContig(), interval.getStart(), interval.getEnd(), false) ) {
								while(iter.hasNext())
									{
									final SAMRecord record = iter.next();
									
									boolean trace = record.getReadName().equals(DEBUG_READ);
									
									if(!SAMRecordDefaultFilter.accept(record, this.min_mapq)) continue;
									
									if(trace) LOG.debug("got1");
									final Cigar cigar = record.getCigar();
									if(cigar==null || cigar.isEmpty()) continue;
									if(trace) LOG.debug("got2");
									
									if(remove_simple_reads) {
										if((!cigar.isClipped() &&
											SAMUtils.getOtherCanonicalAlignments(record).isEmpty() &&
											(record.getReadPairedFlag() && record.getProperPairFlag())
											)) continue;
										}

									
									if( !record.getContig().equals(interval.getContig())) continue;
									if(trace) LOG.debug("got3");
									final Sample.ShortRead shortRead = sample.new  DefaultShortRead(region,record);
									
									if( shortRead.getEnd()  < interval.getStart()) continue;
									if( shortRead.getStart()  > interval.getEnd())continue;
									if(trace) LOG.debug("got4" +record.getSAMString());
									if(!sample.sampleName.equals(shortRead.getSampleName())) continue;
									region.beforePileup.add(shortRead);
									}
								}/* en iterator */
							} /* end loop interval */
						}/* end open SAM */
					}
				this.sampleList.forEach(S->S.regions.forEach(R->R.addSplitReads()));
				this.sampleList.forEach(S->S.regions.forEach(R->R.pileup()));
				
				
				/* sort per sample */
				this.sampleList.sort((A,B)->A.sampleName.compareTo(B.sampleName));
				/* sort per region */
				this.sampleList.stream().forEach(SL->SL.regions.sort((A,B)->{
					final SmartComparator sc = new SmartComparator();
					int i = sc.compare(A.getContig(), B.getContig());
					if(i!=0) return i;
					return A.getStart() -  B.getStart();
					}));
				
				
				buildDocument();
				
				final Transformer tr = TransformerFactory.newInstance().newTransformer();
				final Result result;
				
				if(this.outputFile!=null)
					{
					result = new StreamResult(this.outputFile);
					}
				else
					{
					result = new StreamResult(stdout());
					}
				tr.transform(new DOMSource(this.document),result);
				
				return 0;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(this.indexedFastaSequenceFile);
				CloserUtil.close(this.vcfFileReader);
				this.document = null;
				}
			}
		
	
	public static void main(final String[] args)
		{
		new SvToSVG().instanceMainWithExit(args);
		}
	}

/**
BEGIN_CSS

g.maing {stroke:black;stroke-width:0.5px;fill:none;}
.maintitle {stroke:blue;fill:none;}
rect.frame {stroke:darkgray;fill:none;} 
path.opEQ {stroke:black;fill:gainsboro;} 
path.opX {stroke:black;fill:tomato;}
path.opM {stroke:black;fill:gainsboro;}
path.opS {stroke:black;fill:yellow;} 
path.opH {stroke:black;fill:yellow;} 
line.opN {stroke:black;}
line.opD {stroke:black;}
path.opEQx {stroke:yellow;fill:salmon;opacity:0.8;} 
path.opXx {stroke:yellow;fill:tomato;opacity:0.8;}
path.opMx {stroke:yellow;fill:salmon;opacity:0.8;}
path.opSx {stroke:yellow;fill:yellow;opacity:0.8;} 
path.opHx {stroke:yellow;fill:yellow;opacity:0.8;} 
line.opNx {stroke:yellow;opacity:0.8;}
line.opDx {stroke:yellow;opacity:0.8;}
rect.mismatch {fill:red;opacity:0.7;}
rect.insert {fill:none;stroke:red;stroke-width:1.5px}
text.samplename {stroke:none;fill:black;stroke-width:1px;}
.discordant {stroke:darkred; fill:none;stroke-dasharray:2;}
line.ruler {stroke:darkgray;stroke-dasharray:4;fill:none;stroke-width:1;}
rect.variant  {stroke:none;fill:orange;stroke-width:1px;opacity:0.8;}
rect.variantHOM_REF  {stroke:none;fill:green;stroke-width:1px;opacity:0.8;}
rect.variantHOM_VAR  {stroke:none;fill:red;stroke-width:1px;opacity:0.8;}
rect.variantHET  {stroke:none;fill:orange;stroke-width:1px;opacity:0.8;}
rect.variantNO_CALL  {stroke:none;fill:blue;stroke-width:1px;opacity:0.8;}
path.coverage {stroke:darkslateblue;fill:darkseaGreen}


END_CSS
*/
