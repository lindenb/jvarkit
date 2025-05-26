/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.cmpbams;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineCodec;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.hershey.Hershey;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC

## Example

```
java -jar dist/bammatrix.jar -o out.png \
	--kg "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" \
	-r "chr1:2345-6789" -B cnv.bed --name BX \
	NOVASEQ/Sample/outs/phased_possorted_bam.bam
```

input with two bam files (for comparing mappers)
```
java  -jar dist/bammatrix.jar -r "chr1:234-567" -o out.png sample.markdup.01.bam sample.markdup.02.bam
```

https://twitter.com/yokofakun/status/1142088565326843904

![https://twitter.com/yokofakun/status/1142088565326843904](https://pbs.twimg.com/media/D9mDYo4WsAAOaSK.jpg)


https://twitter.com/yokofakun/status/1038060108373286912

![https://twitter.com/yokofakun/status/1038060108373286912](https://pbs.twimg.com/media/Dmft0cSXoAAp78l.jpg)

END_DOC
*/

@Program(name="bammatrix",
description="Bam matrix, inspired from 10x/loupe ",
keywords={"sam","bam","compare","matrix"},
creationDate="20190620",
modificationDate="20211206",
jvarkit_amalgamion = true,
menu="CNV/SV"
)
public class BamMatrix  extends Launcher
	{
	private static final Logger LOG = Logger.of(BamMatrix.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"-s","--size"},description="matrix size in pixel")
	private int matrix_size = 1_000;
	@Parameter(names={"-r","-r1","--region"},description="first region." + IntervalParser.OPT_DESC,required=true)
	private String region1Str=null;
	@Parameter(names={"-r2","--region2"},description="2nd region. Default: use first region. " + IntervalParser.OPT_DESC)
	private String region2Str=null;
	@Parameter(names={"--name","-name"},description="user read name or use 'BX:Z:'/'MI:i:' attribute from 10x genomics  as the read name. \"Chromium barcode sequence that is error-corrected and confirmed against a list of known-good barcode sequences.\". See https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam")
	private NameExtractor nameExtractor = NameExtractor.READ_NAME;
	@Parameter(names={"-sa","--sa"},description="Use other canonical alignements from the 'SA:Z:*' attribute")
	private boolean use_sa_align = false;
	@Parameter(names={"-su","--supplementary"},description="Use other supplementary alignements")
	private boolean use_suppl_align = false;
	@Parameter(names={"--color-scale"},description="Color scale")
	private ColorScale color_scale = ColorScale.LOG;
	@Parameter(names={"--mapq"},description="minimal mapping quality")
	private int min_mapq = 30;
	@Parameter(names={"--gtf","-g"},description="Optional gtf file to draw the exons. "+GtfReader.OPT_DESC)
	private Path gtfPath = null;
	@Parameter(names={"--higligth","-B"},description="Optional Bed file to hightlight regions of interest")
	private String highlightPath = null;
	@Parameter(names={"-d","--distance"},description="Don't evaluate a point if the distance between the regions is lower than 'd'. Negative: don't consider distance.",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_distance = -1;
	@Parameter(names={"-C","--min-common"},description="Don't print a point if there are less than 'c' common names at the intersection")
	private int min_common_names = 0;
	@Parameter(names={"--no-coverage"},description="Don't print coverage")
	private boolean hide_coverage = false;
	@Parameter(names={"--pixel"},description="pixel size. Each dot at intersection will have the following size")
	private int pixel_size = 1;

	
	
	/* actual SamReader */
	private SamReader samReaderX = null;
	private SamReader samReaderY = null;
	/* actual sam dict */
	private SAMSequenceDictionary dict;
	/* user interval X axis */
	private SimpleInterval userIntervalX = null;
	/* user interval Y axis */
	private SimpleInterval userIntervalY = null;
		
	private enum NameExtractor {
		READ_NAME,BX,MI;
		String getName(final SAMRecord rec) {
			switch(this) {
				case READ_NAME :  return rec.getReadName();
				case BX: return rec.hasAttribute("BX")?String.class.cast(rec.getAttribute("BX")):null;
				case MI: return rec.hasAttribute("MI")?String.valueOf(rec.getAttribute("MI")):null;
				default: throw new IllegalStateException(this.name());
				}
			}
		};
	
	private enum ColorScale {
		LINEAR,LOG
	}
	
	
	/**  get the reads names in a given interval */
	private abstract class ReadCounter
		{
		protected ReadCounter() throws IOException {}
		/** return the names of the Read names in the interval */
		protected abstract Set<String> getNamesMatching(int side,final Locatable r) throws IOException;
		void dispose() throws IOException {}
		}
	
	/** in memory implementation of ReadCounter: the fastest */
	private class MemoryReadCounter extends ReadCounter
		{
		private final IntervalTreeMap<List<Interval>> treeMapX  = new IntervalTreeMap<>();
		private final IntervalTreeMap<List<Interval>> treeMapY;
		MemoryReadCounter()  throws IOException{
			this.treeMapY = (samReaderY==samReaderX?treeMapX:new IntervalTreeMap<>());
			for(int side=0;side<2;++side) {
				@SuppressWarnings("resource")
				final SamReader sr  =(side==0?samReaderX:samReaderY);
				final IntervalTreeMap<List<Interval>> treeMap  =(side==0?treeMapX:treeMapY);
				final QueryInterval[] qArray;
				
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(sr.getFileHeader());
				if(samReaderY==samReaderX && side==0) {
					final QueryInterval qIntervalX = new QueryInterval(dict.getSequenceIndex(userIntervalX.getContig()), userIntervalX.getStart(), userIntervalX.getEnd());
					final QueryInterval qIntervalY = new QueryInterval(dict.getSequenceIndex(userIntervalY.getContig()), userIntervalY.getStart(), userIntervalY.getEnd());
					qArray = QueryInterval.optimizeIntervals(new QueryInterval[] {qIntervalX,qIntervalY});
					}
				else if(side==0) {
					final QueryInterval qIntervalX = new QueryInterval(dict.getSequenceIndex(userIntervalX.getContig()), userIntervalX.getStart(), userIntervalX.getEnd());
					qArray = QueryInterval.optimizeIntervals(new QueryInterval[] {qIntervalX});
					}
				else if(side==1) {
					final QueryInterval qIntervalY = new QueryInterval(dict.getSequenceIndex(userIntervalY.getContig()), userIntervalY.getStart(), userIntervalY.getEnd());
					qArray = QueryInterval.optimizeIntervals(new QueryInterval[] {qIntervalY});
					}
				else{
					break;
					}
				try(final SAMRecordIterator iter= sr.query(qArray, false))
					{
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						for(final Interval r: BamMatrix.this.samRecordToIntervals(rec)) {
							List<Interval> list = treeMap.get(r);
							if(list==null) {
								list=new ArrayList<>();
								treeMap.put(r,list);
								}
							list.add(r);
							}	
						}
					}
				}
			}

		@Override
		protected Set<String> getNamesMatching(final int side,final Locatable r) throws IOException
			{
			if (!(side==0 || side==1)) throw new IllegalArgumentException("not 0/1");
			final IntervalTreeMap<List<Interval>> treeMap  =(side==0?treeMapX:treeMapY);
			return treeMap.getOverlapping(r).
					stream().
					flatMap(L->L.stream()).
					map(R->R.getName()).
					collect(Collectors.toCollection(()->new HashSet<>(10_000)));
			}

		}
	
	
	
	private boolean accept(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return false;
		if(rec.getReadFailsVendorQualityCheckFlag())  return false;
		if(!this.use_suppl_align && rec.getSupplementaryAlignmentFlag()) return false;
		if(rec.isSecondaryAlignment())  return false;
		if(rec.getMappingQuality()  < this.min_mapq) return false;
		return true;
	}
	
	private List<Interval> samRecordToIntervals(final SAMRecord rec) {
		if(!accept(rec)) return Collections.emptyList();
		
		final String name = this.nameExtractor.getName(rec);
		if(StringUtils.isBlank(name)) return Collections.emptyList();
		final List<Interval> list=new ArrayList<>();
		
		for(final AlignmentBlock ab:rec.getAlignmentBlocks())
			{
			final Interval r = new Interval(
					rec.getReferenceName(),
					ab.getReferenceStart(),
					ab.getReferenceStart()+ab.getLength(),
					rec.getReadNegativeStrandFlag(),
					name
					);
			list.add(r);
			}
		if( use_sa_align) {
			for(final SAMRecord rec2:SAMUtils.getOtherCanonicalAlignments(rec))
				{
				if(!(userIntervalX.overlaps(rec2) || userIntervalY.overlaps(rec2)))
					{
					continue;
					}
											
				final Interval r = new Interval(
						rec2.getReferenceName(),
						rec2.getStart(),
						rec2.getStart(),
						rec2.getReadNegativeStrandFlag(),
						name
						);
				list.add(r);
				}
			}
		return list;
		}
	
	private boolean validateDisance(final Interval r1,final Interval r2) {
		if(min_distance<0) return true;
		if(r1.withinDistanceOf(r2, this.min_distance)) return false;
		return true;
		
	}
	
	@Override
	public int doWork(final List<String> args) {
		if(pixel_size<1) {
			LOG.error("pixel size is too small ("+this.pixel_size+")");
			return -1;
			}
		if(StringUtils.isBlank(region2Str)) {
			this.region2Str = region1Str;
		}
		try {
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT)
					;
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
						
			final String inputX;
			final String inputY;
			if(args.size()==1) {
				inputX = args.get(0);
				inputY = null;
				}
			else if(args.size()==2) {
				inputX = args.get(0);
				inputY = args.get(1);
				}
			else {
				LOG.error("illegal number of arguments.");
				return -1;
				}
			this.samReaderX = srf.open(SamInputResource.of(inputX));
			if(!this.samReaderX.hasIndex()) {
				LOG.error("Input "+inputX+" is not indexed");
				return -1;
				}
			this.dict = SequenceDictionaryUtils.extractRequired(this.samReaderX.getFileHeader());
			
			if(inputY!=null) {
				this.samReaderY = srf.open(SamInputResource.of(inputY));
				if(!this.samReaderY.hasIndex()) {
					LOG.error("Input "+inputY+" is not indexed");
					return -1;
					}
				SequenceUtil.assertSequenceDictionariesEqual(
					SequenceDictionaryUtils.extractRequired(this.samReaderY.getFileHeader()),
					this.dict);
				}
			else
				{
				this.samReaderY = this.samReaderX;
				}
			final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(this.dict);
			
			final Function<String,Optional<SimpleInterval>> intervalParser = 
					new IntervalParser(dict).
					enableWholeContig();
			this.userIntervalX = intervalParser.apply(this.region1Str).orElseThrow(IntervalParser.exception(this.region1Str));
			this.userIntervalY = intervalParser.apply(this.region2Str).orElseThrow(IntervalParser.exception(this.region2Str));
			
			
			// adjust intervals so they have the same length
			if(this.userIntervalX.getLengthOnReference() > this.userIntervalY.getLengthOnReference()) {
				final int mid =  this.userIntervalY.getStart()+ this.userIntervalY.getLengthOnReference()/2;
				final int start = Math.max(1,mid- this.userIntervalX.getLengthOnReference()/2);
				
				this.userIntervalY = new SimpleInterval(
						this.userIntervalY.getContig(),
						start,
						start + this.userIntervalX.getLengthOnReference()
						);
				LOG.warn("Adjusting interval Y to "+this.userIntervalY+" so both intervals have the same length");
				}
			else if(this.userIntervalY.getLengthOnReference() > this.userIntervalX.getLengthOnReference()) {
				final int mid =  this.userIntervalX.getStart()+ this.userIntervalX.getLengthOnReference()/2;
				final int start = Math.max(1,mid- this.userIntervalY.getLengthOnReference()/2);

				this.userIntervalX = new SimpleInterval(
						this.userIntervalX.getContig(),
						start,
						start + this.userIntervalY.getLengthOnReference()
						);
				LOG.warn("Adjusting interval X to "+this.userIntervalX+" so both intervals have the same length");
				}
			LOG.info("One pixel is "+(this.userIntervalX.getLengthOnReference()/(double)matrix_size)+" bases");
			
			
			
			final int distance= Math.max(
					this.userIntervalX.getLengthOnReference(),
					this.userIntervalY.getLengthOnReference()
					);
			final double pixel2base = distance/(double)matrix_size;
			short max_count=1;
			final short counts[]=new short[this.matrix_size*this.matrix_size];

			
			final ReadCounter counter  = new MemoryReadCounter();
			
			/* loop over each pixel 1st axis */
			for(int pixY=0;pixY< this.matrix_size;pixY++)
				{				
				final int start1 = (int)(this.userIntervalY.getStart() + pixY * pixel2base);
				final int end1 = start1 + (int)pixel2base;
				final Interval qy = new Interval(this.userIntervalY.getContig(), start1, end1);
				if(!qy.overlaps(this.userIntervalY)) continue;
				final Set<String> set1 = counter.getNamesMatching(1,qy);
				if(set1.isEmpty()) continue;
				
				/* loop over each pixel 2nd axis */
				for(int pixX=0;pixX< this.matrix_size;pixX++)
					{
					final int start2 = (int)(this.userIntervalX.getStart() + pixX* pixel2base);
					final int end2 = start2 + (int)pixel2base;
					final Interval qx = new Interval(this.userIntervalX.getContig(), start2, end2);
					if(!qx.overlaps(this.userIntervalX)) continue;
					if(!validateDisance(qy,qx)) continue;
					
					
					final int count_common;
					if(qx.compareTo(qy)==0) {
						count_common = set1.size();
						}
					else
						{
						final HashSet<String> common = new HashSet<>(set1);
						common.retainAll(counter.getNamesMatching(0,qx));
						count_common = common.size();
						}
					final short count =  count_common>Short.MAX_VALUE?Short.MAX_VALUE:(short)count_common;
					max_count = (short)Math.max(count, max_count);
					counts[pixY*this.matrix_size+pixX] = count;
					}
				}
			counter.dispose();
			final int font_size=10;
			final int cov_height = (this.hide_coverage?0:50);
			final int gene_height = 25;
			final int margin = font_size+ cov_height+ (this.gtfPath==null?0:gene_height);
			final Insets margins = new Insets(margin,margin, 10, 10);
			
			final Dimension drawingAreaDim = new Dimension(
					this.matrix_size+margins.left+margins.right,
					this.matrix_size+margins.top+margins.bottom
					);
			
			final BufferedImage img = new BufferedImage(
					drawingAreaDim.width,
					drawingAreaDim.height,
					BufferedImage.TYPE_INT_RGB
					);
			final Graphics2D g = img.createGraphics();
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, drawingAreaDim.width,drawingAreaDim.height);
			
			// draw sample
			final Hershey herschey = new Hershey();

			final String sampleX = samReaderX.getFileHeader().getReadGroups().stream().map(R->R.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse(inputX);
			final String sampleY = (samReaderX==samReaderY?sampleX:samReaderX.getFileHeader().getReadGroups().stream().map(R->R.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse(inputY));
			final String sample = (sampleX.equals(sampleY)?sampleX:String.join(" ", sampleX,sampleY));
			
			g.setColor(Color.DARK_GRAY);
			herschey.paint(g, sample, new Rectangle2D.Double(0,1,margins.left-1,font_size));
			
			
			// draw highlight regions
			
			for(int side=0;side< 2 && !StringUtils.isBlank(this.highlightPath);++side) {
				final int curr_side=side;
				final SimpleInterval r = (side==0?this.userIntervalX:this.userIntervalY);
				final BedLineCodec bedCodec = new BedLineCodec();
				final Composite oldComposite = g.getComposite();
				g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.3f));
					try(BufferedReader br = IOUtils.openURIForBufferedReading(this.highlightPath)) {
						br.lines().
							filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
							map(L->bedCodec.decode(L)).
							filter(B->B!=null).
							filter(K->converter.apply(K.getContig())!=null && r.getContig().equals(converter.apply(K.getContig()))).
							filter(K->CoordMath.overlaps(K.getStart(), K.getEnd(), r.getStart(), r.getEnd())).
							map(E->new Interval(converter.apply(E.getContig()),E.getStart()+1,E.getEnd())).
							filter(E->CoordMath.overlaps(E.getStart(), E.getEnd(), r.getStart(), r.getEnd())).
							map(E->new Interval(E.getContig(),Math.max(r.getStart(),E.getStart()),Math.min(r.getEnd(),E.getEnd()))).
							forEach(E->{
								double d = ((E.getStart()-r.getStart())/(double)r.getLengthOnReference())*matrix_size;
								double dL  = ((E.getLengthOnReference())/(double)r.getLengthOnReference())*matrix_size;
								g.setColor(Color.YELLOW);
								if(curr_side==0)
									{
									g.fill(new Rectangle2D.Double(d, 0, dL, margins.left));
									}
								else
									{
									g.fill(new Rectangle2D.Double(0, d, margins.top,dL));
									}
								
							});
					}
				g.setComposite(oldComposite);
			}
			
			
			g.translate(margins.left, margins.top);
			
			final double logMaxV =Math.log(max_count); 
			for(int pix1=0;pix1< this.matrix_size;pix1++)
				{
				for(int pix2=0;pix2< this.matrix_size;pix2++)
					{
					final short count = counts[pix1*this.matrix_size+pix2];
					if(count==0 || count < this.min_common_names) continue;
					final int gray;
					switch(color_scale) {
						case LINEAR:
							gray = 255-(int)(255*(count/(double)max_count));
							break;
						case LOG:
							gray = 255-(int)(255*((Math.log(count))/logMaxV));
							break;
						default: throw new IllegalStateException(color_scale.name());
						}
					g.setColor(new Color(gray,0,0));
					g.fill(new Rectangle2D.Double(
						pix1 - pixel_size/2.0,
						pix2 - pixel_size/2.0,
						pixel_size,
						pixel_size
						));
					}
				}
			// draw frame
			g.setColor(Color.GRAY);
			g.drawRect(0, 0, this.matrix_size, this.matrix_size);

			g.translate(-margins.left,-margins.top);
			

			// used to plot depth
			final double coverage[]=new double[matrix_size];
			
			final List<SimpleInterval> exonsList;
			if(this.gtfPath==null) {
				exonsList = Collections.emptyList();
			} else
				{
				try(GtfReader gtfReader = new GtfReader(this.gtfPath)) {
					gtfReader.setContigNameConverter(converter);
					exonsList = gtfReader.
						getAllGenes().
						stream().
						filter(K->K.overlaps(this.userIntervalX) || K.overlaps(this.userIntervalY)).
						flatMap(G->G.getTranscripts().stream()).
						filter(T->T.hasExon()).
						flatMap(K->K.getExons().stream()).
						filter(E->E.overlaps(this.userIntervalX) || E.overlaps(this.userIntervalY)).
						map(E->new SimpleInterval(E)).
						collect(Collectors.toSet()).
						stream().
						collect(Collectors.toList());
					}
				}
			
			
			for(int side=0;side< 2;++side) {
				final SimpleInterval r = (side==0?this.userIntervalX:this.userIntervalY);
				final AffineTransform oldtr = g.getTransform();
				AffineTransform tr;
				if(side==0) {
					//horizonal axis
					tr = AffineTransform.getTranslateInstance(margins.left,1);
				} else {
					// vertical
					tr=  AffineTransform.getTranslateInstance(margins.left,margins.top);
					tr.concatenate(AffineTransform.getRotateInstance(Math.PI/2.0));
				}
				g.setTransform(tr);
				
				// calculate coverage , do this only once if regionX==regionY
				if(!hide_coverage && !(side==1 && this.userIntervalX.equals(this.userIntervalY))) {
					Arrays.fill(coverage,0);
					final int count[]=new int[this.matrix_size];
					
					final IntervalList intervalList = new IntervalList(this.dict);
					intervalList.add(new Interval(r));
					try(final SamLocusIterator sli = new SamLocusIterator(this.samReaderX,intervalList,true)) {
					while(sli.hasNext()) {
						final LocusInfo locusInfo = sli.next();
						final int pos = locusInfo.getPosition();
						if(pos < r.getStart() || pos > r.getEnd()) continue;
						
						final int depth = locusInfo.getRecordAndOffsets().size();
						final int array_index = (int)(((pos-r.getStart())/(double)r.getLengthOnReference())*matrix_size);
						
						coverage[array_index] += depth;
						count[array_index]++;
						}
					}
					
					
					for(int i=0;i< coverage.length;++i) {
						if(count[i]==0) continue;
						coverage[i]/=count[i];
						}
					}
				// draw ruler
				int y = 0;
				
				if(!this.hide_coverage) {
					final double max_cov = Arrays.stream(coverage).max().orElse(1);
					final GeneralPath gp =new GeneralPath();
					gp.moveTo(0, cov_height);
					for(int x=0;x < coverage.length;++x) {
						gp.lineTo(x,y + cov_height-(coverage[x]/max_cov)*cov_height);	
						}
					gp.lineTo(coverage.length, cov_height);
					gp.closePath();
					g.setColor(Color.GRAY);
					g.fill(gp);
					// string for max cov
					String label = StringUtils.niceInt( (int)Arrays.stream(coverage).max().orElse(9));
					g.setColor(Color.DARK_GRAY);
					herschey.paint(g,label,
						new Rectangle2D.Double(
						matrix_size - label.length()*font_size,
						y,
						label.length()*font_size,
						font_size));
					
					y+= cov_height;
					}
				
				//draw label
				g.setColor(Color.DARK_GRAY);
				
				// label is 'start position'
				String label = StringUtils.niceInt(r.getStart());
				herschey.paint(
						g,
						label,
						new Rectangle2D.Double(0, y, label.length()*font_size, font_size)
						);
				
				// label is 'end position'
				label = StringUtils.niceInt(r.getEnd());
				herschey.paint(
						g,
						label,
						new Rectangle2D.Double(matrix_size-(label.length()*font_size),y, label.length()*font_size, font_size)
						);
				
				// label is 'chromosome and length'
				label = r.getContig()+" ( "+ StringUtils.niceInt(r.getLengthOnReference())+" bp )";
				herschey.paint(
						g,
						label,
						new Rectangle2D.Double(matrix_size/2.0-(label.length()*font_size)/2.0,y, label.length()*font_size, font_size)
						);
				
				y+= font_size;
				
				// draw genes
				if(this.gtfPath!=null) {
					final double curr_y = y;
					double midy = y+ gene_height/2.0;
					g.setColor(Color.CYAN);
					g.draw(new Line2D.Double(0,midy,matrix_size,midy));
					exonsList.
							stream().
							filter(E->E.overlaps(r)).
							map(E->new SimpleInterval(E.getContig(),Math.max(r.getStart(),E.getStart()),Math.min(r.getEnd(),E.getEnd()))).
							forEach(E->{
								final double x = ((E.getStart()-r.getStart())/(double)r.getLengthOnReference())*matrix_size;
								final double width  = ((E.getLengthOnReference())/(double)r.getLengthOnReference())*matrix_size;
								g.setColor(Color.BLUE);
								g.fill(new Rectangle2D.Double(x, curr_y, width, gene_height));
							});
					}
				
				
				g.setTransform(oldtr);
				}
			
			
			g.dispose();
			try {
				if(this.outputFile==null)
					{
					ImageIO.write(img,"PNG",stdout());
					}
				else
					{
					ImageIO.write(img,this.outputFile.getName().endsWith(".png")?"PNG":"JPG", this.outputFile);
					}
				} catch(final IOException err)
				{
					throw new RuntimeIOException(err);
				}
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.samReaderX);
			CloserUtil.close(this.samReaderY);
			}
		}
	public static void main(final String[] args) {
		new BamMatrix().instanceMainWithExit(args);
		}
}

