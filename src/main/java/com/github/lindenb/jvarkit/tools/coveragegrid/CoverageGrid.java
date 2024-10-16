/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.coveragegrid;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Arc2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.canvas.Canvas;
import com.github.lindenb.jvarkit.canvas.CanvasFactory;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.SampleSheet;
import com.github.lindenb.jvarkit.io.SampleSheetFactory;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.Average;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.math.MinMaxInteger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Algorithm;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC
## input

input is a tab-delimited samplesheet with the following columns:

| column | required ? | description |
|--------|------------|-------------|
| bam | required | /path/to/indexed/bam+or+cram |
| color | optional | color used for frame (could be a hint to spot cases/controls). Can be empty. Otherwise, use default color |
| sample | optional | sample name. If empty the sample will be extracted from the bam. If it starts with '+=' , the name will appended to the original bam name  |



## example:

```
find dir -type f -name "*bam" > in.list 
java -jar dist/jvarkit.jar coveragegrid -R src/test/resources/rotavirus_rf.fa --region "RF01:100-200" in.list
```


END_DOC 
 */
@Program(
	name="coveragegrid",
	description="Display an image of depth to display any anomaly an intervals+bams as a grid image",
	keywords={"cnv","bam","depth","coverage","svg","postscript"},
	creationDate="20241009",
	modificationDate="20241009",
	jvarkit_amalgamion =  true,
	menu="CNV/SV"
	)
public class CoverageGrid extends Launcher {
	private enum PlotType {COVERAGE,MEDIAN_COVERAGE,PILEUP,PILEUP_PAIR,GRID,DISCORDANT,SUPPL};
	private static final Logger LOG = Logger.build( CoverageGrid.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--regions","--region","--interval"},description = "Interval region",required=true)
	private String intervalStr=null;
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--dimension","--dim"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--extend","-x"},description = "extends the interval x times")
	private double extendFactor= 3.0;
	@Parameter(names= {"--max-y"},description="Max normalized Y, when using --type MEDIAN_COVERAGE")
	private double max_normalized_y = 3.0;
	@Parameter(names={"--threads"},description="number of threads")
	private int nThreads = 1;
	@Parameter(names= {"--format"},description=CanvasFactory.OPT_FORMAT_DESC)
	private CanvasFactory.Format canvas_format=CanvasFactory.Format.SVG;
	@Parameter(names= {"--title"},description="Set Title for graphics")
	private String title="";
	@Parameter(names= {"--vcf"},description="indexed VCF file to show variants")
	private Path vcfFile=null;
	@Parameter(names= {"--gtf"},description="indexed GTF file to show genes")
	private String gtfFile=null;
	@Parameter(names= {"--inversion-flag"},description="for DISCORDANT or SUPPL : arcs as a edge outside interval and the other inside the interval")
	private boolean inversion_flag=false;

	
	
	@Parameter(names= {"--type","--what","--plot"},description="Plot type. "
			+ "COVERAGE: DEPTH of coverage, "
			+ "MEDIAN_COVERAGE: coverage normalize by external boundaries, "
			+ "PILEUP: show reads, "
			+ "PILEUP_PAIR: show paired-end fragments, "
			+ "GRID: matrix of read name co-occurence. Slow and memory consuming.,"
			+ "DISCORDANT: plot discordant reads as arcs overlaping edges, "
			+ "SUPPL: plot sipplementary read pairs as arcs overlaping edges"
			)
	private PlotType plotType = PlotType.MEDIAN_COVERAGE;
	@DynamicParameter(names = "-D", description = "other parameters '-Dkey=value'. "+
			"-Dhide_insertions=true "+
			"-Dinversion_mode=true "+
			"-Doverlap_boundaries=true "
			)
	private Map<String, String> dynaParams = new HashMap<>();

	
	private final Color EXTERNAL_AREA_COLOR =new Color(230,230,230);

	
	private abstract class BamTask implements Runnable {
		/** name in BAM */
		String bamSampleName="";
		/** use provided sample name */
		String sampleName="";
		Color sampleColor = null;
		Path bamPath;
		Path referencePath;
		Locatable userInterval;
		Locatable extendedInterval;
		Rectangle boundaries;
		Path vcfPath;
		String gtfPath;
		byte[] reference_bases;
		Canvas canvas;
		SAMSequenceDictionary bamSequenceDictionary = null;
		
		protected BamTask() {
			}
		
		protected SamReader openSamReader() throws IOException {
			final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					disable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES).
					referenceSequence(this.referencePath).
					validationStringency(ValidationStringency.LENIENT);
			final SamReader sr= samReaderFactory.open(this.bamPath);
			if(!sr.hasIndex()) {
				throw new IOException("index missng for "+bamPath);
				}
			final SAMFileHeader hdr = sr.getFileHeader();
			this.bamSequenceDictionary  = hdr.getSequenceDictionary();
			this.bamSampleName = hdr.getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtils.isBlank(S)).
					findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(this.bamPath))
					;
			if(StringUtils.isBlank(this.sampleName)) {
				this.sampleName = this.bamSampleName ;
				}
			else if(this.sampleName.startsWith("+=")) {
				this.sampleName = this.bamSampleName  + " "+this.sampleName.substring(2).trim();
				}
			return sr;
			}
		protected double trimX(double coordX) {
			return Math.min(this.boundaries.getMaxX(), Math.max(this.boundaries.getX(),coordX));
			}
		
		protected double position2pixel(final int pos) {
			return this.boundaries.getX() +
					((pos-this.extendedInterval.getStart())/(double)this.extendedInterval.getLengthOnReference())*this.boundaries.getWidth();
			}
		
		protected double getGeneMidY() {
			return (float)Math.max(this.boundaries.getY(),this.boundaries.getMaxY()-4);
		}
		
		protected void plotGenes(final Canvas canvas) {
			if(this.gtfPath==null) return;
			final float midy=(float)getGeneMidY();
			final GTFCodec codec=new GTFCodec();
			try(TabixReader tabix = new TabixReader(this.gtfPath)) {
				for(int step=0;step< 3;++step) {
					String expect;
					float featHeight;
					Color color;
					switch(step) {
						case 0: expect="gene";featHeight=0.1f;color=Color.MAGENTA;break;
						case 1: expect="exon";featHeight=0.4f;color=Color.MAGENTA;break;
						case 2: expect="CDS";featHeight=0.6f;color=Color.MAGENTA;break;
						default: expect=null;featHeight=0f;color=Color.MAGENTA;break;
						}
			
					final List<Rectangle2D> rects = new ArrayList<>();
					TabixReader.Iterator iter= tabix.query(
							this.extendedInterval.getContig(),
							this.extendedInterval.getStart(),
							this.extendedInterval.getEnd());
					for(;;) {
						final String line = iter.next();
						if(line==null) break;
						final GTFLine feat=codec.decode(line);
						if(feat==null) continue;
						if(!feat.getType().equals(expect)) continue;
						final double x1= trimX(position2pixel(feat.getStart()));
						final double x2= trimX(position2pixel(feat.getEnd()+1));
						if(x1>=x2) continue;
						if(step>0 /* not a gene */ && (x2-x1)<2 /* too small to be displayed */) continue;
						rects.add(new Rectangle2D.Float((float)x1, midy-featHeight/2f, (float)(x2-x1), featHeight));
						}
					Collections.sort(rects,(A,B)->Double.compare(A.getX(),B.getX()));
					int k=0;
					while(k+1< rects.size()) {
						final Rectangle2D r1 = rects.get(k+0);
						final Rectangle2D r2 = rects.get(k+1);
						if(r2.getX() - r1.getMaxX() < 0.1) {
							final float x1 = (float)Math.min(r1.getX(), r2.getX());
							final float x2 = (float)Math.max(r1.getMaxX(), r2.getMaxX());
							rects.set(k,new Rectangle2D.Float(
								x1,
								(float)r1.getY(),
								(x2-x1),
								(float)r1.getHeight()
								));
							rects.remove(k+1);
							}
						else
							{
							k++;
							}
						}
					for(Rectangle2D rect:rects) {
						// draw feature
						canvas.shape(
							rect,
							FunctionalMap.of(
									Canvas.KEY_FILL, color,
									Canvas.KEY_STROKE,null
									)
							);
						}
					}
				}
			catch(Exception err) {
				LOG.error(err);
				}
			}
		
		protected void plotVariants(final Canvas canvas) {
			if(this.vcfPath==null) return;
			double prev_x=-10000;
			try(VCFReader reader = new VCFFileReader(this.vcfPath, true)) {
				final VCFHeader header= reader.getHeader();
				final Hyperlink hlink =Hyperlink.compile(header.getSequenceDictionary());
				try(CloseableIterator<VariantContext> iter = reader.query(this.extendedInterval)) {
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();

						if(!ctx.isVariant()) continue;

						if(CoordMath.encloses(
							ctx.getStart(), ctx.getEnd(),
							this.extendedInterval.getStart(), this.extendedInterval.getEnd()
							)) {
							continue;
							}
						if(ctx.hasGenotypes() && !StringUtils.isBlank(this.sampleName)) {
							final Genotype g = ctx.getGenotype(this.sampleName);
							if(g==null) continue;//no read group in bam header
							if(!g.hasAltAllele()) continue;
							}
						final double x1 = trimX(position2pixel(ctx.getStart()));
						final double x2 = trimX(position2pixel(ctx.getEnd()+1));
						
						
						if(x1>=x2 || (x2-x1) < 0.01) {
							continue;
							}
						
						if((int)prev_x==(int)x1) {
							continue;
							}
						
						prev_x=x1;
						// draw variant
						canvas.rect(
							x1,
							boundaries.getY(),
							(x2-x1),
							boundaries.getHeight(),
							FunctionalMap.of(
									Canvas.KEY_FILL,Color.MAGENTA,
									Canvas.KEY_STROKE,null,
									Canvas.KEY_HREF, hlink.apply(ctx).orElse(null),
									Canvas.KEY_TITLE, ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString()
									)
							);
						
						}
					}
				}
			catch(Throwable err) {
				LOG.error(err);
				}
			}
		
		
		protected CloseableIterator<SAMRecord> query(final SamReader sr) {
			final int extend = 100;//give a chance to get clipping
			return sr.query(
					this.extendedInterval.getContig(),
					Math.max(1,this.extendedInterval.getStart() - extend),
					this.extendedInterval.getEnd() + extend, 
					false
					);
			}
		
		protected boolean acceptRead(final SAMRecord rec) {
			if(!SAMRecordDefaultFilter.accept(rec, CoverageGrid.this.min_mapq)) return false;
			
			return true;
			}
		
		protected void plotVerticalGuides(final Canvas canvas) {
			// draw original region vertical bounds 
			if(!this.userInterval.equals(this.extendedInterval)) {
				for( int i=0;i<2;i++) {
					final int p = (i==0?userInterval.getStart():userInterval.getEnd()+1);
					final double x = trimX(position2pixel(p));
					canvas.line(
							x,
							this.boundaries.getY(),
							x,
							this.boundaries.getMaxY(),
							FunctionalMap.of(
								Canvas.KEY_STROKE,Color.CYAN,
								Canvas.KEY_TITLE, StringUtils.niceInt(p),
								Canvas.KEY_STROKE_WIDTH,0.5
								)
							);
					}
				}
			}
		
		protected void plotFrame(final Canvas canvas) {
			// draw frame
			canvas.rect(
				boundaries.getX(),
				boundaries.getY(),
				boundaries.getWidth(),
				boundaries.getHeight(),
				FunctionalMap.of(
						Canvas.KEY_FILL,null,
						Canvas.KEY_STROKE,(sampleColor==null?Color.DARK_GRAY:sampleColor),
						Canvas.KEY_STROKE_WIDTH,0.5
						)
				);
			}
		protected void plotBackgroundFrame(Canvas canvas) {
			
			for(int i=0;i< 2;i++) {
				double x1 = trimX(position2pixel(i==0?this.extendedInterval.getStart():this.userInterval.getEnd()+1));
				double x2 = trimX(position2pixel(i==0?this.userInterval.getStart():this.extendedInterval.getEnd()+1));
				if(x1>=x2) continue;
				// draw frame
				canvas.rect(
					x1,
					boundaries.getY(),
					(x2-x1),
					boundaries.getHeight(),
					FunctionalMap.of(
						Canvas.KEY_FILL,EXTERNAL_AREA_COLOR,
						Canvas.KEY_STROKE,null
						)
					);
				}
			}
		
		abstract boolean loadData();
		abstract void disposeData();
		abstract void plot(Canvas canvas);
		@Override
		public final void run() {
			System.err.println("submit "+this.bamPath.getFileName().toString());
			try {
				if(loadData()) {
					synchronized (this.canvas) {
						plot(this.canvas);
						}
					}
				}
			catch(Throwable err) {
				LOG.error(err);
				}
			finally {
				disposeData();
				}
			}
		}	

	/***************************************************************************************************************************/
	/***************************************************************************************************************************/
	/***************************************************************************************************************************/

	private class PileupTask extends BamTask {
		final boolean by_pair;
		PileupTask(boolean by_pair) {
			this.by_pair  = by_pair;
			}
		
		/** simplified read */
		private  class MiniRead implements Locatable {
			final int start;
			final int end;
			final Cigar cigar;
			final int flags;
			Set<Integer> mistmatches = null;
			MiniRead(int start,int end,Cigar cigar,int flags) {
				this.start= start;
				this.end = end;
				this.cigar= cigar;
				this.flags = flags;
				}
			MiniRead(final SAMRecord rec) {
				this(rec.getAlignmentStart(),rec.getAlignmentEnd(),rec.getCigar(),rec.getFlags());
				}
			@Override
			public String getContig() {
				return PileupTask.this.extendedInterval.getContig();
				}
			@Override
			public int getStart() {
				return start;
				}
			@Override
			public int getEnd() {
				return end;
				}
			public Cigar getCigar() {
				return cigar;
				}
			public int getUnclippedStart() {
				return SAMUtils.getUnclippedStart(getStart(), getCigar());
				}
			public int getUnclippedEnd() {
				return SAMUtils.getUnclippedEnd(getEnd(), getCigar());
				}
			private boolean testFlag(SAMFlag flg) {
				return (this.flags & flg.intValue()) != 0;
				}
			public boolean isReadPaired() {
		        return testFlag(SAMFlag.READ_PAIRED);
				}
			}
		
		private class MiniReadPair implements Locatable,Comparable<MiniReadPair>,Iterable<MiniRead> {
			private final MiniRead R1 ;
			private MiniRead R2 = null;
			
			public MiniReadPair(final MiniRead rec)
				{
				this.R1 = rec;
				}
			
			
			@Override
			public int compareTo(final MiniReadPair o) {
				int i= Integer.compare(this.getStart(),o.getStart());
				if(i!=0) return i;
				i= Integer.compare(this.getEnd(),o.getEnd());
				if(i!=0) return i;
				return 0;
				}
			
			@Override
			public String getContig() {
				return R1.getContig();
				}
			@Override
			public int getStart() {
				return R2==null?
					R1.getUnclippedStart():
					Math.min(R1.getUnclippedStart(), R2.getUnclippedStart())
					;
				}
			@Override
			public int getEnd() {
				return R2==null?
						R1.getUnclippedEnd():
						Math.max(R1.getUnclippedEnd(), R2.getUnclippedEnd())
						;
				}
			@Override
			public Iterator<MiniRead> iterator() {
				return (R2==null?Arrays.asList(R1):Arrays.asList(R1,R2)).iterator();
				}
			boolean isStrangelyPaired() {
				if(!R1.isReadPaired()) return false;
				if(!R1.testFlag(SAMFlag.PROPER_PAIR)) return true;
				if(R1.testFlag(SAMFlag.READ_REVERSE_STRAND)==R1.testFlag(SAMFlag.MATE_REVERSE_STRAND)) return true;
				return false;
				}
			}

		
		private final List<List<MiniReadPair>> rows = new ArrayList<>();
		
		@Override
		public boolean loadData() {
			final Map<String,MiniReadPair> readName2pairs = new HashMap<>();
			try(SamReader sr=openSamReader()) {
				 try(CloseableIterator<SAMRecord> iter=this.query(sr)) {
					 while(iter.hasNext()) {
						 final SAMRecord rec0=iter.next();
						 if(!acceptRead(rec0)) continue;
						 final Cigar cigar = rec0.getCigar();
						 if(cigar==null || cigar.isEmpty()) continue;
						 
						 final MiniRead rec = new MiniRead(rec0);
						 
						 /** record mismatch ? */
						 if((this.boundaries.getWidth()/(double)this.extendedInterval.getLengthOnReference())>=3 &&
								 rec0.getReadBases()!=null &&
								 rec0.getReadBases()!=SAMRecord.NULL_SEQUENCE) {
							 	int read0=0;
				   	     		int ref1 = rec0.getAlignmentStart();
				   	     		
							 	final byte bases[]=rec0.getReadBases();
					   	     	for(CigarElement ce: cigar.getCigarElements()) {
					    			if(ref1> this.extendedInterval.getEnd()) {
					    				break;
					    				}
					   	     		

					    			final CigarOperator op=ce.getOperator();
					    			switch(op) {
					    			case P:break;
					    			case H:break;
					    			case D: case N: ref1+=ce.getLength(); break;
					    			case S: case I: read0+=ce.getLength(); break;
					    			case EQ:case M: case X:
					    				{
					    				for(int j=0;j< ce.getLength();j++) {
					    					if(ref1+ j <  this.extendedInterval.getStart()) {
					    						continue;
					    						}
					    					if(ref1+ j > this.extendedInterval.getEnd()) {
					    						break;
					    						}
					    					final int ref_base_idx = ref1-this.extendedInterval.getStart()+j;
					    					final char ctgBase =(char)(ref_base_idx<0 || ref_base_idx>= this.reference_bases.length?'N':Character.toUpperCase(this.reference_bases[ref_base_idx]));
					    					
					    					if(ctgBase=='N') {
					    						continue;
					    						}
					    					final char readBase = (char)(read0<0 || read0>= bases.length?'N':Character.toUpperCase(bases[read0+j]));
					    					if(readBase=='N') {
					    						continue;
					    						}

					    					
					    					if(readBase==ctgBase) {
					    						continue;
					    						}
					    					if(rec.mistmatches==null) rec.mistmatches = new HashSet<>();
					    					rec.mistmatches.add(ref1+j);
					    					}
					    				read0+=ce.getLength();
					    				ref1+=ce.getLength();
					    			
					    				break;
					    				}
					    			default:break;
				    				}
				   	     		}
						 	}
						 
						 final String pair_name = this.by_pair ? rec0.getReadName(): rec0.getPairedReadName();
						 MiniReadPair pair=readName2pairs.get(pair_name);
						 if(pair==null) {
							pair =new MiniReadPair(rec);
							if( this.by_pair &&
								rec0.getReadPairedFlag() &&
								!rec0.getMateUnmappedFlag() && 
								rec0.getReferenceName().equals(rec0.getMateReferenceName())) {
								final int mate_start = rec0.getMateAlignmentStart();
								final int mate_end = getMateEnd(rec0);
								final int mate_flags = SAMFlag.READ_PAIRED.intValue()+
										(rec0.getFirstOfPairFlag()?SAMFlag.SECOND_OF_PAIR.intValue():SAMFlag.FIRST_OF_PAIR.intValue())+
										(rec0.getMateNegativeStrandFlag()?SAMFlag.READ_REVERSE_STRAND.intValue():0) +
										(rec0.getReadNegativeStrandFlag()?SAMFlag.MATE_REVERSE_STRAND.intValue():0)  +
										(rec0.getProperPairFlag()?SAMFlag.PROPER_PAIR.intValue():0) 
										;
								final MiniRead mate = new MiniRead(
										mate_start,
										mate_end,
										new Cigar(Arrays.asList(new CigarElement(1+mate_end-mate_start, CigarOperator.M))),
										mate_flags
										);
								pair.R2 = mate;
								}
							readName2pairs.put(pair_name,pair);
						 	}
						 else {
							if(!rec0.getReadPairedFlag()) LOG.error("reads not paired with same name "+pair_name);
							pair.R2 = rec;
						 	}
					 	 }
					 }//end iterator
				 
				 final Pileup<MiniReadPair> pileup = new Pileup<>((L,R)->position2pixel(L.getEnd()+1) +1  < position2pixel(R.getStart()));
				 for(MiniReadPair pair: readName2pairs.values().stream().sorted().collect(Collectors.toList())) {
					 if(!pair.overlaps(extendedInterval)) continue;
					 pileup.add(pair);
				 	}
				
				 this.rows.addAll(pileup.getRows());
				return true;
				}//end samreader
			catch(IOException err) {
				LOG.error(err);
				return false;
				}
			}
		
		@Override
		void plot(final Canvas canvas) {
		     final double featureHeight= Math.min(20,(this.boundaries.getHeight()/Math.max(1.0,(double)this.rows.size())));
			 final boolean hide_insertions=  CoverageGrid.this.dynaParams.getOrDefault("hide_insertions", "false").equals("true");

		     double y= boundaries.getMaxY()-featureHeight;
		     
		     
		     plotBackgroundFrame(canvas);
		     plotVariants(canvas);
		     plotGenes(canvas);
			 plotVerticalGuides(canvas);

		     for(final List<MiniReadPair> row:this.rows) {
		    	final double h2= Math.max(featureHeight*0.8,featureHeight-2);
		    	for(final MiniReadPair pair: row) {		    
		    		final List<Integer> insertions = new ArrayList<>();
		    		/* draw rec itself */
		    		final double midy=y+h2/2.0;
		    		if(pair.R2!=null || 
		    				/* is there a deletion */ pair.R1.getCigar().getCigarElements().
		    				stream().
		    				map(CE->CE.getOperator()).
		    				anyMatch(O->O.equals(CigarOperator.N) || O.equals(CigarOperator.D))
		    				) {
		    			canvas.line(
		    				trimX(position2pixel(pair.getStart())),
		    				midy,
		    				trimX(position2pixel(pair.getEnd()+1)),
		    				midy,
		    				FunctionalMap.of(
	    						Canvas.KEY_STROKE, pair.isStrangelyPaired()? Color.RED:Color.GRAY,
	    						Canvas.KEY_STROKE_WIDTH, 0.1
	    						)
		    				);
			    		}
		    		
			    	for(MiniRead rec:pair) {
			    		int ref1 = rec.getUnclippedStart();
			    		for(final CigarElement ce: rec.getCigar().getCigarElements()) {
			    			if(ref1> this.extendedInterval.getEnd()) break;
			    			final CigarOperator op=ce.getOperator();
			    			Shape shape = null;
			    			Color fill=null;
			    			Color stroke=Color.DARK_GRAY;
			    			switch(op) {
			    				case P: break;
			    				case M://through
			    				case X://through
			    				case EQ://through
			    				case S: //through
			    				case H: 
			    						final double x1= trimX(position2pixel(ref1));
			    						final double x2 =trimX(position2pixel(ref1+ce.getLength()));
			    						if(x1 < x2) {
				    						shape = new Rectangle2D.Double(
					    						x1,
					    						y,
					    						(x2-x1),
					    						h2
					    						);
				    						}
			    						
			    						ref1+=ce.getLength();
			    						switch(op) {
			    							case H: case S: fill=Color.YELLOW;break;
			    							case X: fill=Color.RED;break;
			    							case EQ: case M: fill=(pair.isStrangelyPaired()?Color.PINK:Color.LIGHT_GRAY);break;
			    							default:break;
			    							}
			    						break;
			    				case N://through
			    				case D: shape=null;fill=null;stroke=null;ref1+=ce.getLength();break;
			    				case I: shape=null;fill=null;stroke=null;if(!hide_insertions) insertions.add(ref1);break;
			    				default: throw new IllegalStateException(""+op);
			    				}
			    			if(ref1 < this.extendedInterval.getStart()) continue;
			    			
			    			if(shape!=null && fill!=null) {
			    				canvas.shape(shape,
		    						FunctionalMap.of(
	    	    						Canvas.KEY_FILL,fill,
	    	    						Canvas.KEY_STROKE,(h2>4?stroke:null),
	    	    						Canvas.KEY_STROKE_WIDTH, 0.1
	    	    						)
		    						);
			    				}
			    			} // end loop cigar
			    		
			    		if(rec.mistmatches!=null) {
			    			for(int pos:rec.mistmatches) {
			    				final double x1 = trimX(position2pixel(pos));
			    				final double x2 = trimX(position2pixel(pos+1));
			    				if(x1>=x2) continue;
				    			canvas.rect(
					    				x1,
					    				y,
					    				(x2-x1),
					    				h2,
					    				FunctionalMap.of(
				    						Canvas.KEY_STROKE, null,
				    						Canvas.KEY_FILL,Color.ORANGE
				    						)
					    				);
				    			}
			    			}
			    		
			    		for(int px:insertions) {
			    			canvas.line(
				    				position2pixel(px),
				    				y-0.5,
				    				position2pixel(px),
				    				y+h2+0.5,
				    				FunctionalMap.of(
			    						Canvas.KEY_STROKE, Color.RED,
			    						Canvas.KEY_STROKE_WIDTH, 0.5
			    						)
				    				);
			    			}
			    		}
		    		}
		    	y-=featureHeight;
		     	}
			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));
		 	canvas.text(
					this.sampleName,
					this.boundaries.getX()+1,
					this.boundaries.getY()+fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath
						)
					); 
		 	plotFrame(canvas);
			}
		@Override
		void disposeData() {
			this.rows.clear();
			this.bamSequenceDictionary = null;
			}
		}
	/***************************************************************************************************************************/
	/***************************************************************************************************************************/
	/***************************************************************************************************************************/

	
	private class CoverageTask extends BamTask {
		double[] coverage = null;
		final Average mean_depth= new Average();
		final boolean normalize_median;
		CoverageTask(boolean normalize_median) {
			this.normalize_median=normalize_median;
		}
		
		@Override
		boolean loadData() {
			final float[] array = new float[this.extendedInterval.getLengthOnReference()];
			Arrays.fill(array, 0);
			try(SamReader sr= openSamReader()) {
				try(CloseableIterator<SAMRecord> iter= this.query(sr)) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!acceptRead(rec)) continue;
						for(AlignmentBlock ab: rec.getAlignmentBlocks()) {
							for(int x=0;x< ab.getLength();++x) {
								final int refpos1 = ab.getReferenceStart()+x;
								final int array_index = refpos1 - this.extendedInterval.getStart();
								if(array_index<0 || array_index>=array.length) continue;
								array[array_index]++;
								}
							}
						}
					}
				this.coverage = new double[boundaries.width];
				
				if(normalize_median) {
					final DiscreteMedian<Integer> discreteMedian= new DiscreteMedian<>();
					for(int i=0;i< array.length;i++) {
						this.mean_depth.accept(array[i]);
						final int refpos1 = this.extendedInterval.getStart() + i;
						if(CoordMath.overlaps(
								refpos1, refpos1,
								this.userInterval.getStart(),this.userInterval.getEnd()
								)) {
							continue;
							}
						discreteMedian.add((int)array[i]);
						}					
					final double median = Math.max(1.0,discreteMedian.getMedian().orElse(0.0));
					for(int i=0;i< array.length;i++) {
						array[i]=(float)(array[i]/median);
						}
					}
				else
					{
					for(int i=0;i< array.length;i++) {
						this.mean_depth.accept(array[i]);
						}
					}
				for(int i=0;i< this.coverage.length;i++) {
					      int x1 = (int)(((i+0)/(double) this.coverage.length)*array.length);
					final int x2 = Math.max((int)(((i+1)/(double) this.coverage.length)*array.length),x1+1);
					final Average avg=new Average();
					while(x1<x2 && x1 < array.length) {
						avg.accept(array[x1]);
						x1++;
						}
					this.coverage[i]=avg.get().orElse(0);
					}
				return true;
				}
			catch(Throwable err) {
				LOG.error(err);
				this.coverage=null;
				return false;
				}
			}
		@Override
		void disposeData() {
			this.coverage = null;
			this.bamSequenceDictionary = null;
			}
		@Override
		void plot(Canvas canvas) {
			final double bottom_y = this.boundaries.getMaxY();
			final double max_y_value = 
					this.normalize_median ?
						CoverageGrid.this.max_normalized_y:
						Math.max(1, this.mean_depth.get().orElse(1.0)*2.0)
						;
			final ToDoubleFunction<Double> cov2y = COV-> bottom_y - (COV/max_y_value)*this.boundaries.getHeight();
			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));

			canvas.comment(this.sampleName);
			plotBackgroundFrame(canvas);
			plotVariants(canvas);
			plotGenes(canvas);
			
			canvas.text(
					this.sampleName+" avgDP:"+(int)this.mean_depth.getAsDouble(),
					this.boundaries.getX()+1,
					this.boundaries.getY()+fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath
						)
					);
			// draw horizonal lines
			for(double v=0.5;v<=1.5 && this.normalize_median ;v+=0.5) {
				
				double y = cov2y.applyAsDouble(v);
				if(y< this.boundaries.getY()) continue;
				canvas.line(
						this.boundaries.getX(),
						y,
						this.boundaries.getMaxX(),
						y,
						FunctionalMap.of(
							Canvas.KEY_STROKE,(v==1.0?Color.BLUE:Color.ORANGE),
							Canvas.KEY_TITLE, String.valueOf(v),
							Canvas.KEY_STROKE_WIDTH,0.5
							)
						);
				}

			

			final List<Point2D> points = new ArrayList<>(2+this.coverage.length);
			points.add(new Point2D.Double(this.boundaries.getX(), bottom_y));
			for(int i=0;i< this.coverage.length;i++) {
				double y = Math.max(this.boundaries.getY(), cov2y.applyAsDouble(this.coverage[i]));
				points.add(new Point2D.Double(this.boundaries.getX()+i, y));
				}
			points.add(new Point2D.Double(this.boundaries.getMaxX(), bottom_y));
			int i=1;
			while(i+2 < points.size())
				{
				final Point2D p1 = points.get(i+0);
				final Point2D p2 = points.get(i+1);
				final Point2D p3 = points.get(i+2);
				if(p1.getY()==p2.getY() && p2.getY()==p3.getY() && p1.getX() < p3.getX()) {
					points.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			canvas.polygon(
					points,
					FunctionalMap.of(
							Canvas.KEY_STROKE, null,
							Canvas.KEY_FILL,Color.DARK_GRAY
							))
					;

			plotVerticalGuides(canvas);
			plotFrame(canvas);
			
			}
		}

	
	/**************************************************************************
	 * 
	 * GridTask : Co-Occurence of read names in memory
	 *
	 */
	private class GridTask extends BamTask {
		private int square_size=1;
		private int[] counts;
		
		
		private int px2pos(int m) {
			return this.extendedInterval.getStart()+(int)((m/(double)this.square_size)*this.extendedInterval.getLengthOnReference());
			}
		private double pos2px(int p) {
			return ((p-this.extendedInterval.getStart())/(double)this.extendedInterval.getLengthOnReference())*this.square_size;
			}
		
		
		private Set<String> getReadNames(final Algorithm<Interval,Integer> algorithm,final List<Interval> intervals,final Interval rgn) {
			final Set<String> names= new HashSet<>();
			final int pos = rgn.getStart();
			int index = algorithm.lower_bound(intervals, pos);
			while(index < intervals.size()) {
				final Interval r = intervals.get(index);
				if(r.getEnd() < rgn.getStart()) {
					index++;
					continue;
					}
				if(r.getStart()> rgn.getEnd()) break;
				names.add(r.getName());
				++index;
				}
			return names;
			}
		@Override
		boolean loadData() {
			try {
				/* do not use IntervalTreeMap, it's too slow */
				final List<Interval> intervals = new ArrayList<>(100_000);
				final Algorithm<Interval,Integer> algorithm = new Algorithm<>(RGN->RGN.getStart(),(A,B)->A.compareTo(B));
				try(SamReader sr= openSamReader()) {
					try(CloseableIterator<SAMRecord> iter= this.query(sr)) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							 if(!acceptRead(rec)) continue;
							 final Interval r= new Interval(rec.getContig(), rec.getUnclippedStart(), rec.getUnclippedEnd(),rec.getReadNegativeStrandFlag(),rec.getReadName());
							 if(!this.extendedInterval.overlaps(r)) continue;
							 intervals.add(r);
						     }
						}
					}
				algorithm.sortIfNeeded(intervals);
				this.square_size = Math.max(Math.min(this.boundaries.width, this.boundaries.height),1);
				this.counts=new int[this.square_size*this.square_size];
				for(int x=0;x<this.square_size;x++) {
					final Interval rx = new Interval(this.extendedInterval.getContig(),px2pos(x),px2pos(x+1));
					final Set<String> setx =  getReadNames(algorithm,intervals,rx);
					if(setx.isEmpty()) continue;
					for(int y=0;y<this.square_size;y++) {
						final Interval ry = new Interval(this.extendedInterval.getContig(),px2pos(y),px2pos(y+1));
						final Set<String> sety =  getReadNames(algorithm,intervals,ry);
						if(sety.isEmpty()) continue;
						sety.retainAll(setx);
						this.counts[y*this.square_size+x] = sety.size();
						}
					}
				return true;
				}
			catch(Throwable err) {
				LOG.error(err);
				return false;
				}
			}

		@Override
		void plot(final Canvas canvas) {
			final int margin_top = this.boundaries.y +  (this.boundaries.height - this.square_size)/2;
			final int margin_left =  this.boundaries.x +   (this.boundaries.width - this.square_size)/2;

			final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));

			canvas.comment(this.sampleName);
			
			canvas.text(
					this.sampleName,
					this.boundaries.getX()+1,
					this.boundaries.getY()+fontSize+3,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_STROKE, null,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath
						)
					);
			// area frames
			for(int side=0;side<2;++side) {
				final int p1 = side==0?this.extendedInterval.getStart():this.userInterval.getEnd();
				final int p2 = side==0?this.userInterval.getStart():this.extendedInterval.getEnd();
				final double v1 = pos2px(p1);
				final double v2 = pos2px(p2);
				if(v1>=v2) continue;
				canvas.rect(
						margin_left+v1,
						margin_top, 
						(v2-v1),
						this.square_size,
						FunctionalMap.of(
							Canvas.KEY_STROKE,null,
							Canvas.KEY_FILL,EXTERNAL_AREA_COLOR
							)
						);
				canvas.rect(
						margin_left,
						margin_top+v1, 
						this.square_size,
						(v2-v1),
						FunctionalMap.of(
							Canvas.KEY_STROKE,null,
							Canvas.KEY_FILL,EXTERNAL_AREA_COLOR
							)
						);
				}
						
			final MinMaxInteger mM = MinMaxInteger.of(Arrays.stream(this.counts));
			if(!mM.isEmpty()) {
				for(int y=0;y<this.square_size;y++) {
					int x=0;
					while(x<this.square_size) {
						if(this.counts[y*this.square_size+x]==0) {
							++x;
							continue;
							}
						final int x0=x;
						while(x<this.square_size && 
							counts[y*this.square_size+x]==counts[y*this.square_size+(x0)] ) {
							++x;
							}
						
						final float f = (float) mM.normalize(this.counts[y*this.square_size+x0]);
						final Color c =  Color.getHSBColor(f, 0.85f, 1.0f);
						canvas.rect(
								margin_left+x0,
								margin_top+(this.square_size-(1.0/square_size))-y, /* inverse y coordinate so origin is at bottom*/
								(x-x0),
								1,
								FunctionalMap.of(
									Canvas.KEY_FILL, c,
									Canvas.KEY_STROKE,null
									)
								);
						}	
					}
				}
			
			// frame
			canvas.rect(
					margin_left,
					margin_top, 
					this.square_size,
					this.square_size,
					FunctionalMap.of(
						Canvas.KEY_FILL, null,
						Canvas.KEY_STROKE,(sampleColor==null?Color.DARK_GRAY:sampleColor),
						Canvas.KEY_STROKE_WIDTH,0.1
						)
					);
			}
		@Override
		void disposeData() {
			this.counts=null;
			this.bamSequenceDictionary = null;
			}
		}
	
	/**************************************************************************
	 * 
	 * DiscordantReadsTask : plot arcs of paired reads
	 *
	 */
	private class DiscordantReadsTask extends BamTask {
		private  class Arc  {
			int p1;
			int p2;
			int readLen;
			String name;
			boolean discordant;
			public double getX1() {
				return position2pixel(p1);
				}
			public double getX2() {
				return position2pixel(p2);
				}
			boolean sameXY(final Arc arc) {
				return this.discordant==arc.discordant && Math.abs(this.getX1() - arc.getX1())<0.1 && Math.abs(this.getX2() - arc.getX2()) < 0.1;
				}
			}
		
		private final List<Arc> arcs =new ArrayList<>();
		DiscordantReadsTask() {
			}
		
		@Override
		boolean loadData() {
			
			try(SamReader sr=openSamReader()) {
				final Map<String,Arc> name2arc = new HashMap<>();
				try(CloseableIterator<SAMRecord> iter= this.query(sr)) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!rec.getReadPairedFlag()) continue;
						if(rec.getMateUnmappedFlag()) continue;
						if(!rec.getMateReferenceName().equals(this.extendedInterval.getContig())) continue;
						if(!acceptRead(rec)) continue;
						final Arc arc = new Arc();
						arc.name = rec.getReadName();
						arc.p1 = Math.min(rec.getUnclippedStart(), rec.getMateAlignmentStart());
						arc.p2 = Math.max(rec.getUnclippedEnd(), getMateEnd(rec));
						/* use length on reference rather than real read length */
						arc.readLen = Math.max(
								CoordMath.getLength(rec.getUnclippedStart(),rec.getUnclippedEnd()),
								CoordMath.getLength(rec.getMateAlignmentStart(),  getMateEnd(rec))
								);
						arc.discordant = !rec.getProperPairFlag() || rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag();
						
						final Arc other= name2arc.getOrDefault(arc.name,null);
						if(other!=null) {
							other.p1 = Math.min(other.p1,arc.p1);
							other.p2 = Math.max(other.p2,arc.p2);
							other.readLen = Math.max(other.readLen, arc.readLen);
							}
						else
							{
							name2arc.put(arc.name,arc);
							}
						}
					this.arcs.addAll(name2arc.values());
					}
				return true;
				}
			catch(Throwable err) {
				LOG.error(err);
				return false;
				}
			}
		protected double getGeneMidY() {
			return this.boundaries.getCenterY();
			}
		@Override
		void plot(final Canvas canvas) {
			 canvas.comment(this.sampleName);
		     final double midy=getGeneMidY();
		     plotBackgroundFrame(canvas);
		     plotVariants(canvas);
		     plotVerticalGuides(canvas);
		     
		     //pileup Arcs so arcs at the same X1-X2 are displayed with a small shift
		     final List<List<Arc>> pileup = new ArrayList<>();
		     while(!arcs.isEmpty()) {
		    	final Arc arc = arcs.remove(arcs.size()-1);
				if(!CoordMath.encloses(
					this.extendedInterval.getStart(), this.extendedInterval.getEnd(),
					arc.p1, arc.p2)) continue;
				
				
				if(!(
						CoordMath.overlaps(arc.p1,arc.p2, this.userInterval.getStart(), this.userInterval.getStart()) ||
						CoordMath.overlaps(arc.p1,arc.p2, this.userInterval.getEnd(), this.userInterval.getEnd())
						)) continue;
		    	 
				
				if(CoverageGrid.this.inversion_flag) {
					// ARC ALL OUTSIDE INV
			    	if(arc.p1 + arc.readLen < this.userInterval.getStart() && 
			    			this.userInterval.getEnd()+ arc.readLen < arc.p2) {
			    		continue;
			    		}
			    	// ARC ALL INSIDE INV
			    	if(this.userInterval.getStart() + arc.readLen < arc.p1 && 
			    		arc.p2 + arc.readLen < this.userInterval.getEnd()) {
			    		continue;
			    		}
		    		}
				
		    	 int y=0;
		    	 for(y=0;y<pileup.size();++y) {
		    		final List<Arc> row =pileup.get(y);
		    		if(row.get(0).sameXY(arc)) {
		    			row.add(arc);
		    			break;
		    			}
		    	 	}
		    	 if(y==pileup.size()) {
		    		 final List<Arc> row= new ArrayList<>();
		    		 row.add(arc);
		    		 pileup.add(row);
		    	 	}
		     }
		     
		     
		     
		     final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));
		     final Hyperlink hyperlink = Hyperlink.compile(this.bamSequenceDictionary);
		     canvas.text(
						this.sampleName+" "+
							pileup.stream().flatMap(T->T.stream()).filter(A->A.discordant).count()+
							" / "+
							pileup.stream().flatMap(T->T.stream()).filter(A->!A.discordant).count(),
						this.boundaries.getX()+1,
						this.boundaries.getY()+fontSize+3,
						FunctionalMap.of(
							Canvas.KEY_FONT_SIZE, fontSize,
							Canvas.KEY_STROKE, null,
							Canvas.KEY_FILL, Color.DARK_GRAY,
							Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath,
							Canvas.KEY_HREF,hyperlink.apply(this.extendedInterval).orElse(null)
							)
						);
		     
		     // horizontal line
		     canvas.line(
		    		 this.boundaries.getX(),
		    		 midy,
		    		 this.boundaries.getMaxX(),
		    		 midy,
	    			FunctionalMap.of(
    					Canvas.KEY_STROKE,Color.DARK_GRAY,
    					Canvas.KEY_STROKE_WIDTH,0.1
    					)
	    			);
		     
		     plotGenes(canvas);
		     
		   
			for(List<Arc> row: pileup) {
				final Arc first = row.get(0);
		    	final double len = first.getX2()-first.getX1();
		    	final double max_height = Math.min(((this.boundaries.getHeight()*0.9)/1.0),len) ;
		    	final double half_height = max_height/2;
				
			     for(int i=0;i< row.size();++i) {
			    	final Arc arc= row.get(i);
			    	
			    	// INV is outside user boundaries
			    	final double height = max_height - ((max_height-half_height)/(double)row.size())*i;
			    			    	
			    	final Arc2D path =new Arc2D.Double(
			    			new Rectangle2D.Double(
				    			arc.getX1(),
				    			midy - height/2,
				    			(arc.getX2()-arc.getX1()),
				    			height
				    			),
			    			(arc.discordant?0:180),
			    			180,
			    			Arc2D.OPEN
			    			);
			    	canvas.shape(path,
		    			FunctionalMap.of(
	    					Canvas.KEY_FILL, null,
	    					Canvas.KEY_STROKE,	Color.getHSBColor( (arc.discordant?0.f:0.3f)+0.1f*(i/(float)row.size()), 0.5f, 0.75f),
	    					Canvas.KEY_STROKE_WIDTH,0.1
	    					)
		    			);
			     	}
			     }
			plotFrame(canvas);
			}
		@Override
		void disposeData() {
			arcs.clear();
			this.bamSequenceDictionary = null;
			}
	}
	
	

	/**************************************************************************
	 * 
	 * SupplementaryReadTask : plot arcs of supplementary reads
	 *
	 */
	private class SupplementaryReadTask extends BamTask {
		private  class ClipArc extends AbstractLocatable {
			final SAMRecord R1;
			final SAMRecord R2;
			ClipArc(SAMRecord rec,SAMRecord suppl) {
				if(rec.getStart()< suppl.getStart()) {
					this.R1=rec;
					this.R2=suppl;
					}
				else
					{
					this.R1=suppl;
					this.R2=rec;
					}
				}
			@Override
			public String getContig() { return R1.getContig();}
			@Override
			public int getStart() {
				return Math.min(R1.getStart(), R2.getStart());
			}
			@Override
			public int getEnd() {
				return Math.max(R1.getEnd(), R2.getEnd());
			}
			boolean sameXY(final ClipArc arc) {
				return Math.abs(this.getX1() - arc.getX1())<0.1 && Math.abs(this.getX2() - arc.getX2()) < 0.1;
				}
			public int getUnclippedStart() { return  Math.min(R1.getUnclippedStart(), R2.getUnclippedStart());}
			public int getUnclippedEnd() { return  Math.max(R1.getUnclippedEnd(), R2.getUnclippedEnd());}
			public double getX1() {
				return position2pixel(getStart());
				}
			public double getX2() {
				return position2pixel(getEnd()+1);
				}
			}
		
		private final List<ClipArc> arcs =new ArrayList<>();
	
		
		@Override
		boolean loadData() {
			try(SamReader sr=openSamReader()) {
				try(CloseableIterator<SAMRecord> iter= this.query(sr)) {
					while(iter.hasNext()) {
						final SAMRecord rec = iter.next();
						if(!acceptRead(rec)) continue;
						for(SAMRecord suppl:SAMUtils.getOtherCanonicalAlignments(rec)) {
							if(!suppl.getContig().equals(this.extendedInterval.getContig())) continue;
							final ClipArc arc = new ClipArc(rec,suppl);
							
							if(!CoordMath.encloses(
									this.extendedInterval.getStart(), this.extendedInterval.getEnd(),
									arc.getStart(), arc.getEnd())) continue;
							
							if(!(
									CoordMath.overlaps(arc.getUnclippedStart(),arc.getUnclippedEnd(), this.userInterval.getStart(), this.userInterval.getStart()) ||
									CoordMath.overlaps(arc.getUnclippedStart(),arc.getUnclippedEnd(), this.userInterval.getEnd(), this.userInterval.getEnd())
									)) continue;
							
							if(CoverageGrid.this.inversion_flag) {
								// ARC ALL OUTSIDE INV
						    	if(arc.R1.getUnclippedEnd() < this.userInterval.getStart() && 
						    			this.userInterval.getEnd()< arc.R2.getUnclippedStart()) {
						    		continue;
						    		}
						    	// ARC ALL INSIDE INV
						    	if(this.userInterval.getStart() < arc.getUnclippedStart() && 
						    		arc.getUnclippedEnd() < this.userInterval.getEnd()) {
						    		continue;
						    		}
					    		}
							
							this.arcs.add(arc);
							}
						}
					}
				return true;
				}
			catch(Throwable err) {
				LOG.error(err);
				return false;
				}
			}
		protected double getGeneMidY() {
			return this.boundaries.getCenterY();
			}
		@Override
		void plot(final Canvas canvas) {
			 canvas.comment(this.sampleName);
		     final double midy=getGeneMidY();
		     plotBackgroundFrame(canvas);
		     plotVariants(canvas);
		     plotVerticalGuides(canvas);
		     
		     //pileup Arcs so arcs at the same X1-X2 are displayed with a small shift
		     final List<List<ClipArc>> pileup = new ArrayList<>();
		     while(!arcs.isEmpty()) {
		    	final ClipArc arc = arcs.remove(arcs.size()-1);
				
		    	 int y=0;
		    	 for(y=0;y<pileup.size();++y) {
		    		final List<ClipArc> row =pileup.get(y);
		    		if(row.get(0).sameXY(arc)) {
		    			row.add(arc);
		    			break;
		    			}
		    	 	}
		    	 if(y==pileup.size()) {
		    		 final List<ClipArc> row= new ArrayList<>();
		    		 row.add(arc);
		    		 pileup.add(row);
		    	 	}
		     }
		     
		     
		     
		     final int fontSize=Math.min(12,(int)(this.boundaries.getHeight()/10.0));
		     final Hyperlink hyperlink = Hyperlink.compile(this.bamSequenceDictionary);
		     canvas.text(
						this.sampleName+" ("+
							pileup.stream().flatMap(T->T.stream()).count()+
							")"
							,
						this.boundaries.getX()+1,
						this.boundaries.getY()+fontSize+3,
						FunctionalMap.of(
							Canvas.KEY_FONT_SIZE, fontSize,
							Canvas.KEY_STROKE, null,
							Canvas.KEY_FILL, Color.DARK_GRAY,
							Canvas.KEY_TITLE,this.sampleName+" "+this.bamPath,
							Canvas.KEY_HREF,hyperlink.apply(this.extendedInterval).orElse(null)
							)
						);
		     
		     // horizontal line
		     canvas.line(
		    		 this.boundaries.getX(),
		    		 midy,
		    		 this.boundaries.getMaxX(),
		    		 midy,
	    			FunctionalMap.of(
    					Canvas.KEY_STROKE,Color.DARK_GRAY,
    					Canvas.KEY_STROKE_WIDTH,0.1
    					)
	    			);
		     
		     plotGenes(canvas);
		     
		     
			for(List<ClipArc> row : pileup) {
				final ClipArc first = row.get(0);
		    	final double len = first.getX2()-first.getX1();
		    	final double max_height = Math.min(((this.boundaries.getHeight()*0.9)/1.0),len) ;
		    	final double half_height = max_height/2;
				
			     for(int i=0;i< row.size();++i) {
			    	final ClipArc arc= row.get(i);
			    	
			    	// INV is outside user boundaries
			    	final double height = max_height - ((max_height-half_height)/(double)row.size())*i;
			    	final boolean same_strand = arc.R1.getReadNegativeStrandFlag()!=arc.R2.getReadNegativeStrandFlag();
			    	final Arc2D path =new Arc2D.Double(
			    			new Rectangle2D.Double(
				    			arc.getX1(),
				    			midy - height/2,
				    			(arc.getX2()-arc.getX1()),
				    			height
				    			),
			    			(same_strand ?0:180),
			    			180,
			    			Arc2D.OPEN
			    			);
			    	canvas.shape(path,
		    			FunctionalMap.of(
	    					Canvas.KEY_FILL, null,
	    					Canvas.KEY_STROKE,	Color.getHSBColor((same_strand?0.f:0.4f)+0.1f*(i/(float)row.size()), 0.5f, 0.75f),
	    					Canvas.KEY_STROKE_WIDTH,0.1
	    					)
		    			);
			     	}
			     }
			plotFrame(canvas);
			}
		@Override
		void disposeData() {
			arcs.clear();
			this.bamSequenceDictionary = null;
			}
	}	
	
/********************************************************************************************************************/
	
/** get mate end position or mate-start */
private static int getMateEnd(final SAMRecord rec0) {
	return rec0.hasAttribute(SAMTag.MC)?SAMUtils.getMateAlignmentEnd(rec0):rec0.getMateAlignmentStart();
	}

private Set<String> getGeneNames(final Locatable loc) {
	if(this.gtfFile==null) return Collections.emptySet();
	final GTFCodec codec=new GTFCodec();
	try(TabixReader tabix = new TabixReader(this.gtfFile)) {
		final Set<String> set = new TreeSet<>();
		TabixReader.Iterator iter= tabix.query(loc.getContig(), loc.getStart(), loc.getEnd());
		for(;;) {
			final String line = iter.next();
			if(line==null) break;
			final GTFLine feat=codec.decode(line);
			if(feat==null) continue;
			if(feat.hasAttribute("gene_name")) {
				set.add(feat.getAttribute("gene_name"));
				continue;
				}
			if(feat.hasAttribute("gene_id")) {
				set.add(feat.getAttribute("gene_id"));
				continue;
				}
			}
		return set;
		}
	catch(Exception err) {
		return Collections.emptySet();
		}
	}

@Override
public int doWork(final List<String> args) {
	try
		{
		if(this.extendFactor < 1.5) {
			System.err.println("extendFactor should be >=1.5");
			return -1;
			}
		if(this.max_normalized_y < 2.0) {
			System.err.println("max Y should be >= 2.0");
			return -1;
			}
		
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		
		
		final Path input =  Paths.get(this.oneAndOnlyOneFile(args));
		final SampleSheet sampleSheet = new SampleSheetFactory().splitter(CharSplitter.TAB).of(input);
		sampleSheet.getHeader().assertColumnExists("bam");
		
		if(sampleSheet.isEmpty()) {
			LOG.error("empty samplesheet "+input);
			return -1;
			}
		
		final SimpleInterval user_interval = IntervalParserFactory.newInstance(dict).
				make().
				apply(this.intervalStr).
				orElseThrow(()->new IllegalArgumentException("Cannot parse interval:\""+this.intervalStr+"\""));
		
		final int mid_position = user_interval.getStart() + user_interval.getLengthOnReference()/2;
		final int new_len = (int)Math.max(
				user_interval.getLengthOnReference()+2,
				Math.ceil(user_interval.getLengthOnReference()* this.extendFactor)
				);
		final int new_start = Math.max(1, mid_position-new_len/2);
		final int new_end = Math.min(dict.getSequence(user_interval.getContig()).getEnd(), mid_position+new_len/2);
		final SimpleInterval extended_interval = new SimpleInterval(user_interval.getContig(),new_start,new_end);

		final List<BamTask> tasks = new ArrayList<>(sampleSheet.size());
		
		final int marginTop=50;

		final int ncols = (int)Math.max(1,Math.floor(Math.sqrt(sampleSheet.size())));
		final int nrows = (int)Math.max(1, Math.ceil(sampleSheet.size()/(double)ncols));
		final double wi  = this.dimension.getWidth()/ncols;
		final double hi  = (this.dimension.getHeight()-marginTop)/nrows;
		
		LOG.info("tile:"+ncols+"x"+nrows+" w:"+(int)wi+" h:"+(int)hi);
		
		if(wi < 10 || hi < 10) {
			LOG.error("dimension is too small");
			return -1;
			}
		
		final byte[] bases;
		int nr=0;
		int nc=0;
		try(ReferenceSequenceFile refseq = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.refPath) ) {
			bases = refseq.getSubsequenceAt(
					extended_interval.getContig(),
					extended_interval.getStart(),
					extended_interval.getEnd()
					).getBases();
			}
		
		for(final FileHeader.RowMap rowMap: sampleSheet) {
			final BamTask task ;
				switch(this.plotType) {
				case GRID: task= new GridTask(); break;
				case PILEUP: task= new PileupTask(false); break;
				case PILEUP_PAIR: task= new PileupTask(true); break;
				case COVERAGE:  task = new CoverageTask(false); break;
				case DISCORDANT: task = new DiscordantReadsTask();break;
				case MEDIAN_COVERAGE: task= new CoverageTask(true);break;
				case SUPPL: task = new SupplementaryReadTask();break;
				default: {
					LOG.error("unedefine "+this.plotType);
					return -1;
					}
				}
			task.reference_bases = bases;
			task.bamPath = Paths.get(rowMap.get("bam"));
			task.sampleName = rowMap.getOrDefault("sample","");
			task.sampleColor  = rowMap.containsKey("color")?new ColorUtils().parse(rowMap.get("color")):null;
			task.vcfPath = this.vcfFile;
			task.gtfPath = this.gtfFile;
			task.referencePath = this.refPath;
			task.userInterval = user_interval;
			task.extendedInterval = extended_interval;
			task.boundaries = new Rectangle(
					(int)(nc*wi), 
					(int)(marginTop+nr*hi),
					(int)wi-1,
					(int)hi-1
					);
			tasks.add(task);
			nc++;
			if(nc==ncols) {
				nc=0;
				nr++;
				}
			}
	   
		
		final String title= SequenceDictionaryUtils.getBuildName(dict).orElse("") + " "+ 
				user_interval.toNiceString()+
			" length:"+ StringUtils.niceInt(user_interval.getLengthOnReference())+
			" extended to "+extended_interval.toNiceString()+
			" length:"+StringUtils.niceInt(extended_interval.getLengthOnReference())+
			" model:"+this.plotType.name()+" "+
			this.title
			;
		try(Canvas canvas = new CanvasFactory().
				setWidth(this.dimension.width).
				setHeight(this.dimension.height).
				setFormat(this.canvas_format).
				open(this.outputFile,FunctionalMap.of(
					Canvas.KEY_TITLE,title,
					Canvas.KEY_FILL, Color.BLACK,
					Canvas.KEY_STROKE,null
					) )
				)
			{
			final int fontSize=Math.min(14,marginTop);
			final String url =Hyperlink.compile(dict).apply(extended_interval).orElse(null);

			for(BamTask task: tasks) {
				task.canvas = canvas;
				}
			
			canvas.text(
					title,
					10,
					fontSize,
					FunctionalMap.of(
						Canvas.KEY_FONT_SIZE, fontSize,
						Canvas.KEY_FILL, Color.DARK_GRAY,
						Canvas.KEY_HREF, url,
						Canvas.KEY_TITLE, title
						)
					);

			final Set<String> genes = getGeneNames(user_interval);
			if(!genes.isEmpty()) {
				canvas.text(
						(genes.size()>10?
								StringUtils.niceInt(genes.size())+" genes":
									String.join(" ",genes)),
						10,
						fontSize*2,
						FunctionalMap.of(
							Canvas.KEY_FONT_SIZE, fontSize,
							Canvas.KEY_FILL, Color.GRAY
							)
						);
				}
			
			 final ExecutorService executorService = Executors.newFixedThreadPool(this.nThreads);
				for(int i=0; i < tasks.size();i+=this.nThreads) {
					final List<BamTask> batch = tasks.subList(i, Math.min(i + this.nThreads, tasks.size()));
					batch.forEach(executorService::execute);
					}
				executorService.shutdown();
				executorService.awaitTermination(365, TimeUnit.DAYS);
				}
			
	
		return 0;
		}
	catch(final Throwable err)
		{
		err.printStackTrace();
		LOG.error(err);
		return -1;
		}

	}

public static void main(final String[] args) {
	new CoverageGrid().instanceMainWithExit(args);
	}

}
