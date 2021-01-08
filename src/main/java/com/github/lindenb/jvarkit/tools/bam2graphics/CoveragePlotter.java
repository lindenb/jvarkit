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
package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.IntToDoubleFunction;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.imageio.ImageIO;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.swing.GraphicsState;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC
input is an interval of a file source of interval (bed, vcf, gtf, interval_list , ,etc...)


```
java -jar dist/coverageplotter.jar -R src/test/resources/rotavirus_rf.fa -B src/test/resources/S1.bam -B src/test/resources/S2.bam "RF01:1-4000" -w 50 | less -r
```


END_DOC 
 */
@Program(
	name="coverageplotter",
	description="Find anomaly of depth in intervals+bams",
	keywords={"cnv","bam","depth","coverage"},
	creationDate="20200605",
	modificationDate="20200618"
	)
public class CoveragePlotter extends Launcher {
	private static final Logger LOG = Logger.build( CoveragePlotter.class).make();
	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"-B","--bams"},description = "list of bams. one file with a '.list' suffix is interpretted a a list of path to the bams",required=true)
	private List<String> bamsPath= new ArrayList<>();
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--max-depth"},description = "ignore position if depth > 'x'")
	private int max_depth=500;
	@Parameter(names={"--dimension"},description = "Image Dimension. " + DimensionConverter.OPT_DESC, converter=DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--extend","-x"},description = "extend original interval by this fraction")
	private double extend=1.0;
	@Parameter(names= {"--gtf"},description="Optional Tabix indexed GTF file. Will be used to retrieve an interval by gene name, or to display gene names in a region.")
	private Path gtfFile = null;
	@Parameter(names= {"--known"},description="Optional Tabix indexed Bed or VCF file containing known CNV. Both types must be indexed.")
	private Path knownCnvFile = null;
	@Parameter(names= {"--prefix"},description="Image File Prefix.")
	private String prefix="";
	@Parameter(names= {"--alpha"},description="line opacity. "+ FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double alpha=1.0;
	@Parameter(names= {"--arc-alpha"},description="arc opacity. "+ FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double alpha_arc=1.0;
	@Parameter(names= {"--rrff"},description="Only display arcs where the strands the the read and its mate are Forward-Forward or Reverse-Reverse")
	private boolean arc_only_rr_ff;
	@Parameter(names= {"--manifest"},description="Optional. Manifest file")
	private Path manifestPath =null;
	@Parameter(names= {"--min-arc"},description="min arc length in bp.",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_invert = 1000;
	@Parameter(names= {"--max-arc"},description="max arc length in bp.",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_invert = 10_000_000;
	@DynamicParameter(names = "-D", description = "style. Undocumented.",hidden=true)
	private Map<String, String> dynaParams = new HashMap<>();
	@Parameter(names = {"--black","--exclude"}, description = "Optional. BED Tabix indexed black-listed region")
	private Path blackListedPath=null;
	@Parameter(names= {"--smooth"},description="sliding window smooth size.",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int window_smooth_size = 250;
	@Parameter(names= {"--skip-center"},description="When calculating the median depth, only consider the extended region, not the original interval.")
	private boolean skip_original_interval_for_median = false;
	@Parameter(names= {"--ignore-known-containing"},description="Ignore known CNV containing the whole region (prevent large known CNV to be displayed) ")
	private boolean ignore_cnv_overlapping = false;
	
	private final ColorUtils colorUtils = new ColorUtils();
	
	private Color getColor(final String key,final Color defColor) {
		if(!this.dynaParams.containsKey(key)) {
			return defColor;
			}
		Color c = colorUtils.parse(this.dynaParams.get(key));
		return c==null?defColor:c;
		}


	private void drawKnownCnv(final Graphics2D g,final Rectangle rectangle,final Locatable region) {
		if(this.knownCnvFile==null) return;
		final String fname=this.knownCnvFile.getFileName().toString();

		final Pileup<Interval> pileup = new Pileup<>();
		final Predicate<Interval> rejectCnv = cnv->(this.ignore_cnv_overlapping && cnv.getStart() < region.getStart() && cnv.getEnd() > region.getEnd());

		if(fname.endsWith(".bed.gz")) {
			try(TabixReader tbr = new TabixReader(this.knownCnvFile.toString())) {
				final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				final String ctg = cvt.apply(region.getContig());
					if(!StringUtils.isBlank(ctg)) {
						final BedLineCodec codec = new BedLineCodec();
						final TabixReader.Iterator iter = tbr.query(ctg,region.getStart(), region.getEnd());
						for(;;) {
							final String line = iter.next();
							if(line==null) break;
							final BedLine bed = codec.decode(line);
							if(bed==null) continue;
							final Interval rgn = new Interval(region.getContig(),bed.getStart(),bed.getEnd(),false,bed.getOrDefault(3, ""));
							if(rejectCnv.test(rgn)) return;
							pileup.add(rgn);
							}
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			}
		else if(FileExtensions.VCF_LIST.stream().anyMatch(X->fname.endsWith(X))) {
			try(VCFReader vcfFileReader= VCFReaderFactory.makeDefault().open(this.knownCnvFile,true)) {
				final ContigNameConverter cvt = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(vcfFileReader.getHeader()));
				final String ctg = cvt.apply(region.getContig());
				if(!StringUtils.isBlank(ctg)) {
					vcfFileReader.query(ctg, region.getStart(), region.getEnd()).
							stream().
							filter(VC->!VC.isSNP()).
							forEach(VC->{
								final List<String> list = new ArrayList<>();
								if(VC.hasID()) list.add(VC.getID());
								if(VC.hasAttribute(VCFConstants.SVTYPE))  list.add(VC.getAttributeAsString(VCFConstants.SVTYPE,"."));
								final Interval rgn= new Interval(region.getContig(),VC.getStart(),VC.getEnd(),false,String.join(";",list));
								if(rejectCnv.test(rgn)) return;
								pileup.add(rgn);
								});
					}
				}
			catch(final Throwable err) {
				LOG.error(err);
				}
			}
		else
			{
			LOG.warn("not a vcf of bed.gz file "+this.knownCnvFile);
			}
		if(!pileup.isEmpty() ) {
			final Composite oldComposite = g.getComposite();
			final Stroke oldStroke = g.getStroke();
			g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.4f));
			g.setColor(getColor("known",Color.MAGENTA));

			final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rectangle.getWidth();
			final double featureHeight = 4.0/pileup.getRowCount();
			for(int row=0;row< pileup.getRowCount();++row) {
				for(final Interval cnv:pileup.getRow(row)) {
					final double y= rectangle.getHeight()-8.0 + row*featureHeight;
					final double x1 = position2pixel.applyAsDouble(cnv.getStart());
					final double x2 = position2pixel.applyAsDouble(cnv.getEnd());
					g.draw(new Rectangle2D.Double(x1, y-1, Math.max(1.0,x2-x1),featureHeight*0.9));
					}
				}
			g.setComposite(oldComposite);
			g.setStroke(oldStroke);
			}
		}

	

	private Stream<GTFLine> getGenes(final Locatable region) {
		if(this.gtfFile==null) return Stream.empty();
		TabixReader tbr = null;
		try {
			tbr= new TabixReader(this.gtfFile.toString());
			
			final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
			final String ctg = cvt.apply(region.getContig());
			if(StringUtils.isBlank(ctg)) {
				tbr.close();
				return Stream.empty();
				}
			
			final GTFCodec codec = new GTFCodec();
			final TabixReader.Iterator iter=tbr.query(ctg,region.getStart(),region.getEnd());
			final TabixReader tbrfinal = tbr;
			final AbstractIterator<GTFLine> iter2= new AbstractIterator<GTFLine>() {
				@Override
				protected GTFLine advance() {
					try {
						for(;;) {
							final String line = iter.next();
							if(line==null) return null;
							if(StringUtils.isBlank(line) ||  line.startsWith("#")) continue;
							final String tokens[]= CharSplitter.TAB.split(line);
							if(tokens.length<9 ) continue;
							tokens[0]=region.getContig();
							final GTFLine gtfline = codec.decode(line);
							if(gtfline==null) continue;
							return gtfline;
							}
						} catch (final IOException e) {
						LOG.error(e);
						return null;
						}
					}
				};
			return StreamSupport.stream(new IterableAdapter<GTFLine>(iter2).spliterator(),false).
					onClose(()->{ tbrfinal.close(); });
			}
		catch(Throwable err) {
			if(tbr!=null) tbr.close();
			return Stream.empty();
			}
	}
	
	
	private void drawGenes(final Graphics2D g,final Rectangle rect,final Locatable region) {
		final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rect.getWidth();
		final double y= rect.getMaxY()-4.0;
		try( GraphicsState state = GraphicsState.of(g)) {
			final int geneSize = 10;
			state.getGraphics().setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.5f));
			g.setFont(new Font(g.getFont().getName(), Font.PLAIN, geneSize));
			g.setColor(getColor("gene",Color.DARK_GRAY));
	
			getGenes(region).forEach(gtfline->{
					final double x1 = position2pixel.applyAsDouble(gtfline.getStart());
					final double x2 = position2pixel.applyAsDouble(gtfline.getEnd());
	
					if(gtfline.getType().equals("gene") ) {
						g.drawString(gtfline.getAttribute("gene_name"),(int)Math.max(x1,1),(int)(y-(geneSize+3)));
						}
					else if(gtfline.getType().equals("exon") ) {
						g.draw(new Rectangle2D.Double(x1, y-1, (x2-x1), 3));
						}
					else if(gtfline.getType().equals("transcript") ) {
						g.draw(new Line2D.Double(x1, y, x2, y));
						}
					});
			}
		}
	
	
	
	

@Override
public int doWork(final List<String> args) {
	ArchiveFactory archive = null;
	PrintWriter manifest = null;
	try
		{
	
		
		if(extend<1.0) {
			LOG.error("extend is lower than 1 :"+this.extend);
			return -1;
			}
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					referenceSequence(CoveragePlotter.this.refPath).
					validationStringency(ValidationStringency.LENIENT)
					;
		
		 final List<Path> inputBams =  IOUtils.unrollPaths(this.bamsPath);
		
		if(inputBams.isEmpty()) {
			LOG.error("input bam file missing.");
			return -1;
			}
		
		 Iterator<? extends Locatable> iter;
		 final String input = oneFileOrNull(args); 
		 if(input==null) {
			 final BedLineCodec codec = new BedLineCodec();
			 final LineIterator liter = new LineIterator(stdin());
			 iter = new AbstractIterator<Locatable>() {
			 	@Override
			 	protected Locatable advance() {
			 		while(liter.hasNext()) {
			 			final String line = liter.next();
			 			final BedLine bed = codec.decode(line);
			 			if(bed==null) {
			 				continue;
			 				}
			 			return bed;
			 			}
			 		liter.close();
			 		return null;
			 		}
			 	};
		 	}
		 else
		 	{
			iter = IntervalListProvider.from(input).dictionary(dict).stream().iterator();
		 	}
		final BufferedImage image = new BufferedImage(this.dimension.width,this.dimension.height,BufferedImage.TYPE_INT_ARGB);
		final BufferedImage offscreen = new BufferedImage(this.dimension.width,this.dimension.height,BufferedImage.TYPE_INT_ARGB);
		final double y_mid = this.dimension.getHeight()/2.0;
		final ToDoubleFunction<Double> normToPixelY = NORM->  this.dimension.getHeight() - NORM*y_mid;
		final DiscreteMedian<Integer> discreteMedian = new DiscreteMedian<>();
		
		manifest = (this.manifestPath==null?new PrintWriter(new NullOuputStream()):IOUtils.openPathForPrintWriter(this.manifestPath));
		archive = ArchiveFactory.open(this.outputFile);
		while(iter.hasNext()) {
			final Locatable the_locatable = iter.next();
			String label="";
			if(the_locatable instanceof BedLine) {
				final BedLine bedline = BedLine.class.cast(the_locatable);
				label= bedline.getOrDefault(3,label);
				}
			
			final SimpleInterval rawRegion = new SimpleInterval(the_locatable);
			final SimpleInterval extendedRegion;
			/* extend interval */
			if(this.extend>1) {
				final int L1 = rawRegion.getLengthOnReference();
				final  int L2 = (int)Math.ceil(L1*this.extend);
				final int mid = rawRegion.getCenter().getPosition();
				int x0 = mid - L2/2;
				if(x0<0) x0=1;
				int x1= mid + L2/2;
				final SAMSequenceRecord ssr = dict.getSequence(rawRegion.getContig());
				if(ssr!=null) x1=Math.min(ssr.getSequenceLength(), x1);
				if(x0>x1) continue;
				extendedRegion = new SimpleInterval(rawRegion.getContig(),x0,x1);
				}
			else
				{
				extendedRegion=rawRegion;
				}
			final ToDoubleFunction<Integer> pos2pixel = POS-> (POS - extendedRegion.getStart())/(double)extendedRegion.getLengthOnReference() * this.dimension.getWidth();

			
			final String theGeneName = this.getGenes(extendedRegion).
					filter(G->G.getType().equals("gene")).
					map(G->G.getAttribute("gene_name")).
					filter(S->!StringUtils.isBlank(S)).findFirst().orElse(null);
			
			final String outputFilename = this.prefix + extendedRegion.getContig()+"_"+extendedRegion.getStart()+"_"+extendedRegion.getEnd()+
					(StringUtils.isBlank(label)?"":"."+label.replaceAll("[^A-Za-z\\-\\.0-9]+", "_")) + 
					(StringUtils.isBlank(theGeneName) || !StringUtils.isBlank(label)?"":"."+theGeneName.replaceAll("[^A-Za-z\\-\\.0-9]+", "_")) + 
					".png";
			
			
			final Graphics2D g2= offscreen.createGraphics();
			g2.setColor(Color.BLACK);
			g2.setComposite(AlphaComposite.getInstance(AlphaComposite.CLEAR));
			g2.fillRect(0, 0, this.dimension.width,this.dimension.height);
			g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER));
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);

			
			final Graphics2D g= image.createGraphics();
			g.setColor(Color.WHITE);
			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
			g.fillRect(0, 0, this.dimension.width,this.dimension.height);
			//draw exons
			getGenes(extendedRegion).filter(G->G.getType().equals("exon")).forEach(EX->{
				g.setColor(new Color(240,240,240));
				double x1 = pos2pixel.applyAsDouble(EX.getStart());
				double x2 = pos2pixel.applyAsDouble(EX.getEnd());
				g.fill(new Rectangle2D.Double(x1,0, (x2-x1),this.dimension.height));
			});
			
			
			int y =(int)(this.dimension.height/2.0);
			g.setColor(Color.BLUE);
			g.drawLine(0, y, image.getWidth(),y );
			y=(int)(this.dimension.height/4.0);
			g.setColor(Color.CYAN);
			g.drawLine(0, y,image.getWidth(),y );
			y=(int)(3.0*this.dimension.height/4.0);
			g.drawLine(0, y,image.getWidth(),y );
			g.setColor(Color.DARK_GRAY);
			g.drawRect(0, 0, this.dimension.width-1,this.dimension.height-1);
			drawGenes(g, new Rectangle(0,0,image.getWidth(),image.getHeight()), extendedRegion);
			drawKnownCnv(g, new Rectangle(0,0,image.getWidth(),image.getHeight()), extendedRegion);
			if(this.extend>1) {
				g.setColor(Color.GREEN);
				int x =(int) pos2pixel.applyAsDouble(rawRegion.getStart());
				g.drawLine(x, 0, x, image.getHeight() );
				x = (int) pos2pixel.applyAsDouble(rawRegion.getEnd());
				g.drawLine(x, 0, x, image.getHeight() );
				}
			
			final int depth[]= new int[extendedRegion.getLengthOnReference()];
			final int copy[]= new int[depth.length];
			final BitSet blackListedPositions = new BitSet(depth.length);
			final Map<String,Point2D> sample2maxPoint = new HashMap<>(inputBams.size());
			boolean drawAbove = false;
			
			// fill black listed regions
			if(this.blackListedPath!=null) {
				try(TabixReader tbr= new TabixReader(this.blackListedPath.toString())) {
					final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
					final String ctg = cvt.apply(extendedRegion.getContig());
					if(!StringUtils.isBlank(ctg)) {
						final BedLineCodec codec = new BedLineCodec();
						final TabixReader.Iterator tbxr = tbr.query(ctg,extendedRegion.getStart(), extendedRegion.getEnd());
						for(;;) {
							final String line = tbxr.next();
							if(line==null) break;
							final BedLine bed = codec.decode(line);
							if(bed==null) continue;
							int p1 = Math.max(bed.getStart(), extendedRegion.getStart());
							while(p1 <= extendedRegion.getEnd()  && p1 <= bed.getEnd()) {
								blackListedPositions.set(p1-extendedRegion.getStart());
								++p1;
								}
							}
						}
					}
				catch(Throwable err) {
					LOG.warn(err);
					}
				}
			
			if(this.skip_original_interval_for_median && this.extend>1) {
				for(int i=the_locatable.getStart();i<=the_locatable.getEnd();i++) {
					blackListedPositions.set(i-extendedRegion.getStart());
					}
				}
			
			
			final Set<String> sampleNames = new TreeSet<>();
			for(final Path path: inputBams) {
				try(SamReader sr = samReaderFactory.open(path)) {
					final SAMFileHeader header= sr.getFileHeader();
					
					final String sample = header.getReadGroups().stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
					
					sampleNames.add(sample);
					/* sample has deletion */
					int count_has_dp_le_0_5 = 0;
					/* sample has dup */
					int count_has_dp_ge_1_5 = 0;
					/* longest arc */
					int longest_arc = 0;
					
					SequenceUtil.assertSequenceDictionariesEqual(dict,header.getSequenceDictionary());
					Arrays.fill(depth, 0);
					try(CloseableIterator<SAMRecord> siter = sr.queryOverlapping(extendedRegion.getContig(), extendedRegion.getStart(), extendedRegion.getEnd())) {
						while(siter.hasNext()) {
							final SAMRecord rec= siter.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(!SAMRecordDefaultFilter.accept(rec, this.min_mapq)) continue;
							int ref=rec.getStart();
							final Cigar cigar = rec.getCigar();
							if(cigar==null) continue;
							
							if(this.alpha_arc > 0.0 &&
								rec.getReadPairedFlag() && 
								!rec.getMateUnmappedFlag() && 
								!rec.getProperPairFlag() && 
								rec.getReferenceIndex().equals(rec.getMateReferenceIndex()) && 
								(!this.arc_only_rr_ff || (rec.getReadNegativeStrandFlag()==rec.getMateNegativeStrandFlag()))
								)
								{
								longest_arc = Math.max(longest_arc, Math.abs(rec.getMateAlignmentStart()-rec.getAlignmentStart()));
								
								final double xstart = pos2pixel.applyAsDouble(rec.getAlignmentStart());
								final double xend = pos2pixel.applyAsDouble(rec.getMateAlignmentStart());
								final double len = (xend - xstart);
								
								if(Math.abs(len)>this.min_invert && Math.abs(len)<this.max_invert ) {
									final double y2 = y_mid + (drawAbove?-1:1)*Math.min(y_mid,Math.abs(len/2.0));
									final Composite oldComposite = g2.getComposite();
									g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)this.alpha_arc));
									g2.setColor(getColor("arc",Color.ORANGE));
									final GeneralPath curve = new GeneralPath();
									curve.moveTo(xstart,y_mid);
									curve.quadTo(xstart + len/2.0, y2, xend, y_mid);
									g2.draw(curve);
									g2.setComposite(oldComposite);
									drawAbove = !drawAbove;
									}
								}
							
							for(CigarElement ce:cigar) {
								final CigarOperator op = ce.getOperator();
								final int len = ce.getLength();
								if(op.consumesReferenceBases()) {
									if(op.consumesReadBases()) {
										for(int i=0;i< len;i++) {
											final int pos = ref+i;
											if(pos < extendedRegion.getStart()) continue;
											if(pos > extendedRegion.getEnd()) break;
											depth[pos-extendedRegion.getStart()]++;
										}
									}
									ref+=len;
								}
							}
						}// loop cigar
					}// end samItere
	
				
				if(extendedRegion.getLengthOnReference()>image.getWidth() &&
					extendedRegion.getLengthOnReference() > this.window_smooth_size) {
					//smooth
					final int bases_per_pixel = this.window_smooth_size;
					System.arraycopy(depth, 0, copy, 0, depth.length);
					for(int i=0;i< depth.length && bases_per_pixel>1;i++) {
						double t=0;
						int count=0;
						for(int j=Math.max(0, i-bases_per_pixel);j<=i+bases_per_pixel && j< depth.length;j++) {
							t+=copy[j];
							count++;
							}
						if(count==0) continue;
						depth[i]=(int)(t/count);
						}
					}
				
				discreteMedian.clear();
				for(int i=0;i< depth.length;i++) {
					if(depth[i]>this.max_depth) continue;
					if(blackListedPositions.get(i)) continue;
					discreteMedian.add(depth[i]);
				}

				final double median = discreteMedian.getMedian().orElse(0);
				if(median<=0) {
					final String msg = "Skipping "+sample +" "+extendedRegion+" because median is 0";
					LOG.warning(msg);
					manifest.println("#"+msg);
					continue;
				}
				
				Point2D max_position=null;
				double max_distance_to_1=0.0;
				final GeneralPath line = new GeneralPath();
				
				for(int x=0;x< image.getWidth();x++) {
					discreteMedian.clear();
					int pos1= (int)Math.floor(((x+0)/(double)image.getWidth())*depth.length);
					final int pos2= (int)Math.ceil(((x+0)/(double)image.getWidth())*depth.length);
					while(pos1 <= pos2 && pos1 < depth.length) {
						discreteMedian.add(depth[pos1]);
						pos1++;
						} 
					final double average = discreteMedian.getMedian().orElse(0);
					final double normDepth = (average/median);
					
					if(normDepth<=0.5) count_has_dp_le_0_5++;
					if(normDepth>=1.5) count_has_dp_ge_1_5++;
					
					final double y2 = normToPixelY.applyAsDouble(normDepth);
					double distance_to_1 = Math.abs(normDepth-1.0);
					if(distance_to_1 > 0.3 && (max_position==null || distance_to_1 > max_distance_to_1)) {
						max_distance_to_1 = distance_to_1;
						max_position = new Point2D.Double(x,y2);
						}
					if(x==0) {
						line.moveTo(x, y2);
						}
					else
						{
						line.lineTo(x, y2);
						}
					}
				
				g.setColor(max_distance_to_1<0.1?Color.lightGray:Color.GRAY);
				final Stroke oldStroke = g.getStroke();
				final Composite oldComposite = g.getComposite();
				g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,(float)this.alpha));
				g.setStroke(new BasicStroke(0.5f,BasicStroke.CAP_BUTT,BasicStroke.JOIN_ROUND));
				g.draw(line);
				g.setStroke(oldStroke);
				g.setComposite(oldComposite);
				
				if(max_position!=null) sample2maxPoint.put(sample,max_position);
				
				
				// fill manifest
				manifest.print(rawRegion.getContig());
				manifest.print("\t");
				manifest.print(rawRegion.getStart()-1);
				manifest.print("\t");
				manifest.print(rawRegion.getEnd());
				manifest.print("\t");
				manifest.print(extendedRegion.getStart()-1);
				manifest.print("\t");
				manifest.print(extendedRegion.getEnd());
				manifest.print("\t");
				manifest.print(outputFilename);
				manifest.print("\t");
				manifest.print(StringUtils.isBlank(theGeneName)?".":theGeneName);
				manifest.print("\t");
				manifest.print(sample);
				manifest.print("\t");
				manifest.print(count_has_dp_le_0_5);
				manifest.print("\t");
				manifest.print((int)((count_has_dp_le_0_5/(double)extendedRegion.getLengthOnReference())*100.0));
				manifest.print("\t");
				manifest.print(count_has_dp_ge_1_5);
				manifest.print("\t");
				manifest.print((int)((count_has_dp_ge_1_5/(double)extendedRegion.getLengthOnReference())*100.0));
				manifest.print("\t");
				manifest.print(median);
				manifest.print("\t");
				manifest.print(discreteMedian.getStandardDeviation().orElse(-99999));
				manifest.print("\t");
				manifest.print(longest_arc);
				manifest.println();
				}//end loop over samples
			
				
				
				
			}
		
			

		g2.dispose();
		g.drawImage(offscreen,0,0,null);
				
			
			
		g.setColor(Color.BLACK);
		
		final String label_samples = sampleNames.size()>10?"N="+StringUtils.niceInt(sampleNames.size()):String.join(";",sampleNames);
		final Set<String> all_genes = getGenes(extendedRegion).
				filter(G->G.getType().equals("gene")).
				map(G->G.getAttribute("gene_name")).
				filter(F->!StringUtils.isBlank(F)).
				collect(Collectors.toCollection(TreeSet::new))
				;

		g.drawString(extendedRegion.toNiceString()+" Length:"+StringUtils.niceInt(extendedRegion.getLengthOnReference())+
				" Sample(s):"+label_samples+" "+label+" Gene(s):"+
				(all_genes.size()>20?"N="+StringUtils.niceInt(all_genes.size()):String.join(";", all_genes)),
				10,10
				);

		if(!sample2maxPoint.isEmpty())
			{
			/** draw sample names */
			g.setColor(getColor("sample",Color.BLUE));
			final int sampleFontSize = 7;
			final Font oldFont = g.getFont();
			g.setFont(new Font(oldFont.getName(), Font.PLAIN, sampleFontSize));
			for(final String sample:sample2maxPoint.keySet()) {
				final Point2D pt = sample2maxPoint.get(sample);
				double sny = pt.getY();
				if(sny>y_mid) sny+=sampleFontSize;
				g.drawString(sample, (int)pt.getX(),
						(int)Math.min(this.dimension.height-sampleFontSize,Math.max(sampleFontSize,sny))
						);
				}
			g.setFont(oldFont);
			}
		
		g.dispose();
		
		
		
		try(OutputStream out=archive.openOuputStream(outputFilename)){
			ImageIO.write(image, "PNG", out);
			out.flush();
			}
		
		}// end while iter
		archive.close();
		archive=null;
		manifest.flush();
		manifest.close();
		manifest=null;
		return 0;
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(manifest);
		CloserUtil.close(archive);
		}
	}

public static void main(final String[] args) {
	new CoveragePlotter().instanceMainWithExit(args);
	}

}
