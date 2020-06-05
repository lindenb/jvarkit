/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
import java.awt.Color;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.function.IntToDoubleFunction;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import javax.imageio.ImageIO;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
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
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
/**
BEGIN_DOC
input is an interval of a file source of interval (bed, vcf, gtf, interval_list , ,etc...)


```
java -jar dist/depthanomaly.jar -R src/test/resources/rotavirus_rf.fa -B src/test/resources/S1.bam -B src/test/resources/S2.bam "RF01:1-4000" -w 50 | less -r
```


END_DOC 
 */
@Program(
	name="coverageplotter",
	description="Find anomaly of depth in intervals+bams",
	keywords={"cnv","bam","depth","coverage"},
	creationDate="20200605",
	modificationDate="20200605",
	generate_doc=false
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
	@Parameter(names={"--dimension"},description = "Image Dimension. " + com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter.OPT_DESC, converter=com.github.lindenb.jvarkit.jcommander.converter.DimensionConverter.StringConverter.class,splitter=NoSplitter.class)
	private Dimension dimension = new Dimension(1000,300);
	@Parameter(names={"--extend","-x"},description = FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double extend=1.0;

	@Parameter(names= {"--gtf"},description="Optional Tabix indexed GTF file. Will be used to retrieve an interval by gene name, or to display gene names in a region.")
	private Path gtfFile = null;
	@Parameter(names= {"--known"},description="Optional Tabix indexed Bed or VCF file containing known CNV. Both types must be indexed.")
	private Path knownCnvFile = null;
	@Parameter(names= {"--prefix"},description="Image File Prefix.")
	private String prefix="";

	
	/** return a stream of interval of the known CNV overlapping the region */
	private Stream<Interval> getKnownCnv(final Locatable region) {
		if(this.knownCnvFile==null) return Stream.empty();
		final String fname=this.knownCnvFile.getFileName().toString();
		if(fname.endsWith(".bed.gz")) {
			try(TabixReader tbr = new TabixReader(this.knownCnvFile.toString())) {
				final TabixReader tbrfinal = tbr;
				final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
				final String ctg = cvt.apply(region.getContig());
				if(StringUtils.isBlank(ctg)) {
					return Stream.empty();
					}
				
				final TabixReader.Iterator iter = tbr.query(ctg,region.getStart(), region.getEnd());
				final BedLineCodec codec = new BedLineCodec();
				final AbstractIterator<Interval> iter2= new AbstractIterator<Interval>() {
					@Override
					protected Interval advance() {
						try {
							for(;;) {
								final String line = iter.next();
								if(line==null) return null;
								final BedLine bed = codec.decode(line);
								if(bed==null) continue;
								return new Interval(region.getContig(),bed.getStart(),bed.getEnd(),false,bed.getOrDefault(3, ""));
								}
							} catch (final IOException e) {
							LOG.error(e);
							return null;
							}
						}
					};
				return StreamSupport.stream(new IterableAdapter<Interval>(iter2).spliterator(),false).onClose(()->{
					tbrfinal.close();
					});
				}
			catch(final Throwable err) {
				LOG.error(err);
				return Stream.empty();
				}
			}
		else if(FileExtensions.VCF_LIST.stream().anyMatch(X->fname.endsWith(X))) {
			VCFFileReader vcfFileReader= null;
			try {
				vcfFileReader = VCFReaderFactory.makeDefault().open(this.knownCnvFile,true);
				
				final ContigNameConverter cvt = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(vcfFileReader.getFileHeader()));
				final String ctg = cvt.apply(region.getContig());
				if(StringUtils.isBlank(ctg)) {
					vcfFileReader.close();
					return Stream.empty();
					}
				final VCFFileReader vcfFileReaderFinal = vcfFileReader;
				return vcfFileReader.query(ctg, region.getStart(), region.getEnd()).
						stream().
						filter(VC->!VC.isSNP()).
						map(VC->{
							final List<String> list = new ArrayList<>();
							if(VC.hasID()) list.add(VC.getID());
							if(VC.hasAttribute(VCFConstants.SVTYPE))  list.add(VC.getAttributeAsString(VCFConstants.SVTYPE,"."));
							return new Interval(region.getContig(),VC.getStart(),VC.getEnd(),false,String.join(";",list));}).
						onClose(()->vcfFileReaderFinal.close());
				}
			catch(final Throwable err) {
				if(vcfFileReader!=null) vcfFileReader.close();
				LOG.error(err);
				return Stream.empty();
				}
		}
		else
		{
			LOG.warn("not a vcf of bed.gz file "+this.knownCnvFile);
			return Stream.empty();
		}
	}
	

	
	
	private Stream<GTFLine> getGenes(final Locatable region) {
		if(this.gtfFile==null) return Stream.empty();
		try (TabixReader tbr = new TabixReader(this.gtfFile.toString())) {
			
			final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
			final String ctg = cvt.apply(region.getContig());
			if(StringUtils.isBlank(ctg)) {
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
			return Stream.empty();
			}
	}
	
	private void writeGenes(final Graphics2D g,final Rectangle rect,final Locatable region) {
		if(this.gtfFile==null) return;
		final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rect.getWidth();
		final Composite oldComposite = g.getComposite();
		final Stroke oldStroke = g.getStroke();
		g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.4f));
		g.setColor(Color.ORANGE);
		final double y= rect.getMaxY()-4.0;
		try {
			getGenes(region).
				filter(G->G.getType().equals("exon") || G.getType().equals("transcript")).
				forEach(feature->{
					final double x1 = position2pixel.applyAsDouble(feature.getStart());
					final double x2 = position2pixel.applyAsDouble(feature.getEnd());
					if(feature.getType().equals("exon") ) {
						g.draw(new Rectangle2D.Double(x1, y-1, (x2-x1), 3));
						}
					else if(feature.getType().equals("transcript") ) {
						g.draw(new Line2D.Double(x1, y, x2, y));
						}
					});

			}
		catch(Throwable err) {
			
			}
		finally {
			g.setComposite(oldComposite);
			g.setStroke(oldStroke);
			}
		}
	
	private void writeKnownCnv(final Graphics2D g,final Rectangle rectangle,final Locatable region) {
		if(this.knownCnvFile==null) return;
		final IntToDoubleFunction position2pixel = X->((X-region.getStart())/(double)region.getLengthOnReference())*rectangle.getWidth();
		final Composite oldComposite = g.getComposite();
		final Stroke oldStroke = g.getStroke();
		g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.4f));
		g.setColor(Color.MAGENTA);
		final double y= rectangle.getHeight()-8.0;
		try {
			this.getKnownCnv(region).forEach(R->{
				final double x1 = position2pixel.applyAsDouble(R.getStart());
				final double x2 = position2pixel.applyAsDouble(R.getEnd());
				g.draw(new Rectangle2D.Double(x1, y-1, Math.max(1.0,x2-x1), 3));
				});

			}
		catch(Throwable err) {
			}
		finally {
			g.setComposite(oldComposite);
			g.setStroke(oldStroke);
			}
		}
	
	
	private double median(final int array[]) {
		int len = array.length;
		Arrays.sort(array,0,len);
		while(len>0 && array[len-1]>=this.max_depth) {
			len--;
			}
		if(len==0) return 0;
		int mid_x = len/2;
		if(len%2==0) {
			return (array[mid_x-1]+array[mid_x])/2.0;
		} else {
			return array[mid_x];
		}	
	}

@Override
public int doWork(final List<String> args) {
	ArchiveFactory archive = null;
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
		final BufferedImage image = new BufferedImage(this.dimension.width,this.dimension.height,BufferedImage.TYPE_INT_RGB);
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
			
			final Graphics2D g= image.createGraphics();
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, this.dimension.width,this.dimension.height);
			int y =(int)(this.dimension.height/2.0);
			g.setColor(Color.BLUE);
			g.drawLine(0, y, image.getWidth(),y );
			y=(int)(this.dimension.height/4.0);
			g.setColor(Color.ORANGE);
			g.drawLine(0, y,image.getWidth(),y );
			y=(int)(3.0*this.dimension.height/4.0);
			g.drawLine(0, y,image.getWidth(),y );
			g.setColor(Color.DARK_GRAY);
			g.drawRect(0, 0, this.dimension.width-1,this.dimension.height-1);
			writeGenes(g, new Rectangle(0,0,image.getWidth(),image.getHeight()), extendedRegion);
			writeKnownCnv(g, new Rectangle(0,0,image.getWidth(),image.getHeight()), extendedRegion);
			if(this.extend>1) {
				g.setColor(Color.GREEN);
				int x = (int)(((rawRegion.getStart()-extendedRegion.getStart())/(double)extendedRegion.getLengthOnReference())*image.getWidth());
				g.drawLine(x, 0, x, image.getHeight() );
				x = (int)(((rawRegion.getEnd()-extendedRegion.getStart())/(double)extendedRegion.getLengthOnReference())*image.getWidth());
				g.drawLine(x, 0, x, image.getHeight() );
				}
			
			final int depth[]= new int[extendedRegion.getLengthOnReference()];
			final int copy[]= new int[depth.length];
			// smooth
			final int bases_per_pixel = (int)Math.ceil((double)depth.length/image.getWidth());


			for(final Path path: inputBams) {
				try(SamReader sr = samReaderFactory.open(path)) {
					final SAMFileHeader header= sr.getFileHeader();
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
	
				if(bases_per_pixel>1) {
					//smooth
					System.arraycopy(depth, 0, copy, 0, depth.length);
					for(int i=0;i< depth.length;i++) {
						double t=0;
						int count=0;
						for(int j=i-bases_per_pixel;j<=i+bases_per_pixel && j< depth.length;j++) {
							if(j<0) continue;
							t+=copy[j];
							count++;
							}
						if(count==0) continue;
						depth[i]=(int)(t/count);
						}
					}
				


				System.arraycopy(depth, 0, copy, 0, depth.length);
				final double median = median(copy);
				if(median<=0) continue;
				g.setColor(Color.GRAY);
				
				
				for(int i=0;i+1<depth.length;i++) {
					double x1 = (((i+0)/(double)depth.length))*image.getWidth();
					double x2 = (((i+1)/(double)depth.length))*image.getWidth();
					double y1= image.getHeight() - (depth[i+0]/median)*(image.getHeight()/2.0);
					double y2= image.getHeight() - (depth[i+1]/median)*(image.getHeight()/2.0);
					g.draw(new Line2D.Double(x1, y1, x2, y2));
					}
				}
			}
		
		g.setColor(Color.DARK_GRAY);
		g.drawString(extendedRegion.toNiceString()+" Length:"+StringUtils.niceInt(extendedRegion.getLengthOnReference())+" "+label, 10, 10);
		
		g.dispose();
			
		final String fname=prefix + extendedRegion.getContig()+"_"+extendedRegion.getStart()+"_"+extendedRegion.getEnd()+
				(StringUtils.isBlank(label)?"":"."+label.replaceAll("[^A-Za-z\\-\\.0-9]+", "_"))+".png";
		try(OutputStream out=archive.openOuputStream(fname)){
			ImageIO.write(image, "PNG", out);
			out.flush();
			}
		}// end while iter
		archive.close();
		archive=null;
		return 0;
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(archive);
		}
	}

public static void main(final String[] args) {
	new CoveragePlotter().instanceMainWithExit(args);
	}

}
