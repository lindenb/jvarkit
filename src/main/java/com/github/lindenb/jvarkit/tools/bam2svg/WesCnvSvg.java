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
package com.github.lindenb.jvarkit.tools.bam2svg;

import java.awt.Dimension;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.RatioConverter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.XLINK;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.svg.SVG;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Input

input is a set of bam file or a file with the '*.list' suffix containing the path to the bam files.

## Example

```
$ find dir -name "*.bam"  > bam.list
$ java -jar dist/wescnvsvg.jar -R ref.fasta -B cnv.bed bam.list > cnv.svg 
```

## Screenshots

https://twitter.com/yokofakun/status/1022503372669300738 : 

![ScreenShot](https://pbs.twimg.com/media/DjCpKcYXgAAq4fw.jpg:large)

https://twitter.com/yokofakun/status/1022805656905150464

![ScreenShot](https://pbs.twimg.com/media/DjG8Do0XsAA4U46.jpg:large)

https://twitter.com/yokofakun/status/1023953180927963137

![ScreenShot](https://pbs.twimg.com/media/DjXP8_4X0AAtSQZ.jpg)

https://twitter.com/yokofakun/status/1024315948847849472

![ScreenShot](https://pbs.twimg.com/media/DjcZ6vNXcAEYjEt.jpg)

https://twitter.com/yokofakun/status/1024315948847849472

![ScreenShot](https://pbs.twimg.com/media/DjcZ6vNXcAEYjEt.jpg)

https://twitter.com/yokofakun/status/1025330308193779712

![ScreenShot](https://pbs.twimg.com/media/Djq0Se-W0AEAbyR.jpg)

https://twitter.com/yokofakun/status/1040592885786263554

![ScreenShot](https://pbs.twimg.com/media/DnDttNgX4AAtxax.jpg)

https://twitter.com/yokofakun/status/1040577235856580608

![ScreenShot](https://pbs.twimg.com/media/DnDfaGLXcAArg0P.jpg)

https://twitter.com/yokofakun/status/1057625407913111557

![ScreenShot](https://pbs.twimg.com/media/Dq1whOTX0AAzkZc.jpg)

https://twitter.com/yokofakun/status/1180046139502059521

![ScreenShot](https://pbs.twimg.com/media/EGBdtO7WoAEHjEE?format=jpg&name=small)


END_DOC
 */
@Program(name="wescnvsvg",
description="SVG visualization of bam DEPTH for multiple regions",
keywords={"bam","alignment","graphics","visualization","svg","wes","bed","capture","exome"},
modificationDate="20200622",
creationDate="20180726"
)
public class WesCnvSvg  extends Launcher {
	private static final Logger LOG = Logger.build(WesCnvSvg.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-B","--bed","-b","--capture","-r","-rgn","--region","--interval"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();
	@Parameter(names={"-R","--ref","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-w","--width"},description="Page width")
	private int drawinAreaWidth = 1000 ;
	@Parameter(names={"-height","--height"},description="Sample Track height")
	private int sampleTrackHeight = 100 ;
	@Parameter(names={"-smooth","--smooth"},description="Smoothing pixel window size. Negative=don't smooth")
	private int pixSmoothSize = 100 ;
	@Parameter(names={"-cap","--cap"},description="Cap coverage to this value. Negative=don't set any limit")
	private int capMaxDepth = -1 ;
	@Parameter(names={"--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter samRecordFilter = SamRecordFilterFactory.ACCEPT_ALL;
	@Parameter(names={"--title"},description="document title")
	private String domSvgTitle=WesCnvSvg.class.getSimpleName();
	@Parameter(names={"-u","--url","--hyperlink"},description= "creates a hyperlink an area is 'clicked'. " + Hyperlink.OPT_DESC,converter=Hyperlink.StringConverter.class,splitter=NoSplitter.class)
	private Hyperlink hyperlinkType = Hyperlink.empty();
	@Parameter(names={"-p","-percentile","--percentile"},description="How to compute the percentil of a region")
	private Percentile.Type percentile = Percentile.Type.AVERAGE;
	@Parameter(names={"-css","--css"},description="custom svg css stylesheet")
	private File cssFile = null;
	@Parameter(names={"-x","--extend"},description="Extend each region by this factor. 100bp + 150% -> 150bp." + RatioConverter.OPT_DESC,converter=RatioConverter.class,splitter=NoSplitter.class)
	private double extendFactor= 1.0;
	
	
	private class BamInput implements Closeable
		{
		int index;
		Path bamPath;
		SamReader samReader=null;
		String sample;
		ContigNameConverter contigNameConverter;
		@SuppressWarnings("unused")
		double minDepth=0;
		double maxDepth=0;
		@Override
		public void close() throws IOException {
			CloserUtil.close(samReader);
			}
		double getPixelHeight() {
			return WesCnvSvg.this.sampleTrackHeight;
		}
		}
	private static class SampleInfo
		{
		double pixel_coverage[];
		double pixel_clipping[];
		}
	
	private class CaptureInterval
		implements Locatable
		{
		final QueryInterval queryInterval;
		final List<SampleInfo> sampleInfos = new ArrayList<>(bamInputs.size());
		double pixelx=0.0;
		CaptureInterval(final QueryInterval queryInterval)
			{
			this.queryInterval = queryInterval;
			}
		@Override
		public int getStart() {
			return queryInterval.start;
			}
		@Override
		public int getEnd() {
			return queryInterval.end;
			}
		
		public int getBaseLength() {
			return 1+(getEnd()-getStart());
		}
		
		@Override
		public String getContig() {
			return WesCnvSvg.this.refDict.getSequence(this.queryInterval.referenceIndex).getSequenceName();
			}
		public double getPixelX1() {
			return pixelx;
			}
		
		public double getPixelWidth(){
			return (this.getBaseLength()/(double)WesCnvSvg.this.countBasesToBeDisplayed)*WesCnvSvg.this.drawinAreaWidth;
		}
		
		String getName() {
			return this.getContig()+":"+StringUtils.niceInt(this.getStart())+"-"+StringUtils.niceInt(this.getEnd());
			}
		String getId() {
			return String.valueOf(this.queryInterval.referenceIndex)+"_"+this.getStart()+"_"+this.getEnd();
			}
		}
	
	private final List<CaptureInterval> intervals = new ArrayList<>();
	private final List<BamInput> bamInputs = new ArrayList<>();
	private ReferenceSequenceFile indexedFastaSequenceFile;
	private SAMSequenceDictionary refDict; 
	private DecimalFormat decimalFormater = new DecimalFormat("##.##");
	private double globalMaxDepth = 0.0;
	private double globalMaxClip = 0.0;
	private int countBasesToBeDisplayed = 0;
	private final int gc_win=100;

	private final DoubleUnaryOperator capDepthValue = (v)->
		 this.capMaxDepth<1 ?v:Math.min(this.capMaxDepth, v);
		 
	
	/** convert double to string */
	private String format(double v)
		{
		return this.decimalFormater.format(v);
		}
	
	private double getGcPercent(final GenomicSequence seq,int chromStart,int chromEnd)
		{
		final SAMSequenceRecord ssr = seq.getSAMSequenceRecord();
		while(chromEnd-chromStart+1< this.gc_win)
			{
			chromEnd++;
			chromStart--;
			}
		chromStart = Math.max(chromStart,1);
		chromEnd = Math.min(ssr.getSequenceLength(),chromEnd);
		if(chromStart>chromEnd) return 0.0;
		
		final GenomicSequence.GCPercent gcpercent=seq.getGCPercent(chromStart-1, chromEnd);
		if(gcpercent.isEmpty()) return 0.0;
		return gcpercent.getGCPercent();
		}	
	
	@Override
	public int doWork(final List<String> args) {
		XMLStreamWriter w = null;
		BufferedReader r = null;
		OutputStream fout=null;
		if(this.extendFactor<1.0) {
			LOG.error("extend factor <1.0");
			return -1;
			}
		try
			{
			
			this.indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			this.refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);			
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(refDict);

			
			final List<SimpleInterval> userIntervals = this.intervalListProvider.
					stream().
					map(loc->contigNameConverter.convertToSimpleInterval(loc).orElseThrow(()->new JvarkitException.ContigNotFoundInDictionary(loc.getContig(),refDict))).
					sorted(new ContigDictComparator(this.refDict).createLocatableComparator()).
					collect(Collectors.toCollection(ArrayList::new))
					;
		
			if(userIntervals.isEmpty())
				{
				LOG.error("no interval or bed defined");
				return -1;
				}
			
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.faidx);
			for(final Path bamFile:IOUtils.unrollPaths(args)) {
				final BamInput bi = new BamInput();
				bi.index = this.bamInputs.size();
				bi.bamPath = bamFile;
				bi.samReader = srf.open(bamFile);
				final SAMSequenceDictionary samDict = bi.samReader.getFileHeader().getSequenceDictionary();
				bi.contigNameConverter = samDict==null?ContigNameConverter.getIdentity():ContigNameConverter.fromOneDictionary(samDict);
				
				JvarkitException.BamHasIndex.verify(bi.samReader);
				final SAMSequenceDictionary dict2= SequenceDictionaryUtils.extractRequired(bi.samReader.getFileHeader());
				if(!SequenceUtil.areSequenceDictionariesEqual(this.refDict, dict2))
					{
					LOG.warn("Not the same dictionaries ! REF/BAM "+ JvarkitException.DictionariesAreNotTheSame.getMessage(this.refDict, dict2));
					}
				
				bi.sample = bi.samReader.getFileHeader().getReadGroups().stream().
					map(V->V.getSample()).
					filter(S->!StringUtil.isBlank(S)).
					findFirst().
					orElse(bamFile.getFileName().toString());
				this.bamInputs.add(bi);
				}
			if(this.bamInputs.isEmpty()) {
				LOG.error("no bam input");
				return -1;
			}
			
			
			final List<QueryInterval> listQueryIntervals = new ArrayList<>(userIntervals.size());
			for(SimpleInterval interval: userIntervals)
				{
				final int tid = this.refDict.getSequenceIndex(interval.getContig());
				if(tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(interval.getContig(), refDict);
				
				//user asked to extend the regions
				if(this.extendFactor>1.0)
					{
					final SAMSequenceRecord ssr = Objects.requireNonNull(this.refDict.getSequence(tid),"??");
					final int mid = interval.getStart()+interval.getLengthOnReference()/2;
					final int len2= (int)(interval.getLengthOnReference()*this.extendFactor);
					final int start2= Math.max(1,mid-len2/2);
					final int end2 = Math.min(ssr.getSequenceLength(),mid+len2/2);

					interval = new SimpleInterval(interval.getContig(),start2,end2);
					}
				
				final QueryInterval qInterval = new QueryInterval(tid, interval.getStart(), interval.getEnd());
				listQueryIntervals.add(qInterval);
				}
			this.intervals.addAll( Arrays.stream(QueryInterval.optimizeIntervals(listQueryIntervals.toArray(new QueryInterval[listQueryIntervals.size()]))).
					map(Q->new CaptureInterval(Q)).
					collect(Collectors.toList()));
			this.countBasesToBeDisplayed = intervals.stream().
					mapToInt(R->1+(R.getEnd()-R.getStart())).
					sum();
			if(this.countBasesToBeDisplayed<1) {
				LOG.error("Nothing to display. BED count==0");
				return -1;
				}
			else
				{
				double x1=0;
				for(int i=0;i< this.intervals.size();++i) {
					final CaptureInterval ci = this.intervals.get(i);
					ci.pixelx = x1;
					x1+=ci.getPixelWidth();
					}
				}
			final Percentile thePercentile = Percentile.of(this.percentile);
			
			for(final CaptureInterval ci:this.intervals) {
				for(final BamInput bi:this.bamInputs)
					{
					final SampleInfo si = new SampleInfo();
					si.pixel_coverage = new double[(int)ci.getPixelWidth()];
					Arrays.fill(si.pixel_coverage, 0.0);
					si.pixel_clipping = new double[(int)ci.getPixelWidth()];
					Arrays.fill(si.pixel_clipping, 0.0);
					LOG.info("get cov "+ci.getName()+" for "+bi.bamPath);
					ci.sampleInfos.add(si);
					final int base_coverage[] = new int[ci.getBaseLength()];
					Arrays.fill(base_coverage, 0);
					final int clip_coverage[] = new int[ci.getBaseLength()];
					Arrays.fill(clip_coverage, 0);
					final String newContig = bi.contigNameConverter.apply(ci.getContig());
					if(newContig==null) {
						LOG.error("cannot find contig "+ci.getContig()+" in "+bi.bamPath);
						return -1;
						}
					final SAMRecordIterator iter=bi.samReader.queryOverlapping(newContig,ci.getStart(),ci.getEnd());
					while(iter.hasNext())
						{
						final SAMRecord rec = iter.next();
						if(rec.getReadUnmappedFlag()) continue;
						if(this.samRecordFilter.filterOut(rec)) continue;
						final Cigar cigar=rec.getCigar();
						if(cigar==null || cigar.isEmpty()) continue;
						int ref1=rec.getUnclippedStart();
						
						
						for(final CigarElement ce:cigar) {
							final CigarOperator op = ce.getOperator();
							
							
							if(op.isClipping())
								{
								for(int x=0;x< ce.getLength();++x){
									final int pos=ref1+x;
									if(pos< ci.getStart()) continue;
									if(pos> ci.getEnd()) break;
									clip_coverage[pos-ci.getStart()]++;
									}
								ref1 +=  ce.getLength();
								continue;
								}
							
							if(op.consumesReferenceBases())
								{
								if(op.consumesReadBases()){
									for(int x=0;x< ce.getLength();++x){
										final int pos=ref1+x;
										if(pos< ci.getStart()) continue;
										if(pos> ci.getEnd()) break;
										base_coverage[pos-ci.getStart()]++;
										}
									}
								ref1+=ce.getLength();
								}
							}
						}
					iter.close();
					
					
					for(int x=0;x< si.pixel_coverage.length;x++) {
						final int pos0 = Math.min(base_coverage.length, (int)(((x+0)/ci.getPixelWidth())*ci.getBaseLength()));
						final int pos1 = Math.min(base_coverage.length, (int)Math.ceil(((x+1)/ci.getPixelWidth())*ci.getBaseLength()));
						if(pos0>=pos1) continue;
						si.pixel_coverage[x] = thePercentile.evaluate(base_coverage,pos0,(pos1-pos0)).getAsDouble();
						si.pixel_clipping[x] = thePercentile.evaluate(clip_coverage,pos0,(pos1-pos0)).getAsDouble();
					}
					
					if(this.pixSmoothSize>0)
						{
						final double newcov[]=new double[si.pixel_coverage.length];
						for(int x=0;x<si.pixel_coverage.length;++x) {
							double sum=0;
							int count=0;
							for(int y=-this.pixSmoothSize;y<=this.pixSmoothSize;++y)
								{
								int array_index = x+y;
								if(array_index<0 || array_index>= si.pixel_coverage.length) continue;
								sum+=si.pixel_coverage[array_index];
								count++;
								}
							newcov[x]=(sum/count);
							}
						System.arraycopy(newcov, 0, si.pixel_coverage, 0, newcov.length);
						}
					/* debug
						{
						int max_index=-1;
						for(int z=0;z < si.coverage.length;z++)
							{
							if(max_index==-1 || si.coverage[z]>si.coverage[max_index])
								{
								max_index=z;
								}
							}
						System.err.println(bi.sample+" "+ci.getName()+" / "+
								Arrays.stream(si.coverage).max().getAsDouble()+" at"+(
										ci.getStart()+max_index)+"="+si.coverage[max_index]
												);
						}*/
					}
				}
			
			
			// compute min/max depth for each sample
			for(final BamInput bi:this.bamInputs)
				{
				bi.minDepth = this.intervals.stream().flatMapToDouble(CI->
						DoubleStream.of(CI.sampleInfos.get(bi.index).pixel_coverage)
						).
					map(capDepthValue).
					min().orElse(0);
				bi.maxDepth = this.intervals.stream().flatMapToDouble(CI->
					DoubleStream.of(CI.sampleInfos.get(bi.index).pixel_coverage)
					).
					map(capDepthValue).
					max().orElse(1);
				LOG.debug("Sample "+bi.sample+" Max-Depth:"+bi.maxDepth);
				}
			
			this.globalMaxDepth = Math.max(1.0,this.bamInputs.stream().
					mapToDouble(B->B.maxDepth).
					map(this.capDepthValue).
					max().orElse(0));
			LOG.debug("global max depth "+this.globalMaxDepth);
			
			this.globalMaxClip = Math.max(1.0,this.intervals.stream().
					flatMap(R->R.sampleInfos.stream()).
					flatMapToDouble(R->Arrays.stream(R.pixel_clipping)).
					map(this.capDepthValue).
					max().orElse(0));
			LOG.debug("global max clip "+this.globalMaxDepth);
			
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w=xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				fout =Files.newOutputStream(this.outputFile);
				w=xof.createXMLStreamWriter(fout, "UTF-8");
				}
			
			final Function<List<Point2D.Double>,String> points2str = (L)->
				L.stream().map(S->format(S.getX())+","+format(S.getY())).
				collect(Collectors.joining(" "));
				
			final Consumer<List<Point2D.Double>> simplifyPoints = (L)->{
				for(int z=0;z+1< L.size();++z) {
					if(L.get(z).getY()==L.get(z+1).getY())
						{
						L.get(z).x = L.get(z+1).x;
						L.remove(z+1);
						}
					}
			};
			
			w.writeStartDocument("UTF-8", "1.0");
			
			
			final Dimension dim=new Dimension(this.drawinAreaWidth,0);
			final int bed_header_height=20;
			dim.height += bed_header_height;
			dim.height+=(int)this.bamInputs.stream().
					mapToDouble(B->B.getPixelHeight()).
					sum();
			
			LOG.debug("drawing area: "+dim);
			
			w.writeStartElement("svg");
			w.writeAttribute("width", String.valueOf(dim.width));
			w.writeAttribute("height", String.valueOf(dim.height));
			w.writeDefaultNamespace(SVG.NS);
			w.writeNamespace("xlink", XLINK.NS);
			
			// https://stackoverflow.com/questions/15717970
			w.writeStartElement("style");
			if(this.cssFile!=null) {
				w.writeCharacters(IOUtil.slurp(this.cssFile));
				}
			else {
				w.writeCharacters(
						"g.maing {stroke:black;stroke-width:0.5px;fill:whitesmoke;font-size:10pt;}\n"+
						"text.sampleLabel {stroke:none;stroke-width:0.5px;fill:blue;}" +
						"text.captureLabel {stroke:none;stroke-width:0.5px;fill:slategrey;text-anchor:middle;}" +
						"polygon.area {stroke:darkgray;stroke-width:0.5px;fill:lightgray;}" +
						"line.linedp {stroke:darkcyan;stroke-width:0.3px;opacity:0.4;}" +
						"text.linedp {fill-opacity:0.6;font-size:7px;stroke:none;stroke-width:0.5px;fill:darkcyan;}" +
						"rect.sampleFrame { fill:none;stroke:slategray;stroke-width:0.3px;}" +
						"rect.clickRgn {fill:none;stroke:none;pointer-events:all;}" +
						"polyline.gc {stroke:lightcoral;stroke-width:0.3px;fill:none;}"+
						"polyline.clipping {stroke:orange;stroke-width:0.8px;fill:none;}"
						);
				}
			w.writeEndElement();//style
			
			w.writeStartElement("title");
			w.writeCharacters(this.domSvgTitle);
			w.writeEndElement();

			w.writeStartElement("defs");
			// gc percent
			for(final CaptureInterval ci:this.intervals)
				{
				final GenomicSequence genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile,ci.getContig());
				final int gc_percent_width= (int)ci.getPixelWidth();
				final List<Point2D.Double> points= new ArrayList<>(gc_percent_width);
				for(int x=0;x< gc_percent_width;++x)
					{
					int pos1= ci.getStart()+(int)(((x+0)/ci.getPixelWidth())*ci.getBaseLength());
					int pos2= ci.getStart()+(int)(((x+1)/ci.getPixelWidth())*ci.getBaseLength());
					double gc_percent = getGcPercent(genomicSequence,pos1,pos2);
					double y = this.sampleTrackHeight - this.sampleTrackHeight*gc_percent;
					
					points.add(new Point2D.Double(x, y));
					}
				simplifyPoints.accept(points);
				
				w.writeStartElement("polyline");
				w.writeAttribute("class","gc");
				w.writeAttribute("id","z"+ci.getId());
				w.writeAttribute("points",points2str.apply(points));
				w.writeStartElement("title");
				w.writeCharacters("GC %");
				w.writeEndElement();
				w.writeEndElement();
				}
			
			w.writeEndElement();//defs
			
			
			w.writeStartElement("script");
			
			final StringBuilder openBrowserFunction = new StringBuilder(
					"function openGenomeBrowser(contig,chromStart,chromEnd) {\n"
					);
			if(!this.hyperlinkType.isEmpty())
				{
				openBrowserFunction.append("var url=\""+this.hyperlinkType.getPattern()+"\".replace(/__CHROM__/g,contig).replace(/__START__/g,chromStart).replace(/__END__/g,chromEnd);\n");
				openBrowserFunction.append("window.open(url,'_blank');\n");

				}
			else
				{
				//nothing
				}
			openBrowserFunction.append("}\n");
			
			w.writeCData(
				openBrowserFunction.toString() +
				"function clicked(evt,contig,chromStart,chromEnd){\n" +
			    "    var e = evt.target;\n" +
			    "    var dim = e.getBoundingClientRect();\n" +
			    "    var x = 1.0 * evt.clientX - dim.left;\n" + 
			    "    var cLen = 1.0* (chromEnd - chromStart); if(cLen<1) cLen=1.0;\n" + 
			    "    var pos1 = chromStart + parseInt(((x+0)/dim.width)*cLen);\n" +
			    "    var pos2 = chromStart + parseInt(((x+1)/dim.width)*cLen);\n" +
			    "   openGenomeBrowser(contig,pos1,pos2);\n" +
			    "}\n");                
			w.writeEndElement();//script
			
			
			w.writeStartElement("g");
			w.writeAttribute("class", "maing");
			
			int y=0;
			w.writeStartElement("g");
			w.writeComment("interval background");
			for(final CaptureInterval ci:this.intervals)
				{
				w.writeStartElement("text");
					w.writeAttribute("class", "captureLabel");
					w.writeAttribute("textLength",String.valueOf(ci.getPixelWidth()*0.8));
					w.writeAttribute("lengthAdjust", "spacing");
					w.writeAttribute("x",String.valueOf(ci.getPixelX1()+ci.getPixelWidth()/2.0));
					w.writeAttribute("y",String.valueOf(bed_header_height-2));
					w.writeCharacters(ci.getName());
				w.writeEndElement();//text

				
				w.writeStartElement("rect");
					w.writeAttribute("style","stroke:black;fill:none;");
					w.writeAttribute("x",String.valueOf(ci.getPixelX1()));
					w.writeAttribute("y", "0");
					w.writeAttribute("width",String.valueOf(ci.getPixelWidth()));
					w.writeAttribute("height",String.valueOf(dim.height));
					w.writeStartElement("title");
						w.writeCharacters(ci.getName());
					w.writeEndElement();
				w.writeEndElement();
				}
			w.writeEndElement();//interval background
			
			y+= bed_header_height;
			
			for(final BamInput bi:this.bamInputs)
				{	
				w.writeComment(bi.bamPath.toString());
				w.writeStartElement("g");
				w.writeAttribute("transform","translate(0,"+y+")");
				
				
				
				
				
				for(final CaptureInterval ci:this.intervals)
					{
					final String clickedAttribute = "clicked(evt,\""+ci.getContig()+"\","+ci.getStart()+","+ci.getEnd()+")";
					final SampleInfo si = ci.sampleInfos.get(bi.index);
					final double leftx =ci.getPixelX1();
					w.writeStartElement("g");
					w.writeAttribute("transform","translate("+leftx+",0)");
					
					final int segment_width = (int)ci.getPixelWidth();

					//coverage
					{
					final List<Point2D.Double> points = new ArrayList<>(segment_width);
					points.add(new Point2D.Double(0,bi.getPixelHeight()));

					for(int px=0;px< si.pixel_coverage.length;px++)
						{
						final double y_avg_cov= this.capDepthValue.applyAsDouble(si.pixel_coverage[px]);
						final double new_y = bi.getPixelHeight()-(y_avg_cov/this.globalMaxDepth)*bi.getPixelHeight();
						points.add(new Point2D.Double(px,new_y));
						}
					points.add(new Point2D.Double(ci.getPixelWidth(),bi.getPixelHeight()));
					simplifyPoints.accept(points);				
					points.add(new Point2D.Double(leftx,bi.getPixelHeight()));//close
					w.writeEmptyElement("polygon");
					w.writeAttribute("class","area");
					//w.writeAttribute("onclick", clickedAttribute);
					w.writeAttribute("points",points2str.apply(points));
					}
					//w.writeEndElement();//g
					
					//clipping
					if(this.globalMaxClip>0)
					{
						final List<Point2D.Double> points = new ArrayList<>(segment_width);
						points.clear();
						points.add(new Point2D.Double(0,bi.getPixelHeight()));
						
						for(int px=0;px< si.pixel_clipping.length;px++)
							{
							final double y_avg_cov= this.capDepthValue.applyAsDouble(si.pixel_clipping[px]);
							final double new_y = bi.getPixelHeight()-(y_avg_cov/this.globalMaxClip)*bi.getPixelHeight();
							points.add(new Point2D.Double(px,new_y));
							}
						simplifyPoints.accept(points);
						points.add(new Point2D.Double(ci.getPixelWidth(),bi.getPixelHeight()));
						
						w.writeEmptyElement("polyline");
						w.writeAttribute("class","clipping");
						//w.writeAttribute("onclick", clickedAttribute);
						w.writeAttribute("points",points2str.apply(points));
						}
					
					
					
					int depthshift=1;
					for(;;) {
						final int numdiv =(int) (this.globalMaxDepth/depthshift);
						if(numdiv<=10) break;
						depthshift*=10;
						}
					
					
					int depth=depthshift;
					while(depth< bi.maxDepth)
						{
						double new_y = bi.getPixelHeight()-(depth/this.globalMaxDepth)*bi.getPixelHeight();
						
						w.writeStartElement("text");
							w.writeAttribute("class", "linedp");
							w.writeAttribute("x","1");
							w.writeAttribute("y",String.valueOf(new_y+1));
							w.writeCharacters(String.valueOf(depth));
						w.writeEndElement();//text
						
						w.writeStartElement("line");
							w.writeAttribute("class","linedp");
							w.writeAttribute("stroke-dasharray","4");
							w.writeAttribute("x1", "0");
							w.writeAttribute("x2", String.valueOf(ci.getPixelWidth()));
							w.writeAttribute("y1", String.valueOf(new_y));
							w.writeAttribute("y2", String.valueOf(new_y));
							w.writeStartElement("title");
								w.writeCharacters(String.valueOf(depth));
							w.writeEndElement();
						w.writeEndElement();//line
						depth+=depthshift;
						
						

						}
					// polyline
					w.writeEmptyElement("use");
					w.writeAttribute("xlink",XLINK.NS,"href","#z"+ci.getId());
					w.writeAttribute("x","0");
					w.writeAttribute("y","0");

					
					//click
					w.writeStartElement("rect");
					w.writeAttribute("class","clickRgn");
					w.writeAttribute("onclick", clickedAttribute);
					w.writeAttribute("x", "0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width", String.valueOf(ci.getPixelWidth()));
					w.writeAttribute("height", String.valueOf(bi.getPixelHeight()));
					w.writeEndElement();
					
					
					w.writeEndElement();//g
					}
				
				//frame for this sample
				w.writeStartElement("rect");
					w.writeAttribute("class","sampleFrame");
					w.writeAttribute("x","0");
					w.writeAttribute("y", "0");
					w.writeAttribute("width",String.valueOf(dim.width));
					w.writeAttribute("height",String.valueOf(bi.getPixelHeight()));
				w.writeEndElement();//rect

				
				w.writeStartElement("text");
					w.writeAttribute("class", "sampleLabel");
					w.writeAttribute("x","5");
					w.writeAttribute("y","12");
					w.writeStartElement("title");
					w.writeCharacters(bi.bamPath.toString());
					w.writeEndElement();
					w.writeCharacters(bi.sample);
				w.writeEndElement();//text

				
				w.writeEndElement();//g
				y+=bi.getPixelHeight();
				}
			
			w.writeStartElement("g");
			w.writeComment("interval lines");
			for(int n=0;n<=this.intervals.size();n++)
				{
				w.writeEmptyElement("line");
				String x1=
						n<this.intervals.size()?
						String.valueOf(this.intervals.get(n).getPixelX1()):
						String.valueOf(dim.width)
						;
				w.writeAttribute("x1", x1);
				w.writeAttribute("y1", "0");
				w.writeAttribute("x2", x1);
				w.writeAttribute("y2", String.valueOf(dim.height));
				}
			w.writeEndElement();//interval lines

			
			w.writeEndElement();//g
			w.writeEndElement();//svg
			w.writeEndDocument();
			w.flush();
			w.close();
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(fout);
			CloserUtil.close(r);
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(this.bamInputs);
			}
		}

	
public static void main(final String[] args)
	{
	new WesCnvSvg().instanceMainWithExit(args);
	}
}
