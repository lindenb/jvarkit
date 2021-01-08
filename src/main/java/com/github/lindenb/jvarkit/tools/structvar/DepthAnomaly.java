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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.ansi.AnsiUtils;
import com.github.lindenb.jvarkit.ansi.AnsiUtils.AnsiColor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.IndexCovUtils;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
/**
BEGIN_DOC
input is an interval of a file source of interval (bed, vcf, gtf, interval_list , ,etc...)

output is a text file that may contains unicode characters.

```
java -jar dist/depthanomaly.jar -R src/test/resources/rotavirus_rf.fa -B src/test/resources/S1.bam -B src/test/resources/S2.bam "RF01:1-4000" -w 50 | less -r
```

Note to self : use `less -r` to pipe

END_DOC 
 */
@Program(
	name="depthanomaly",
	description="Find anomaly of depth in intervals+bams",
	keywords={"cnv","bam","depth","coverage"},
	creationDate="20200605",
	modificationDate="20200605",
	generate_doc=false
	)
public class DepthAnomaly extends Launcher {
	private static final Logger LOG = Logger.build( DepthAnomaly.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"-B","--bams"},description = "list of bams. one file with a '.list' suffix is interpretted a a list of path to the bams",required=true)
	private List<String> bamsPath= new ArrayList<>();
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--length","-n"},description = "min anomaly length. Input interval must have length>= 'x'*3.",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_anomaly_length=100;
	@Parameter(names={"--treshold","-t"},description = IndexCovUtils.TRESHOLD_OPT_DESC)
	private double treshold = IndexCovUtils.DEFAULT_TRESHOLD;
	@Parameter(names={"--screen-width","-w"},description = "screen width. Disabled if <=0.")
	private int screen_width= 100;
	@Parameter(names={"--force-screen"},description = "For screen display even if there is no cnv.")
	private boolean force_screen = false;
	@Parameter(names={"--merge"},description = "merge intervals if they are withing '--length' bases.")
	private boolean merge_intervals = false;
	@Parameter(names={"--max-depth"},description = "ignore position if depth > 'x'")
	private int max_depth=500;

	private class CovInterval extends SimpleInterval {
		final List<Double> depths;
		final IndexCovUtils.SvType svtype;
			CovInterval(final String contig,int start,int end,IndexCovUtils.SvType svtype,List<Double> depths) {
				super(contig,start,end);
				this.svtype = svtype;
				this.depths=new ArrayList<>(depths); 
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
	PrintWriter pw = null;
	try
		{
		final IndexCovUtils indexCovUtils = new IndexCovUtils(this.treshold);
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
		final SamReaderFactory samReaderFactory = SamReaderFactory.
					makeDefault().
					referenceSequence(DepthAnomaly.this.refPath).
					
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
		 
		pw = super.openPathOrStdoutAsPrintWriter(this.outputFile);
		while(iter.hasNext()) {
			final SimpleInterval locatable = new SimpleInterval(iter.next());
			boolean found_anomaly_here = false;
			if(this.min_anomaly_length*3>=locatable.getLengthOnReference())  {
				LOG.warning("interval "+ locatable.toNiceString()+" is too short. skipping.");
				continue;
				}
			final int depth[]= new int[locatable.getLengthOnReference()];
			final int copy[]= new int[depth.length];
			for(final Path path: inputBams) {
				try(SamReader sr = samReaderFactory.open(path)) {
					final SAMFileHeader header= sr.getFileHeader();
					final String sample = header.getReadGroups().stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
					SequenceUtil.assertSequenceDictionariesEqual(dict,header.getSequenceDictionary());
					Arrays.fill(depth, 0);
					try(CloseableIterator<SAMRecord> siter = sr.queryOverlapping(locatable.getContig(), locatable.getStart(), locatable.getEnd())) {
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
											if(pos < locatable.getStart()) continue;
											if(pos > locatable.getEnd()) break;
											depth[pos-locatable.getStart()]++;
										}
									}
									ref+=len;
								}
							}
						}// loop cigar
					}// end samItere
				System.arraycopy(depth, 0, copy, 0, depth.length);
				final double median = median(copy);
				final List<CovInterval> anomalies = new ArrayList<>();
				//int minDepth = Arrays.stream(depth).min().orElse(0);	
				int x0=0;
				while(x0 < depth.length && median >0.0) {
					final int xi = x0;
					double total=0;
					double count=0;
					IndexCovUtils.SvType prevType = null;
					while(x0< depth.length ) {
						final IndexCovUtils.SvType type;
						final  int depthHere = depth[x0];
						final double normDepth = depthHere/(median==0?1.0:median);
						if(depthHere > this.max_depth) {
							type  = null;
							}
						else
							{
							type = indexCovUtils.getType(normDepth);
							}
						x0++;
						if(type==null || !type.isVariant()) break;
						if(prevType!=null &&  !type.equals(prevType)) break;
						if(prevType==null) prevType = type;
						total+= depthHere;
						count++;
						}
					if(prevType!=null && count  >= this.min_anomaly_length) {
						anomalies.add(new CovInterval(locatable.getContig(),
								locatable.getStart()+xi,
								locatable.getStart()+x0-1,
								prevType,
								Collections.singletonList(total/count)
								));
						}
					}
				if(!anomalies.isEmpty() || force_screen) {
					int i=0;
					while(i +1 < anomalies.size() && this.merge_intervals) {
						final CovInterval loc1= anomalies.get(i);
						final CovInterval loc2= anomalies.get(i+1);
						if(loc1.svtype.equals(loc2.svtype) && loc1.withinDistanceOf(loc2,this.min_anomaly_length)) {
							final List<Double> newdepths = new ArrayList<>(loc1.depths);
							newdepths.addAll(loc2.depths);
							anomalies.set(i, new CovInterval(loc1.getContig(),loc1.getStart(),loc2.getEnd(),loc1.svtype,newdepths));
							anomalies.remove(i+1);
							}
						else
							{
							i++;
							}
						}
					if(!found_anomaly_here) {
						pw.println(">>> "+ locatable.toNiceString()+" length:"+StringUtils.niceInt(locatable.getLengthOnReference()));
						found_anomaly_here = true;
						}
					if(screen_width>0) {
						pw.print("#");
						pw.print(String.format("%-15s", sample));
						pw.print("[");
						for(i=0;i< screen_width;i++) {
							double t=0;
							double n=0;
							final int x1= (int) (((i+0)/(double)screen_width)*depth.length);
							final int x2= (int) (((i+1)/(double)screen_width)*depth.length);
							for(int x3=x1;x3<=x2 && x3 < depth.length;++x3) {
								t+= depth[x1];
								n++;
								}
							
							double normDepth =  t/=n;
							if(median>0) normDepth/= median;
							normDepth/=2.0; //centered
							final boolean is_anomaly = anomalies.stream().
									anyMatch(R->CoordMath.overlaps(
											R.getStart(), R.getEnd(),
											locatable.getStart()+x1,locatable.getStart()+x2)
											);
							final AnsiUtils.AnsiColor color = is_anomaly?AnsiColor.RED:null;
							if(color!=null) pw.print(color.begin());
							pw.print(AnsiUtils.getHistogram(normDepth));
							if(color!=null) pw.print(color.end());
							}
						pw.print("]");
						pw.println();
						}
					for(i=0;i< anomalies.size();i++) {
						final CovInterval anomalie = anomalies.get(i);
						pw.print(anomalie.getContig());
						pw.print("\t");
						pw.print(anomalie.getStart()-1);
						pw.print("\t");
						pw.print(anomalie.getEnd());
						pw.print("\t");
						pw.print(anomalie.getLengthOnReference());
						pw.print("\t");
						pw.print(anomalie.svtype.name());
						pw.print("\t");
						pw.print(sample);
						pw.print("\t");
						pw.print(path);
						pw.print("\t");
						pw.print(i+1);
						pw.print("\t");
						pw.print(anomalies.size());
						pw.print("\t");
						pw.print(locatable.toNiceString());
						pw.print("\t");
						pw.print((int)median);
						pw.print("\t");
						pw.print((int)anomalie.depths.stream().mapToDouble(X->X.doubleValue()).average().orElse(0));
						pw.print("\t");
						pw.println();
						}
					}
				}
			}
			
			if(found_anomaly_here) {
				pw.println("<<< "+ locatable.toNiceString()+" length:"+
						StringUtils.niceInt(locatable.getLengthOnReference()));
				pw.println();
			}
		}// end while iter
		pw.flush();
		pw.close();
		pw=null;
		return 0;
		}
	catch(final Throwable err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(pw);
		}
	}

public static void main(final String[] args) {
	new DepthAnomaly().instanceMainWithExit(args);
	}

}
