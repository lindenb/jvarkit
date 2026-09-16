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
package com.github.lindenb.jvarkit.tools.bam2bedgraph;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.SampleSheet;
import com.github.lindenb.jvarkit.io.SampleSheetFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StopWatch;

/**
 BEGIN_DOC
 
```
$ java -jar dist/jvarkit.jar bams2bedgraph -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 2> /dev/null   | head
RF01	0	7	0.2
RF01	7	10	0.4
RF01	10	11	1.2
RF01	11	19	1.4
RF01	19	25	2.2
RF01	25	26	2.4
RF01	26	31	2.6
RF01	31	37	2.8
RF01	37	43	3.6
RF01	43	54	4.0

$ java -jar dist/jvarkit.jar bams2bedgraph --interval "RF02:1-10" -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 
RF02	0	6	0.0
RF02	6	9	0.4
RF02	9	10	0.6

$ find src/test/resources/ -type f -name "S*.bam" > jeter.list
$ java -jar dist/jvarkit.jar bams2bedgraph  -R src/test/resources/rotavirus_rf.fa  jeter.list 2> /dev/null | head
RF01	0	7	0.2
RF01	7	10	0.4
RF01	10	11	1.2
RF01	11	19	1.4
RF01	19	25	2.2
RF01	25	26	2.4
RF01	26	31	2.6
RF01	31	37	2.8
RF01	37	43	3.6
RF01	43	54	4.0

# with a samplesheet with a header and a depth column to normalize on the mean depth

$ cat jeter.csv 
bam,depth
src/test/resources/S2.bam,1.0
src/test/resources/S5.bam,2.0
src/test/resources/S3.bam,10.0
src/test/resources/S1.bam,5.3
src/test/resources/S4.bam,2.8

$ java -jar dist/jvarkit.jar bams2bedgraph  -R src/test/resources/rotavirus_rf.fa  jeter.csv 2> /dev/null  |head
RF01	0	7	0.03773585
RF01	7	10	0.0754717
RF01	10	11	0.43320754
RF01	11	19	0.4709434
RF01	19	25	0.9109434
RF01	25	26	0.982372
RF01	26	31	1.0201079
RF01	31	37	1.0915364
RF01	37	43	1.5315363
RF01	43	54	1.6070081

```


 END_DOC 
 */

@Program(name="bams2bedgraph",
description="conver multiple BAM coverage to BedGraph",
keywords={"bam","alignment","bed","bedgraphs"},
creationDate = "20260916",
modificationDate = "20260916",
jvarkit_amalgamion =  true
)
public class BamsToBedGraph extends Launcher {
    private static final Logger LOG = Logger.of(BamsToBedGraph.class);
	@Parameter(names={"-B","--bed"},description="restrict to that input bed (or use --interval or reference dict)")
	protected Path bedPath = null;
	@Parameter(names={"--interval","--regions"},description="restrict to that interval chrom:start-end (or use --bed or reference dict)")
	protected String regionStr = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	protected Path referenceFilePath = null;
	@Parameter(names={"--mapq","-Q"},description="min mapping quality")
	protected int mapq=1;
	@Parameter(names={"--per_base"},description="output every base in each interval, do not group positions with the same depth")
	protected boolean per_base=false;
	@Parameter(names={"--no-normalize"},description="do not divide by the number of bases (cumulative number of bases)")
	protected boolean do_not_normalize=false;
	@Parameter(names={"--output","-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile =null;

	private static class BamInput {
		final Path bam;
		double depth=1.0;
		BamInput(final Path bam) {
			this.bam=bam;
		}
	}
	
	
@Override
public int doWork(final List<String> args) {

	try {
		final List<BamInput> bam_inputs;

		if(args.isEmpty() || (args.size()==1 && args.get(0).equals("-"))) {
			try(BufferedReader br= IOUtils.openStreamForBufferedReader(stdin())) {
				bam_inputs = br.lines().filter(L->!(L.startsWith("#") || StringUtils.isBlank(L)))
					.map(L->new BamInput(Paths.get(L)))
					.collect(Collectors.toList());
				}
			}
		else  if(args.size()==1 && args.get(0).endsWith(".list")) {
			try(BufferedReader br= IOUtils.openPathForBufferedReading(Paths.get(args.get(0)))) {
				bam_inputs = br.lines().filter(L->!(L.startsWith("#") || StringUtils.isBlank(L)))
						.map(L->new BamInput(Paths.get(L)))
					.collect(Collectors.toList());
				}
			}
		else  if(args.size()==1 && (args.get(0).endsWith(".csv") || args.get(0).endsWith(".tsv"))) {
			final Path sheetPath = Paths.get(oneAndOnlyOneFile(args));
			final SampleSheet samplesheet ;
			
			samplesheet = new SampleSheetFactory()
						.splitter(CharSplitter.forPath(sheetPath))
						.of(sheetPath);
			samplesheet.getHeader().assertColumnExists("bam");
			bam_inputs=new ArrayList<>(samplesheet.size());
			for(FileHeader.RowMap row: samplesheet) {
				final BamInput bi=new BamInput(Paths.get(row.get("bam")));
				if(row.containsKey("depth")) {
					bi.depth = Double.parseDouble(row.get("depth"));
					}
				bam_inputs.add(bi);
				}
			}
		else
			{
			bam_inputs=new ArrayList<>(args.size());
			for(String arg: args ) {
				bam_inputs.add(new BamInput(Paths.get(arg)));
				}
			}
		
		if(bam_inputs.isEmpty()) {
			LOG.error("no bam defined");
			return -1;
			}
		
		for(BamInput bi:bam_inputs) {
			IOUtil.assertFileIsReadable(bi.bam);
            if(bi.depth<=0) throw new IllegalArgumentException(bi.bam+" "+bi.depth);
			}
		
		
		final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.referenceFilePath);

		final List<Locatable> intervals = new ArrayList<>();
		
		if(this.bedPath!=null && !StringUtils.isBlank(this.regionStr)) {
			LOG.error("both regions and interval defined");
			return -1;
			}
		else if(!StringUtils.isBlank(this.regionStr)) {
			final Locatable interval =  new IntervalParser(refDict).apply(this.regionStr).orElseThrow(IntervalParser.exception(this.regionStr));
			intervals.add(interval);
			}
		else if(this.bedPath!=null) {
			try(BedLineReader br=new BedLineReader(this.bedPath)) {
				br.setContigNameConverter(ContigNameConverter.fromOneDictionary(refDict));
				while(br.hasNext()) {
					intervals.add(br.next());
					}
				}
			intervals.sort(new ContigDictComparator(refDict).createLocatableComparator());
			int i=0;
			while(i+1< intervals.size()) {
				if(intervals.get(i).overlaps(intervals.get(i+1))) {
					intervals.set(i, new SimpleInterval(
						intervals.get(i).getContig(),
						intervals.get(i).getStart(),
						Math.max(intervals.get(i).getEnd(),intervals.get(i+1).getEnd())
						));
					intervals.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			}
		else
			{
			intervals.addAll(refDict.getSequences());
			}
		
		if(intervals.isEmpty())  {
			LOG.error("not interval defined");
			return -1;
		}
		

		final SamReaderFactory srf = super.createSamReaderFactory();
		srf.referenceSequence(this.referenceFilePath);
		final long n_expected = intervals.size()*(long)bam_inputs.size();
		long n_count=0;
		try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {

			for(Locatable interval: intervals) {
				System.gc();
				final double[] coverage= new double[interval.getLengthOnReference()];
				Arrays.fill(coverage, 0.0);
				final int[] tmp_coverage= new int[interval.getLengthOnReference()];
				
		        final StopWatch stopWatch = new StopWatch();
		        stopWatch.start();
		
				for(BamInput row :bam_inputs) {
					Arrays.fill(tmp_coverage, 0);
		            LOG.info(row.bam+" "+ interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd() +" "+ StringUtils.niceInt(n_count)+"/"+StringUtils.niceInt(n_expected)+
		                    " remains "+(n_count==0?"xxxxx":""+(int)(((n_expected-n_count)*(stopWatch.getElapsedTimeSecs()/(double)n_count))))+" secs");
		            
		            try(SamReader sr= srf.open(row.bam)) {
		            	try(CloseableIterator<SAMRecord> iter= sr.query(interval.getContig(), interval.getStart(), interval.getEnd(), false)) {
		            		while(iter.hasNext()) {
		            			final SAMRecord rec= iter.next();
		            			if(!SAMRecordDefaultFilter.accept(rec,this.mapq)) continue;
		            			for(AlignmentBlock ab: rec.getAlignmentBlocks()) {
		            				for(int x=0;x<ab.getLength();++x) {
		            					final int pos = (ab.getReferenceStart()+x)- interval.getStart();
		            					if(pos<0) continue;
		            					if(pos>= tmp_coverage.length) break;
		            					tmp_coverage[pos]++;
		            					}
		            				}
		            			}
		            		}
		            	}
		            for(int x=0;x< tmp_coverage.length;++x) {
		            	coverage[x] +=  (tmp_coverage[x] / row.depth);
		            	}
		            ++n_count;
					}
				if(!this.do_not_normalize) {
					for(int x=0;x< coverage.length;++x) {
						coverage[x]/=bam_inputs.size();
			        	}
					}
				
				int i=0;
				while(i< coverage.length) {
					int j=i+1;
					while(!this.per_base && j< coverage.length && (float)coverage[i]==(float)coverage[j]) {
						++j;
						}
					pw.print(interval.getContig());
					pw.print("\t");
					pw.print((interval.getStart()-1)+i);// bed coordinate
					pw.print("\t");
					pw.print((interval.getStart()-1)+j);
					pw.print("\t");
					pw.print((float)coverage[i]);
					pw.println();
					i=j;
					}
					
				}
			pw.flush();
			}//end of writer
		return 0;
		}
	catch(final Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}
public static void main(String[] args) {
	new BamsToBedGraph().instanceMainWithExit(args);
}

}
