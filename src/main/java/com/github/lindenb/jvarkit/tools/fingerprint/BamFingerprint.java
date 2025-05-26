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
package com.github.lindenb.jvarkit.tools.fingerprint;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.PrefixSuffixWriter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StopWatch;

/**
BEGIN_DOC



## Example

```
find /path -type f -name "*.bam" > bams.list
java -jar dist/jvarkit.jar tssenrich -R ref.fa --gtf jeter.gtf bams.list > output.txt
```

## Output

output is a multipart text file. Each part can be isolated using awk. For example the following
commands plot the normalized peaks using R

```
$ awk '$1=="NORMALIZED"' output.txt | cut -f 2- > jeter.txt
$ awk '$1=="R_PLOT"' output.txt  | cut -f 2-  | sed 's/__INPUT__/jeter.txt/;s/__OUTPUT__/jeter.svg/' > jeter.R
$ R --vanilla < jeter.R
```

END_DOC
 
*/

@Program(name="bamfingerprint",
description="Bam fingerprint",
keywords={"bam","atacseq"},
modificationDate="20240202",
creationDate="20240202"
)
public class BamFingerprint extends Launcher {
	private static final Logger LOG = Logger.of(BamFingerprint.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path referenceFile = null;
    @Parameter(names={"--bin"},description="bin genomic size." + DistanceParser.OPT_DESCRIPTION, splitter = NoSplitter.class, converter  = DistanceParser.StringConverter.class)
    private int bin_size = 500;
    @Parameter(names={"--bin-reads"},description="bin read count. Group count of reads per group of 'n' reads")
    private int bin_read_count = 10;

    @Parameter(names={"--blacklist","--blacklist-bed"},description="Bed for regions that should be ignored.")
    private Path blackListBed = null;

    @Parameter(names={"--mapq"},description="min mapping quality")
    private int mapq = 1;
    @Parameter(names={"--threads"},description="number of parallel jobs")
    private int nThreads = 1;
    @Parameter(names={"--contig-regex"},description="use contigs matching this regex")
    private String contigRegexStr  = "(chr)?[0-9XY]+" ;

        
	
	private static class Summary {
		Throwable error=null;
		String sampleName;
		final Counter<Integer> numReadsPerBin = new Counter<>();
		
		}
	
	
	private static class Scanner implements Callable<Summary>  {
		protected final Path bamPath;
		protected final Path fastaPath;
		protected final Predicate<SAMRecord> samFilter;
		private final IntervalTreeMap<Boolean> blackList;
		private final int bin_window;
		private final int bin_reads;
		Scanner(
			final Path bamPath,
			final Path fastaPath,
			final Predicate<SAMRecord> samFilter,
			final IntervalTreeMap<Boolean> blackList,
			final int bin_window,
			final int bin_reads
			) {
			this.bamPath = bamPath;
			this.fastaPath = fastaPath;
			this.samFilter = samFilter;
			this.blackList = blackList;
			this.bin_window = bin_window;
			this.bin_reads = bin_reads;
			}
		
		protected SamReader open()  throws IOException {
			final SamReaderFactory srf = SamReaderFactory.makeDefault();
			srf.validationStringency(ValidationStringency.LENIENT);
			srf.referenceSequence(this.fastaPath);
			return srf.open(this.bamPath);
			}
		protected String getSample(final SAMFileHeader header) {
			return  header.getReadGroups().
					stream().
					map(RG->RG.getSample()).
					filter(S->!StringUtils.isBlank(S)).
					findFirst().
					orElse(IOUtils.getFilenameWithoutCommonSuffixes(this.bamPath));
			}

		private int normalize_read_count(int nReads) {
			return bin_reads*(int)(nReads/(double)bin_reads);
			}
		
		@Override
		public Summary call() throws Exception {
			final Summary summary = new Summary();
			final StopWatch stopWatch = new StopWatch();
			stopWatch.start();
			final SortedMap<Integer,Integer> bin2nreads=new TreeMap<>();
			
			String prevContig= SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
			try {
				try (SamReader sr = open()){
					final SAMFileHeader header =  sr.getFileHeader();
					summary.sampleName = getSample(header);
					try(CloseableIterator<SAMRecord> iter= sr.iterator()) {
						while(iter.hasNext()) {
							final SAMRecord record = iter.next();
							if(!this.samFilter.test(record)) continue;
							if(blackList!=null && blackList.containsOverlapping(record)) continue;
							
							if(!record.getContig().equals(prevContig)) {
								for(Integer count_reads : bin2nreads.values()) {
									summary.numReadsPerBin.incr(normalize_read_count(count_reads));
									}
								bin2nreads.clear();
								prevContig= record.getContig();
								}
							
							for(Iterator<Integer> iterator = bin2nreads.keySet().iterator(); iterator.hasNext(); ) {
								final int bin_idx = iterator.next();
								if(bin_idx*bin_window + bin_window >= record.getStart()) break;
								summary.numReadsPerBin.incr(normalize_read_count(bin2nreads.get(bin_idx)));
								iterator.remove();
								}
							
							Integer prev_bin = null;
							for(AlignmentBlock ab: record.getAlignmentBlocks()) {
								for(int n=0;n<ab.getLength();++n) {
									final int rec_pos = ab.getReferenceStart()+n-1/* 0 based */;
									final int bin_idx = rec_pos/this.bin_window;
									if(prev_bin!=null && prev_bin.intValue()==bin_idx) {
										//already filled
										continue;
										}
									prev_bin = bin_idx;
									
									Integer count = bin2nreads.get(bin_idx);
									bin2nreads.put(bin_idx, count==null?1:count.intValue()+1);
									}
								} // end align block
							}
						}
					for(Integer count_reads : bin2nreads.values()) {
						summary.numReadsPerBin.incr(normalize_read_count(count_reads));
						}
					}
				}
			catch(final Throwable error) {
				summary.error = error;
				}
			stopWatch.stop();
			LOG.info( summary.sampleName+" : That took:" + com.github.lindenb.jvarkit.lang.StringUtils.niceDuration(stopWatch.getElapsedTime()) );
			return summary;
			}
		}

	

	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final IntervalTreeMap<Boolean> blackListedMap;
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			final Pattern contigRegex = Pattern.compile(contigRegexStr);
			final Predicate<SAMRecord> samFilter = REC->SAMRecordDefaultFilter.accept(REC,this.mapq) && contigRegex.matcher(REC.getContig()).matches();

			
			if(bamPaths.isEmpty()) {
				LOG.error("input missing");
				return -1;
				}
			
			if(this.blackListBed==null) {
				blackListedMap = null;
				}
			else
				{
				try(final BedLineReader blr = new BedLineReader(this.blackListBed)) {
					if(this.referenceFile!=null) {
						blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(this.referenceFile)));
						}
					blackListedMap = new IntervalTreeMap<>();
					blr.stream().
						map(B->new Interval(B)).
						filter(R->contigRegex.matcher(R.getContig()).matches()).
						forEach(I->blackListedMap.put(I, Boolean.TRUE));
					}
				}
			


			
			final List<Scanner> bams = bamPaths.stream().
					map(B->new Scanner(B,referenceFile,samFilter,blackListedMap,this.bin_size,this.bin_read_count)).
					collect(Collectors.toList());
			
            final ExecutorService executorService = Executors.newFixedThreadPool(Math.min(this.nThreads,bams.size()));
			final List<Future<Summary>> future_summaries = executorService.invokeAll(bams);
            executorService.shutdown();
            executorService.awaitTermination(365, TimeUnit.DAYS);
            
			final List<Summary> summaries = new ArrayList<>(bamPaths.size());

            for(Future<Summary> f: future_summaries) {
            	if(!f.isDone() || f.isCancelled() || f.get()==null) {
            		throw new IllegalStateException("future not available ??");
            		}
            	
            	if(f.get().error!=null) {
            		LOG.error("got an error for "+ f.get().sampleName +" "+f.get().error);
            		return -1;
            		}
            	
            	if(f.get().numReadsPerBin.isEmpty()) {
            		LOG.warning("No read was found for "+ f.get().sampleName);
            		continue;
            		}
            	summaries.add(f.get());
            	}
            Collections.sort(summaries,(A,B)->A.sampleName.compareTo(B.sampleName));
			
            
            try(PrefixSuffixWriter pw = new PrefixSuffixWriter(super.openPathOrStdoutAsPrintWriter(output))) {
            	pw.setNoPrefix();
            	pw.println("# This is the output of "+getProgramName());
            	pw.setPrefix("GENERAL\t");
            	pw.println("date\t"+StringUtils.now());
            	pw.println("version\t"+getVersion());
            	pw.println("bin_size\t"+this.bin_size);
            	pw.println("fasta\t"+this.referenceFile);
            	pw.println("contig_regex\t"+contigRegex.pattern());
            	pw.println("exclude_bed\t"+this.blackListBed);
            	pw.println("command\t"+this.getProgramCommandLine());
            	
            	pw.setNoPrefix();
            	pw.println("# List of bams");
            	pw.setPrefix("FILE\t");
            	for(Path bamPath: bamPaths) {
            		pw.println(bamPath.toString());
            		}
            	
            	//
            	pw.setNoPrefix();
            	pw.println("# Raw outputs");
            	pw.setPrefix("RAW\t");
            	pw.println("sample\tcount_reads_min\tcount_reads_max\tcount_bins\thistogram");
            	for(final Summary summary: summaries) {
            		final long max_count= summary.numReadsPerBin.getMaxCount().orElse(1L);
	            	for(Map.Entry<Integer,Long> kv :  summary.numReadsPerBin.stream().
	            			sorted((A,B)->Integer.compare(A.getKey(), B.getKey())).
	            			collect(Collectors.toList())) {
	            		pw.print(summary.sampleName);
	            		pw.print("\t");
	            		pw.print(String.valueOf(kv.getKey()*this.bin_read_count));
	            		pw.print("\t");
	            		pw.print(String.valueOf((kv.getKey()+1)*this.bin_read_count));
	            		pw.print("\t");
	            		pw.print(String.valueOf(kv.getValue()));
	            		pw.print("\t");
	            		pw.print(StringUtils.repeat((int)(kv.getValue()/(double)max_count)*80, '#'));
	            		pw.println();
	            		}
	            	}
            	
            	// normalisation
            	pw.setNoPrefix();
            	pw.println("# normalized outputs");
            	pw.setPrefix("NORMALIZED\t");
            	final Set<Integer> all_count_or_reads = 
            			summaries.stream().
            			flatMap(S->S.numReadsPerBin.keySet().stream()).
            			collect(Collectors.toCollection(TreeSet::new));
            	pw.print("sample");
            	for(Integer n_reads: all_count_or_reads) {
            		pw.print("\t");
            		pw.print("["+(n_reads*this.bin_read_count) + "-"+((n_reads+1)*this.bin_read_count+"["));
            		}
            	pw.println();
            	for(final Summary summary: summaries) {
            		pw.print( summary.sampleName);
            		for(Integer n_reads: all_count_or_reads) {
                		pw.print("\t");
                		pw.print(String.valueOf(summary.numReadsPerBin.count(n_reads)));
                		}
            		pw.println();
            		}
            	
            	// R PLOT
            	pw.setNoPrefix();
            	pw.println("#\tR script used to plot the normalized results");
            	pw.println("#\t$ awk '$1==\"TABLE\"' output.txt | cut -f 2- > jeter.txt");
            	pw.println("#\t$ awk '$1==\"R_PLOT\"' output.txt  | cut -f 2-  | sed 's/__INPUT__/jeter.txt/;s/__OUTPUT__/jeter.svg/' > jeter.R");
            	pw.println("#\t$ R --vanilla < jeter.R");
            	pw.setPrefix("R_PLOT\t");
            	pw.write("data <- read.table(\"__INPUT__\", header = TRUE, sep = \"\\t\")\n");
            	pw.write("data_values <- data[, -1]\n");
            	pw.write("sample_name <- data[, 1]\n");
            	pw.write("y_max <- max(data_values)\n");
            	pw.write("sncolors <- rainbow(length(sample_name))\n");
            	pw.write("svg(\"__OUTPUT__\")\n");
            	pw.write("matplot(t(data_values), type = \"l\", col = sncolors, lty = 1, lwd = 2, ylim = c(0, y_max),\n");
            	pw.write("        xlab = \"Distance to TSS\", ylab = \"Score\",  xaxt='n',las=2,main=\"TSS Enrichment Score Summary\",sub=\""+this.referenceFile.getFileName().toString()+"\")\n");
            	pw.write("abline(v=ncol(data_values)/2.0,col=\"gray\")\n");
            	pw.write("legend(\"topright\", legend = sample_name, col = sncolors, lty = 1, lwd = 2, cex = 0.8)\n");
            	pw.write("dev.off()\n");

            	
            	pw.setPrefix("HIST\t");
            	pw.println("sample\tindex\tf\tstatus\thistogram");
            	for(final Summary summary: summaries) {
            		
            		for(Map.Entry<Integer,Long> kv :  summary.numReadsPerBin.stream().
	            			sorted((A,B)->Integer.compare(A.getKey(), B.getKey())).
	            			collect(Collectors.toList())) {
							pw.print(summary.sampleName+ "\t"+((i-array2.length/2)*winsize)+"\t"+v+"\t"+  score2status.apply(v) + "\t");

            				}
            		
            		final double genome_mean_cov = summary.sum_depth/(double)(adjusted_reference_length);
            		final double[] array2 = bin_array(summary.array,winsize);
	            	final double max=Arrays.stream(array2).map(V->V/genome_mean_cov).max().orElse(1.0);
					for(int i=0;i< array2.length;i++) {
						final double v = (array2[i]/genome_mean_cov);
						pw.print(summary.sampleName+ "\t"+((i-array2.length/2)*winsize)+"\t"+v+"\t"+  score2status.apply(v) + "\t");
						for(int x=0;x< (v/max)*80;++x) {
							pw.print("#");
							}
						pw.print("\n");
						}
	            	}
            	
            	pw.setNoPrefix();
            	pw.println("# work in progress do not use.");
            	pw.setPrefix("MULTIQC\t");
            	pw.println("plot_type: \"linegraph\"");
            	pw.println("description: \"multiqc output is not ready\"");
            	pw.println("pconfig:");
    			pw.println("  title: \"TSS Enrichment Score Summary\"");
    			pw.println("  xlab: \"Distance to TSS\"");
    			pw.println("  ylab: \"Score\"");
            	pw.println("data:");
            	for(final Summary summary: summaries) {
            		pw.println("  \"" + summary.sampleName+"\":");
            		}
            	
            	
            	pw.setNoPrefix();
            	pw.println("# EOF");
            	pw.flush();
            	}
            
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new BamFingerprint().instanceMainWithExit(args);

	}

}
