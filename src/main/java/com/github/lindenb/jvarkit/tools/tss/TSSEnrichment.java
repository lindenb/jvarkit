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
package com.github.lindenb.jvarkit.tools.tss;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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
import com.github.lindenb.jvarkit.samtools.util.AbstractLocatable;
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
import htsjdk.samtools.util.StopWatch;

/**
BEGIN_DOC

## Transcription Start Site (TSS) Enrichment Score

https://www.encodeproject.org/data-standards/terms/#enrichment

> The TSS enrichment calculation is a signal to noise calculation.
> The reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 2000 bp in either direction (for a total of 4000bp).
> This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth. This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome) there should be an increase in signal up to a peak in the middle. 


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

@Program(name="tssenrich",
description="Transcription Start Site (TSS) Enrichment Score calculation",
keywords={"bam","atacseq","peak","tss"},
modificationDate="20240206",
creationDate="20240130"
)
public class TSSEnrichment extends Launcher {
	private static final Logger LOG = Logger.of(TSSEnrichment.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION, required = true)
	private Path referenceFile = null;
    @Parameter(names={"--extend","-x"},description="extend tss site by 'x' bases." + DistanceParser.OPT_DESCRIPTION, splitter = NoSplitter.class, converter  = DistanceParser.StringConverter.class)
    private int extend_tss = 2_000;
    @Parameter(names={"--gtf","-gtf"},description="GFF file to be used as a source of TSS start sites")
    private Path gtfPath = null;
    @Parameter(names={"--mapq"},description="min mapping quality")
    private int mapq = 1;
    @Parameter(names={"--use-transcript"},description="use 'transcript' type in GTF instead of 'gene'")
    private boolean use_transcript_type = false;
    @Parameter(names={"--strand"},description="use strand TSS information to swap the before/after position of the read according to the TSS strand.")
    private boolean use_strand_info = false;
    @Parameter(names={"--threads"},description="number of parallel jobs")
    private int nThreads = 1;
    @Parameter(names={"--bins"},description="in the graphical output, normalize coverage over TSS using 'x' regions")
    private int num_bins = 100;
    @Parameter(names={"--treshold"},description="comma separated tresholds to set a status concerning,acceptable/ideal see https://www.encodeproject.org/atac-seq/#standards")
    private String treshold_str  ="6,10";
    @Parameter(names={"--contig-regex"},description="use contigs matching this regex")
    private String contigRegexStr  = "(chr)?[0-9XY]+" ;    
    
    
	private static class TSS extends AbstractLocatable {
		private final String contig;
		private final int position;
		private final int extend;
		private final int strand_flag ;
		TSS(final String c,int pos,int strand_flag,final int extend) {
			this.contig = c;
			this.position = pos;
			this.extend = extend;
			this.strand_flag = strand_flag;
			}
		@Override
		public String getContig() {
			return contig;
			}
		public int getPosition() {
			return position;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof TSS)) return false;
			final TSS other= TSS.class.cast(obj);
			return getPosition()==other.getPosition() && getContig().equals(other.getContig());
			}
		@Override
		public int hashCode() {
			return getContig().hashCode()*31 + Integer.hashCode(getPosition());
			}
		
		@Override
		public int getStart() {
			return Math.max(1,getPosition() - extend);
			}
		
		@Override
		public int getEnd() {
			return getPosition() + extend;
			}
		
		int distanceTo(final SAMRecord rec) {
			if(this.overlaps(rec)) return 0;
			return Math.min(
					Math.abs(rec.getStart()-this.getPosition()),
					Math.abs(rec.getEnd()-this.getPosition())
					);
			}
		}
	
	private static class Summary {
		Throwable error=null;
		String sampleName;
		
		
		long count_reads = 0L;
		long count_reads_in_tss = 0L;
		double[] array;
		
		// coverage
		long sum_depth_not_in_tss = 0L;
		int count_tss_with_reads = 0;
		
		// max score, will be used to show the quality
		double max_ratio = 0;
		}
		
	private static class TSSScanner implements Callable<Summary>  {
		protected final Path bamPath;
		protected final Path fastaPath;
		protected final Predicate<SAMRecord> samFilter;
		private final IntervalTreeMap<TSS> intervalTreeMap;
		private final int window_size;
		TSSScanner(
			final Path bamPath,
			final Path fastaPath,
			final Predicate<SAMRecord> samFilter,
			final IntervalTreeMap<TSS> intervalTreeMap,
			final int window_size
			) {
			this.bamPath = bamPath;
			this.fastaPath = fastaPath;
			this.samFilter = samFilter;
			this.intervalTreeMap = intervalTreeMap;
			this.window_size = window_size;
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

		
		@Override
		public Summary call() throws Exception {
			final Summary summary = new Summary();
			summary.array = new double[window_size*2];
			Arrays.fill(summary.array, 0.0);
			final StopWatch stopWatch = new StopWatch();
			stopWatch.start();
			
			try {
				final Set<TSS> found_tss = new HashSet<>();
				try (SamReader sr = open()){
					final SAMFileHeader header =  sr.getFileHeader();
					summary.sampleName = getSample(header);
					try(CloseableIterator<SAMRecord> iter= sr.iterator()) {
						while(iter.hasNext()) {
							final SAMRecord record = iter.next();
							if(!this.samFilter.test(record)) continue;

							summary.count_reads++;
							TSS closest = null;
							for(TSS tss:this.intervalTreeMap.getOverlapping(record)) {
								if(closest==null || tss.distanceTo(record) < closest.distanceTo(record)) {
									closest = tss;
									}
								}
							if(closest!=null) {
								summary.count_reads_in_tss++;
								found_tss.add(closest);
								for(AlignmentBlock ab: record.getAlignmentBlocks()) {
									for(int n=0;n<ab.getLength();++n) {
										final int rec_pos = ab.getReferenceStart()+n;
										final int distance_to_tss = rec_pos - closest.getPosition() *  /* ok negative if read is before TSS */
												(closest.strand_flag==-1?-1:1) /* inverse sign , only when strand is taken into account*/
												;
										
										final int array_index = distance_to_tss + summary.array.length/2  /* middle of array */ ;
										if(array_index >=0 && array_index<summary.array.length ) {
											summary.array[array_index]++;
											}
										}
									} // end align block
								} // end closest!=null
							else 
								{
								// no read overlapping
								for(AlignmentBlock ab: record.getAlignmentBlocks()) {
									summary.sum_depth_not_in_tss  += ab.getLength();
									}
								}
							}
						}
					summary.count_tss_with_reads = found_tss.size();
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

	
	private double[] bin_array(double[] array,final int winsize) {
		final double[] array2 = new double[array.length/winsize];
		for(int i=0;i< array.length;i+=winsize) {
			int count=0;
			double sum=0;
			for(int x=0;x< winsize && i+x < array.length;++x) {
				count++;
				sum+= array[i+x];
				}
			array2[i/winsize]=sum/count;
			}
		return array2;
		}

	
	private static long count_N_in_reference(final Path ref,Pattern contigRegex) throws IOException {
		if(ref==null) return 0L;
		long n=0L;
		boolean in_autosome=false;
		try(BufferedReader br = IOUtils.openPathForBufferedReading(ref)) {
			String line;
			while((line=br.readLine())!=null) {
				if(line.startsWith(">")) {
					line=line.substring(1).split("[ \t]")[0];
					in_autosome = contigRegex.matcher(line).matches();
					continue;
					}
				for(int i=0;in_autosome && i< line.length();++i) {
					switch(line.charAt(i)) {
						case 'a':case 'A':
						case 't':case 'T':
						case 'g':case 'G':
						case 'c':case 'C':
						case '\r': case ' ':
							break;
						default: n++; break;
						}
					}
				}
			}
		return n;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final IntervalTreeMap<TSS> tssMap = new IntervalTreeMap<>();
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			
			final Pattern contigRegex = Pattern.compile(contigRegexStr);
			final Predicate<SAMRecord> samFilter = REC->SAMRecordDefaultFilter.accept(REC,this.mapq) && contigRegex.matcher(REC.getContig()).matches();

			
			if(bamPaths.isEmpty()) {
				LOG.error("input missing");
				return -1;
				}
			
			if(gtfPath==null) {
				LOG.error("GTF file missing");
				return -1;
				}
			
			final double treshold_low;
			final double treshold_high;
			final Function<Double, String> score2status;
			
			{
				final double[] tokens = Arrays.stream(CharSplitter.COMMA.split(treshold_str)).
						filter(S->!StringUtils.isBlank(S)).
						mapToDouble(Double::parseDouble).
						sorted().toArray();
				if(tokens.length!=2) {
					LOG.error("expected two int "+this.treshold_str);
					return -1;
					}
				treshold_low = tokens[0];
				treshold_high = tokens[1];
				score2status = (V)->{
					if(V> treshold_high) {
        				return ("IDEAL");
        				}
        			else if(V < treshold_low) {
        				return ("CONCERNING");
        				}
        			else
        				{
        				return ("ACCEPTABLE");
        				}
					};
			}
			
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceFile);
			if(SequenceDictionaryUtils.isGRCh38(dict) && (treshold_low!=5 || treshold_high!=7)) {
				LOG.warning("for hg38, recommandation for treshold are 5,7 (you: "+treshold_low+","+treshold_high+")");
				}
			else if(SequenceDictionaryUtils.isGRCh37(dict) && (treshold_low!=6 || treshold_high!=10)) {
				LOG.warning("for hg19, recommandation for treshold are 6,10 (you: "+treshold_low+","+treshold_high+")");
				}
			
			final long raw_autosome_length = dict.getSequences().stream().
					filter(SR->contigRegex.matcher(SR.getSequenceName()).matches()).
					mapToLong(SR->SR.getSequenceLength()).
					sum();
			if(raw_autosome_length==0) {
				LOG.info("no autosome found in "+this.referenceFile+" using regex "+ contigRegex.pattern());
				return -1;
				}
			
			try(BufferedReader br =  IOUtils.openPathForBufferedReading(this.gtfPath)) {
				String line;
				final GTFCodec gtfCodec = new GTFCodec();
				final Function<String,String> ctgConverter =  ContigNameConverter.fromOneDictionary(dict);
				while((line=br.readLine())!=null) {
					if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
					final GTFLine entry= gtfCodec.decode(line);
					if(entry==null || !entry.getType().equals(use_transcript_type?"transcript":"gene")) continue;
					final String ctg = ctgConverter.apply(entry.getContig());
					if(StringUtils.isBlank(ctg)) continue;
					if(!contigRegex.matcher(ctg).matches()) continue;
					final int pos;
					if(entry.isPostiveStrand()) {
						pos = entry.getStart();
						}
					else
						{
						pos = entry.getEnd();
						}
					final TSS tss = new TSS(ctg, pos,(entry.isPostiveStrand() || !this.use_strand_info?1:-1),this.extend_tss);
					tssMap.put(new Interval(tss),tss);
					}
				}
			if(tssMap.isEmpty()) {
				LOG.error("no TSS found");
				return -1;
				}
			
            final long N_in_reference = count_N_in_reference(this.referenceFile,contigRegex);
			
            
            final long adjusted_reference_length = raw_autosome_length - (N_in_reference);

            if(adjusted_reference_length <=0L) {
				LOG.error("adjusted_reference_length <=0");
				return -1;
            	}
			
			final List<TSSScanner> bams = bamPaths.stream().
					map(B->new TSSScanner(B,referenceFile,samFilter,tssMap,extend_tss)).
					collect(Collectors.toList());
			
            final ExecutorService executorService = Executors.newFixedThreadPool(Math.min(this.nThreads,bams.size()));
			final List<Future<Summary>> future_summaries = executorService.invokeAll(bams);
            executorService.shutdown();
            executorService.awaitTermination(365, TimeUnit.DAYS);
            
			final List<Summary> tss_summaries = new ArrayList<>(bamPaths.size());

            for(Future<Summary> f: future_summaries) {
            	if(!f.isDone() || f.isCancelled() || f.get()==null) {
            		throw new IllegalStateException("future not available ??");
            		}
            	
            	if(f.get().error!=null) {
            		LOG.error("got an error for "+ f.get().sampleName +" "+f.get().error);
            		return -1;
            		}
            	
            	tss_summaries.add(f.get());
            		
            	}
            Collections.sort(tss_summaries,(A,B)->A.sampleName.compareTo(B.sampleName));
			
            
            try(PrefixSuffixWriter pw = new PrefixSuffixWriter(super.openPathOrStdoutAsPrintWriter(output))) {
            	pw.setNoPrefix();
            	pw.println("# This is the output of "+getProgramName());
            	pw.setPrefix("GENERAL\t");
            	pw.println("date\t"+StringUtils.now());
            	pw.println("version\t"+getVersion());
            	pw.println("extend_tss\t"+this.extend_tss);
            	pw.println("fasta\t"+this.referenceFile);
            	pw.println("contig_regex\t"+contigRegex.pattern());
            	pw.println("raw_autosome_length\t"+raw_autosome_length);
            	pw.println("number of N in fasta\t"+N_in_reference);
            	pw.println("adjusted_genome_length\t"+adjusted_reference_length);
            	pw.println("GTF\t"+this.gtfPath);
            	pw.println("number of TSS\t"+tssMap.size());
            	pw.println("treshold_low\t"+treshold_low);
            	pw.println("treshold_high\t"+treshold_high);
            	pw.println("use_strand\t"+this.use_strand_info);
            	pw.println("use_transcript_type\t"+this.use_transcript_type);
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
            	pw.write("sample\ttotal_reads\ttotal_reads_in_tss\tmean_depth\tcount_tss_with_reads");
            	for(int i=0;i< this.extend_tss*2;i++) {
            		pw.write("\t");
            		pw.write(String.valueOf(i-this.extend_tss));
            		}
            	pw.write("\n");
            	for(Summary summary: tss_summaries) {

            		pw.write(summary.sampleName);
            		pw.write("\t");
            		pw.write(String.valueOf(summary.count_reads));
            		pw.write("\t");
            		pw.write(String.valueOf(summary.count_reads_in_tss));
            		pw.write("\t");
            		pw.write(String.valueOf(summary.sum_depth_not_in_tss/(double)(adjusted_reference_length)));
            		pw.write("\t");
            		pw.write(String.valueOf(summary.count_tss_with_reads));

            		for(int i=0;i< summary.array.length;i++) {
                		pw.write("\t");
                		pw.write(String.valueOf((long)summary.array[i]));
                		}
            		pw.write("\n");
            		}
            	
            	// normalisation
            	pw.setNoPrefix();
            	pw.println("# normalized outputs");
            	pw.setPrefix("NORMALIZED\t");
            	final int winsize = Math.max(1, (this.extend_tss*2)/this.num_bins);
            	pw.write("sample");
            	for(int i=0;i< this.extend_tss*2;i+=winsize) {
            		pw.write("\t");
            		pw.write("["+(i-this.extend_tss)+" - "+(i-this.extend_tss+winsize)+"[");
            		}
            	pw.write("\n");
            	for(final Summary summary: tss_summaries) {
            		// normalize on read count in tss
            		if(summary.count_tss_with_reads > 0) {
						for(int i=0;i< summary.array.length;++i) {
							summary.array[i]/=summary.count_tss_with_reads;
							}
						}
            		
            		pw.write(summary.sampleName);
            		final double[] array2 = bin_array(summary.array,winsize);
            		final double genome_mean_cov = summary.sum_depth_not_in_tss/(double)adjusted_reference_length;
            		for(int i=0;i< array2.length;i++) {
            			pw.write("\t");
            			final double score= array2[i]/(genome_mean_cov);
                		pw.write(String.valueOf(score));
                		if(i==array2.length/2) summary.max_ratio =score;
                		}
            		pw.write("\n");
            		}
            	
            	// R PLOT
            	pw.setNoPrefix();
            	pw.println("#\tR script used to plot the normalized results");
            	pw.println("#\t$ awk '$1==\"NORMALIZED\"' output.txt | cut -f 2- > jeter.txt");
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
            	pw.write("abline(h="+treshold_low+",col=\"red\")\n");
            	pw.write("abline(h="+treshold_high+",col=\"green\")\n");
            	pw.write("abline(v=ncol(data_values)/2.0,col=\"gray\")\n");
            	pw.write("axis(1, at = c(0,ncol(data_values)/2.0,ncol(data_values) ),\n");
            	pw.write("     labels = c(\""+(-this.extend_tss)+"\",\"0\",\""+(this.extend_tss)+"\"))\n");
            	pw.write("legend(\"topright\", legend = sample_name, col = sncolors, lty = 1, lwd = 2, cex = 0.8)\n");
            	pw.write("dev.off()\n");

            	
            	pw.setPrefix("HIST\t");
            	pw.println("sample\tindex\tf\tstatus\thistogram");
            	for(final Summary summary: tss_summaries) {
            		final double genome_mean_cov = summary.sum_depth_not_in_tss/(double)(adjusted_reference_length);
            		final double[] array2 = bin_array(summary.array,winsize);
	            	final double max=Arrays.stream(array2).map(V->V/genome_mean_cov).max().orElse(1.0);
					for(int i=0;i< array2.length;i++) {
						final double v = (array2[i]/genome_mean_cov);
						pw.print(summary.sampleName+ "\t"+((i-array2.length/2)*winsize)+"\t"+v+"\t"+  score2status.apply(v) + "\t");
						pw.print(StringUtils.repeat((int)((v/max)*80), '#'));
						pw.println();
						}
	            	}
            	
            	pw.setPrefix("STATUS\t");
            	pw.println("sample\tscore\tstatus");
            	for(final Summary summary: tss_summaries) {
            		pw.print(summary.sampleName);
        			pw.print("\t");
            		pw.print(String.valueOf(summary.max_ratio));
        			pw.print("\t");
        			pw.print(score2status.apply(summary.max_ratio));
        			pw.println();
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
            	for(final Summary summary: tss_summaries) {
            		pw.println("  \"" + summary.sampleName+"\":");
            		final double[] array2 = bin_array(summary.array,winsize);
        			final double genome_mean_cov = summary.sum_depth_not_in_tss/(double)(adjusted_reference_length);

            		for(int i=0;i< this.extend_tss*2;i+=winsize) {
            			pw.println("    \""+ (i-this.extend_tss)+"\": "+array2[i/winsize]/genome_mean_cov);
                		}
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
		new TSSEnrichment().instanceMainWithExit(args);

	}

}
