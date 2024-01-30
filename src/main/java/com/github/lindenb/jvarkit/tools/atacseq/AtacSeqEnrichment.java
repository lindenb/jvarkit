package com.github.lindenb.jvarkit.tools.atacseq;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.log.Logger;

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

public class AtacSeqEnrichment extends Launcher {
	private static final Logger LOG = Logger.build(AtacSeqEnrichment.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	protected Path referenceFile = null;
    @Parameter(names={"--extend","-x"},description="extend tss site by 'x' bases." + DistanceParser.OPT_DESCRIPTION, splitter = NoSplitter.class, converter  = DistanceParser.StringConverter.class)
    private int extend_tss = 2_000;
    @Parameter(names={"--gtf","-gtf"},description="GFF file to be used as a source of TSS start sites")
    private Path gtfPath = null;


    @Parameter(names={"--threads"},description="number of threads")
    private int nThreads = 1;

	private static class TSS implements Locatable {
		private final String contig;
		private final int position;
		private final int extend;
		TSS(final String c,int pos,final int extend) {
			this.contig = c;
			this.position = pos;
			this.extend = extend;
			}
		@Override
		public String getContig() {
			return contig;
			}
		public int getPosition() {
			return position;
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
			return Math.min(
					Math.abs(rec.getStart()-this.getPosition()),
					Math.abs(rec.getEnd()-this.getPosition())
					);
			}
		}
	
	private static class Summary {
		Throwable error=null;
		String sampleName;
		long total = 0L;
		int[] array;
		}
	
	private static class BamScanner implements Callable<Summary> {
		private final Path bamPath;
		private final Path fastaPath;
		private final Predicate<SAMRecord> samFilter;
		private final IntervalTreeMap<TSS> intervalTreeMap;
		private final int window_size;
		BamScanner(
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
		@Override
		public Summary call() throws Exception {
			final Summary summary = new Summary();
			summary.array = new int[window_size*2];
			Arrays.fill(summary.array, 0);
			StopWatch stopWatch = new StopWatch();
			stopWatch.start();
			try {
				final SamReaderFactory srf = SamReaderFactory.makeDefault();
				srf.validationStringency(ValidationStringency.LENIENT);
				srf.referenceSequence(this.fastaPath);
				try (SamReader sr = srf.open(this.bamPath)){
					final SAMFileHeader header = sr.getFileHeader();
					summary.sampleName = header.getReadGroups().
							stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(this.bamPath));
					
					try(CloseableIterator<SAMRecord> iter= sr.iterator()) {
						while(iter.hasNext()) {
							final SAMRecord record = iter.next();
							if(!this.samFilter.test(record)) continue;
							summary.total++;
							TSS closest = null;
							for(TSS tss:this.intervalTreeMap.getOverlapping(record)) {
								if(closest==null || tss.distanceTo(record) < closest.distanceTo(record)) {
									closest = tss;
									}
								}
							if(closest!=null) {
								final int mid_record = record.getStart() + record.getReadLength()/2;
								final int distance_to_tss = mid_record - closest.getPosition() ; /* ok negative if read is before TSS */
								final int array_index = summary.array.length/2  /* middle of array */ + distance_to_tss;
								if(array_index >=0 && array_index<summary.array.length ) {
									summary.array[array_index]++;
									}
								}
							}
						}
					}
				
				}
			catch(Throwable error) {
				summary.error = error;
				}
			stopWatch.start();
			LOG.info("Generation : That took:" + com.github.lindenb.jvarkit.lang.StringUtils.niceDuration(stopWatch.getElapsedTime()) );
			return summary;
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final IntervalTreeMap<TSS> tssMap = new IntervalTreeMap<>();
			final Predicate<SAMRecord> samFilter = REC->SAMRecordDefaultFilter.accept(REC);
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			if(bamPaths.isEmpty()) {
				LOG.error("input missing");
				return -1;
				}
			
			if(gtfPath==null) {
				LOG.error("GTF file missing");
				return -1;
				}
			
			try(BufferedReader br =  IOUtils.openPathForBufferedReading(this.gtfPath)) {
				String line;
				final GTFCodec gtfCodec = new GTFCodec();
				final Function<String,String> ctgConverter;
				if(referenceFile!=null) {
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(referenceFile);
					ctgConverter = ContigNameConverter.fromOneDictionary(dict);
					}
				else
					{
					ctgConverter = S->S;
					}
				while((line=br.readLine())!=null) {
					if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
					final GTFLine entry= gtfCodec.decode(line);
					if(entry==null || !entry.getType().equals("transcript")) continue;
					final String ctg = ctgConverter.apply(entry.getContig());
					if(StringUtils.isBlank(ctg)) continue;
					int pos;
					if(entry.isPostiveStrand()) {
						pos = entry.getStart();
						}
					else
						{
						pos = entry.getEnd();
						}
					final TSS tss = new TSS(ctg, pos,this.extend_tss);
					tssMap.put(new Interval(tss),tss);
					}
				}
			if(tssMap.isEmpty()) {
				LOG.error("no TSS found");
				return -1;
				}
			final List<BamScanner> bams = bamPaths.stream().
					map(B->new BamScanner(B,referenceFile,samFilter,tssMap,extend_tss)).
					collect(Collectors.toList());
			
            final ExecutorService executorService = Executors.newFixedThreadPool(this.nThreads);
			final List<Future<Summary>> summaries = executorService.invokeAll(bams);
            executorService.shutdown();
            executorService.awaitTermination(365, TimeUnit.DAYS);
            
            

			
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new AtacSeqEnrichment().instanceMainWithExit(args);

	}

}
