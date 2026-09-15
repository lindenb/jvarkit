package com.github.lindenb.jvarkit.tools.atacdiff;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.SampleSheet;
import com.github.lindenb.jvarkit.io.SampleSheetFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;

public class AtacSeqDiff extends Launcher {
	private static final Logger LOG = Logger.of(AtacSeqDiff.class);
	
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private Path fasta=null;

	@Parameter(names={"--peaks","--bed"},description="Bed of all peaks",required = true)
	private Path bedForPeaks;
	@Parameter(names={"--mapq","-Q"},description="min mapq")
	private int mapq = 1;
	@Parameter(names={"--fraction","-f"},description="in an interval, take the x fraction of higest coverages")
	private double fraction=0.1;
	@Parameter(names={"--fold"},description="fold")
	private double foldFactor=2.0;
	@Parameter(names={"--output","-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output;

	
	private static class InputBam  {
		Path bam;
		boolean is_case;
		double median_depth;
		SamReader samReader;
		double local_norm_median_cov;
		
	}
	
	private int[] getCoverage(SamReader sr,final Locatable r) {
		final int[] coverage=new int[r.getLengthOnReference()];
		Arrays.fill(coverage, 0);
		try(CloseableIterator<SAMRecord> iter = sr.query(r.getContig(), r.getStart(),r.getEnd(),false)) {
			while(iter.hasNext()) {
				final SAMRecord rec= iter.next();
				if(!SAMRecordDefaultFilter.accept(rec, mapq)) continue;
				for(AlignmentBlock bloc : rec.getAlignmentBlocks()) {
					for(int x=0;x< bloc.getLength();++x) {
						int p = (bloc.getReferenceStart()+x)- r.getStart();
						if(p< 0) continue;
						if(p>= coverage.length) break;
						coverage[p]++;
						}
					}
				}//end iter
			}
		return coverage;
		}
	
	
	@Override
	public int doWork(List<String> args) {
		
		try {
			final List<Interval> intervals;
			try(BedLineReader blr=new BedLineReader(this.bedForPeaks)) {
				intervals = blr.stream().map(B->new Interval(B)).collect(Collectors.toList());
				}
			
			final IntervalTreeMap<Interval> intervalMap = new IntervalTreeMap<Interval>();
			for(Interval r: intervals) {
				intervalMap.put(r, r);
				}
			final String sheetPath = oneAndOnlyOneFile(args);
			final SampleSheet samplesheet ;
			
			samplesheet = new SampleSheetFactory().of(Paths.get(sheetPath));
				
			
			samplesheet.getHeader().assertColumnExists("bam");
			samplesheet.getHeader().assertColumnExists("status");
			
			final SamReaderFactory srf = super.createSamReaderFactory().referenceSequence(this.fasta);

			final List<InputBam> inputBams = new ArrayList<>(samplesheet.size());

			for(FileHeader.RowMap row: samplesheet ) {
				InputBam input = new InputBam();
				input.bam = Paths.get(row.get("bam"));
				final String status = row.get("status");
				if(status.equalsIgnoreCase("case")) {
					input.is_case = true;
					}
				else if(status.equalsIgnoreCase("control") || status.equalsIgnoreCase("ctrl")) {
					input.is_case = false;
					}
				else
					{
					LOG.error("illegal status in "+row);
					return -1;
					}
				input.samReader = srf.open(input.bam);
				inputBams.add(input);
				}
			
			
			for(InputBam inputBam : inputBams) {
				LOG.info("coverage for "+inputBam.bam);
				final DiscreteMedian<Integer> median = new DiscreteMedian<>();
				for(final Interval r: intervals) {
					if(!r.getContig().matches("(chr)?[0-9]+")) continue;
					final int[] coverage=getCoverage(inputBam.samReader,r);
					//higest cov at the end
					Arrays.sort(coverage);
					// take best_idx coverage from the end
					int best_idx= (int)(coverage.length* this.fraction);
					// at least one position
					if(best_idx<=0) best_idx=1;
					for(int i=0;i< best_idx && i< coverage.length;i++) {
						median.add(coverage[(coverage.length-1)-i]);
						}
					}
				inputBam.median_depth = median.getMedian().getAsDouble();
				LOG.info("coverage for "+inputBam.bam+" median cov:"+inputBam.median_depth);
				if(inputBam.median_depth<=0 || Double.isNaN(inputBam.median_depth)) {
					LOG.error("bad median depth for "+inputBam.bam);
					}
				}
			
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.output)) {
				
				pw.print("#chrom");
				pw.print("\t");
				pw.print("start");
				pw.print("\t");
				pw.print("end");
				
				for(InputBam inputBam : inputBams) {
					pw.print("\t");
					pw.print(inputBam.bam);
					}
				pw.print("\t");
				pw.print("control_depth");
				pw.print("\t");
				pw.print("fisher_high");
				pw.print("\t");
				pw.print("fisher_low");
				pw.print("\t");
				pw.print("fisher_best");
				pw.println();
				
				for(final Interval r: intervals) {
					pw.print(r.getContig());
					pw.print("\t");
					pw.print(r.getStart()-1);
					pw.print("\t");
					pw.print(r.getEnd());
					

					for(InputBam inputBam : inputBams) {
						final List<Double> vals = new ArrayList<>();
						final int[] int_coverage=getCoverage(inputBam.samReader,r);
						final double[] norm_coverage =  new double[int_coverage.length];
						for(int i=0;i< norm_coverage.length;++i) {
							norm_coverage[i] = int_coverage[i]/inputBam.median_depth;
							}
						Arrays.sort(norm_coverage);
						// take best_idx coverage from the end
						int best_idx= (int)(norm_coverage.length* this.fraction);
						// at least one position
						if(best_idx<=0) best_idx=1;
						for(int i=0;i< best_idx && i< norm_coverage.length;i++) {
							vals.add(norm_coverage[(norm_coverage.length-1)-i]);
							}
						inputBam.local_norm_median_cov  = new Median().evaluate(vals.stream().mapToDouble(v->v).toArray());
						pw.print("\t");
						pw.print(inputBam.local_norm_median_cov );
						}
					final double control_depth = new Median().evaluate(inputBams.stream()
								.filter(I->I.is_case==false)
								.mapToDouble(I->I.local_norm_median_cov)
								.sorted()
								.toArray()
								);
					pw.print("\t");
					pw.print(control_depth);

					final double high_dp= control_depth*this.foldFactor;
					int count_ctrl_ref = 0;
					int count_ctrl_alt = 0;
					int count_case_ref = 0;
					int count_case_alt = 0;
					for(InputBam inputBam : inputBams) {
						if(inputBam.is_case) {
							if(inputBam.local_norm_median_cov>=high_dp) {
								count_case_alt++;
								}
							else
								{
								count_case_ref++;
								}
							}
						else 
							{
							if(inputBam.local_norm_median_cov>=high_dp) {
								count_ctrl_alt++;
								}
							else
								{
								count_ctrl_ref++;
								}
							}
						}
					double fisher1 = FisherExactTest.compute(
							count_case_ref,count_case_alt,
							count_ctrl_ref,count_ctrl_alt
							).getAsDouble();
					pw.print("\t");
					pw.print(fisher1);

					
					
					final double low_dp= control_depth/this.foldFactor;
					count_ctrl_ref = 0;
					count_ctrl_alt = 0;
					count_case_ref = 0;
					count_case_alt = 0;
					for(InputBam inputBam : inputBams) {
						if(inputBam.is_case) {
							if(inputBam.local_norm_median_cov<=low_dp) {
								count_case_alt++;
								}
							else
								{
								count_case_ref++;
								}
							}
						else 
							{
							if(inputBam.local_norm_median_cov<=low_dp) {
								count_ctrl_alt++;
								}
							else
								{
								count_ctrl_ref++;
								}
							}
						}
					double fisher2 = FisherExactTest.compute(
							count_case_ref,count_case_alt,
							count_ctrl_ref,count_ctrl_alt
							).getAsDouble();
					pw.print("\t");
					pw.print(fisher2);
					pw.print("\t");
					pw.print(Math.min(fisher1, fisher2));
					pw.println();
					}
			pw.flush();
			}
			
			for(InputBam inputBam : inputBams) {
				inputBam.samReader.close();
			}
			
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new AtacSeqDiff().instanceMainWithExit(args);
	}

}
