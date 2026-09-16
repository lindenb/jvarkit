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
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.SampleSheet;
import com.github.lindenb.jvarkit.io.SampleSheetFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.wig.BigWigReader;

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
	@Parameter(names={"--fraction","-f"},description="in an interval, take the x fraction of higest coverages")
	private double fraction=0.1;
	@Parameter(names={"--fold"},description="fold")
	private double foldFactor=2.0;
	@Parameter(names={"--output","-o"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output;
	@Parameter(names={"--min-dp","-D"},description="min norm DEPTH")
	private double min_depth=0.2;

	
	private static class InputBigWig {
		String path;
		String sampleName;
		boolean is_case;
		BigWigReader bigWigReader;
		double local_norm_median_cov;
		
	}
	
	private double[] getCoverage(final BigWigReader reader,final Locatable loc) {
		final double[] coverage=new double[loc.getLengthOnReference()];
		Arrays.fill(coverage, 0.0);
		try(CloseableIterator<BigWigReader.WigItem> iter = reader.query(loc)) {
			while(iter.hasNext()) {
				final BigWigReader.WigItem rec= iter.next();
				for(int x=0;x< rec.getLengthOnReference();++x) {
						final int idx = (rec.getStart()+x)-loc.getStart();
						if(idx<0) break;
						if(idx>=coverage.length) continue;
						coverage[idx]= rec.getValue();
						}
					}
			}
		return coverage;
		}
	
	private double[] bestOf(final double[] array,double fraction) {
		Arrays.sort(array);
		int n= Math.max((int)(array.length*fraction),1);
		final double[] a2 = new double[n];
		System.arraycopy(array, array.length-n, a2, 0, n);
		return a2;
		}
	
	@Override
	public int doWork(final List<String> args) {
		
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
				
			
			samplesheet.getHeader().assertColumnExists("bigWig");
			samplesheet.getHeader().assertColumnExists("status");
			

			final List<InputBigWig> inputBigWigs = new ArrayList<>(samplesheet.size());

			for(FileHeader.RowMap row: samplesheet ) {
				final InputBigWig input = new InputBigWig();
				input.path = row.get("bigWig");
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
				input.bigWigReader = new BigWigReader(input.path);
				if(row.containsKey("sample")) {
					input.sampleName = row.get("sample");
					}
				else
					{
					input.sampleName = IOUtils.getFilenameWithoutCommonSuffixes(Paths.get(input.path));
					}
				inputBigWigs.add(input);
				}
			
			
			
			
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.output)) {
				
				pw.print("#chrom");
				pw.print("\t");
				pw.print("start");
				pw.print("\t");
				pw.print("end");
				
				for(InputBigWig inputWig : inputBigWigs) {
					pw.print("\t");
					pw.print(inputWig.sampleName);
					pw.print(inputWig.is_case?"[case]":"[control]");
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

					for(InputBigWig inputWig : inputBigWigs) {
						final double[] coverage=  bestOf( getCoverage(inputWig.bigWigReader,r), this.fraction);
						
						inputWig.local_norm_median_cov  = new Median().evaluate(coverage);
						}
					
					if(inputBigWigs.stream().allMatch(B->B.local_norm_median_cov< AtacSeqDiff.this.min_depth)) {
						continue;
						}
					
					final double control_depth = new Median().evaluate(inputBigWigs.stream()
								.filter(I->I.is_case==false)
								.mapToDouble(I->I.local_norm_median_cov)
								.sorted()
								.toArray()
								);
				
			
						
						pw.print(r.getContig());
						pw.print("\t");
						pw.print(r.getStart()-1);
						pw.print("\t");
						pw.print(r.getEnd());

						for(InputBigWig inputWig : inputBigWigs) {
							pw.print("\t");
							pw.print((float)inputWig.local_norm_median_cov);
							}
						pw.print("\t");
						pw.print(control_depth);

						
						final double high_dp= control_depth*this.foldFactor;
						int count_ctrl_ref = 0;
						int count_ctrl_alt = 0;
						int count_case_ref = 0;
						int count_case_alt = 0;
						for(InputBigWig inputWig : inputBigWigs) {
							if(inputWig.is_case) {
								if(inputWig.local_norm_median_cov>=high_dp) {
									count_case_alt++;
									}
								else
									{
									count_case_ref++;
									}
								}
							else 
								{
								if(inputWig.local_norm_median_cov>=high_dp) {
									count_ctrl_alt++;
									}
								else
									{
									count_ctrl_ref++;
									}
								}
							}
						final double fisher1 = FisherExactTest.compute(
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
						for(InputBigWig inputWig : inputBigWigs) {
							if(inputWig.is_case) {
								if(inputWig.local_norm_median_cov<=low_dp) {
									count_case_alt++;
									}
								else
									{
									count_case_ref++;
									}
								}
							else 
								{
								if(inputWig.local_norm_median_cov<=low_dp) {
									count_ctrl_alt++;
									}
								else
									{
									count_ctrl_ref++;
									}
								}
							}
						final double fisher2 = FisherExactTest.compute(
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
			
			for(InputBigWig inputBam : inputBigWigs) {
				inputBam.bigWigReader.close();
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
