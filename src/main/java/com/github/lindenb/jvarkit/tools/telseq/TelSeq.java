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
package com.github.lindenb.jvarkit.tools.telseq;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;

/**
BEGIN_DOC 
# Java implementation of https://github.com/zd1/telseq/
 * Ding Z, Mangino M, Aviv A, Spector T, Durbin R; UK10K Consortium. Estimating telomere length from whole genome sequence data. Nucleic Acids Res. 2014;42(9):e75. doi:10.1093/nar/gku181
 * @author lindenb
END_DOC
 */
@Program(name="telseq",
description="Java implementation of https://github.com/zd1/telseq/ (Ding Z & al.).",
keywords={"depth","bam","sam","telomere"},
creationDate="20200907",
modificationDate="20200907",
generate_doc=false
)
public class TelSeq extends Launcher {
	private static final Logger LOG = Logger.build( TelSeq.class).make();
	
	/** default We defined reads as telomeric if they contained k or more TTAGGG repeats, with a default threshold value of k = 7  */
	private static final int TEL_MOTIF_CUTOFF = 7;
	private static final float GC_LOWERBOUND = 0.4f;
	private static final float GC_UPPERBOUND = 0.6f;
	private static final float GC_BINSIZE = 0.02f;

	
	/** "we define s as a fraction of all reads with gas chromatography (GC) composition between 48 and 52%. " */
	private static final float GC_TELOMERIC_LOWERBOUND = 0.48f;
	private static final float GC_TELOMERIC_UPPERBOUND = 0.52f;

	/** convert to Kilobases */
	private final int LENGTH_UNIT = 1_000;
	/** number of telomeres in the human genome 23 chromosomes *2 */
	private static final int TELOMERE_ENDS = 46;
	/** estimation of human telomere length */
	private static final long GENOME_LENGTH_AT_TEL_GC =  332_720_800;
	
	private static class Result {
		//long numMapped = 0L;
		//long duplicate = 0L;
		long numTotal = 0L;
		final List<Long> telcounts = new ArrayList<>();
		final List<Long> gccounts = new ArrayList<>();
		
		void visitKmerTimes(int times) {
			while(telcounts.size()<=times) telcounts.add(0L);
			telcounts.set(times, telcounts.get(times)+1L);
		}
		void visitGCTimes(int gc) {
			while(gccounts.size()<=gc) gccounts.add(0L);
			gccounts.set(gc, gccounts.get(gc)+1L);
		}
	}
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path refPath = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;
	@Parameter(names={"--pattern"},description="pattern")
	private String pattern = "TTAGGG";
	@Parameter(names={"-K","--min-occurence"},description="We defined reads as telomeric if they contained 'k' or more repeats")
	private int tel_k = TEL_MOTIF_CUTOFF;
	@Parameter(names={"--mapped"},description="By default, we query only the unmapped reads. This option uses all the reads.")
	private boolean query_mapped_and_unmapped  = false;

	
	
	
	private double calcGC(final byte[] seq)
	{
	    double num_gc = 0.0;
	    for(int i = 0; i < seq.length; ++i)
	    {
	    	switch(seq[i]) {
	    	case 'c': case 'g':
	    	case 'C': case 'G':  ++num_gc;  break;
	    	default:  break;
	    	}

	    }
	    return num_gc / (double)seq.length;
	}
	
	int countMotif(final byte[] read,final byte[] pattern){

		int motifcount1 = 0;
		int p=0;
		for(;;) {
			if( p+pattern.length > read.length) break;
			int i=0;
			for(i=0;i<pattern.length;i++) {
				if(read[p+i]!=pattern[i]) break;
				}
			
			if(i==pattern.length) {
				p+=pattern.length;
				motifcount1++;
				}
			else
				{
				p++;
				}
		}
		return motifcount1;
	}
	
	
	
	/**
	 * We defined reads as telomeric if they contained 'k' or more TTAGGG repeats
	 * with a default threshold value of k = 7. 
	 * These can then be translated into an estimate of the physical length via a size factor s 
	 * and a constant length c in l = tksc, 
	 * where l is the length estimate, 
	 * tk is the abundance of telomeric reads at threshold k 
	 * and c is a constant for genome length divided by number of telomere ends 46 (23 × 2).
	 * @param results
	 * @return
	 */
	private OptionalDouble calcTelLength(final Result results){

		double acc = 0;
		for(int i = this.tel_k ; i < results.telcounts.size(); ++i){
			acc += results.telcounts.get(i);
		}

		float gc_tel = 0;
		for(int i=0 ; i < results.gccounts.size(); ++i){
			float gc1=GC_LOWERBOUND + GC_BINSIZE*i;
			float gc2=GC_LOWERBOUND +  GC_BINSIZE*(i+1);
			if(gc1 >= GC_TELOMERIC_LOWERBOUND && gc2 <=  GC_TELOMERIC_UPPERBOUND ){
				gc_tel += results.gccounts.get(i);
			}
		}

		if(gc_tel == 0f){
			return OptionalDouble.empty();
		}

		/** This fraction is then converted to a mean telomere length estimate in kilobases by multiplying by the cumulative length of genomic regions with the same GC composition c. */
		
		return OptionalDouble.of((acc/gc_tel)*(float)(GENOME_LENGTH_AT_TEL_GC)/LENGTH_UNIT/TELOMERE_ENDS);
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final byte[] PATTERN=this.pattern.toUpperCase().getBytes();
			final byte[] PATTERN_REVCOMP=AcidNucleics.reverseComplement(this.pattern.toUpperCase()).getBytes();

			
			
			final List<Path> bamPaths = IOUtils.unrollPaths(args);
			if(bamPaths.isEmpty()) {
				LOG.error("no bam was defined");
				return -1;
				}
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.refPath!=null) srf.referenceSequence(this.refPath);
			
			
			try(PrintWriter w = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				w.println("#bam\t"+this.groupBy.name()+"\tcount-reads\ttelomere_size_kb");
				

			
				for(final Path bamPath: bamPaths) {
					try(SamReader r=srf.open(bamPath)) {
						final Map<String, Result> rgid2result = new HashMap<>();
						final SAMFileHeader hdr = r.getFileHeader();
						this.groupBy.getPartitions(hdr.getReadGroups()).forEach(R->{
							rgid2result.putIfAbsent(R, new Result());
						});
						
						try(CloseableIterator<SAMRecord> iter = (this.query_mapped_and_unmapped?r.iterator():r.queryUnmapped())) {
							while(iter.hasNext()) {
								final SAMRecord rec = iter.next();
								if(rec.isSecondaryOrSupplementary()) continue;
								final SAMReadGroupRecord rg = rec.getReadGroup();
								if(rg==null) continue;
								final String rgid = this.groupBy.apply(rg,null);
								if(StringUtils.isBlank(rgid)) continue;
								
								final Result result = rgid2result.get(rgid);
								if(result==null) {
									LOG.error("groupid "+rg.getId() +" is not defined in the header ? " +String.join(", ", rgid2result.keySet()));
									continue;
								}
								result.numTotal++ ;
								//if(rec.getReadUnmappedFlag()) result.numMapped++;
								//if(rec.getDuplicateReadFlag()) result.duplicate++;
								final byte[] bases=rec.getReadBases();
								if(bases==null || bases.length==0 || bases==SAMRecord.NULL_SEQUENCE ) continue;
								final double gc = calcGC(bases);
								
								if(gc >= GC_LOWERBOUND && gc <= GC_UPPERBOUND){
									  // get index for GC bin.
									  final int idx =(int) Math.floor((gc- GC_LOWERBOUND)/GC_BINSIZE);
									  result.visitGCTimes(idx);
									}
								
								final int ptn_count = Math.max(countMotif(bases,PATTERN) , countMotif(bases,PATTERN_REVCOMP));
								
								result.visitKmerTimes(ptn_count);
								}
						}
			
					for(final String rgid: rgid2result.keySet()) {
						final Result result = rgid2result.get(rgid);
						w.print(bamPath);
						w.print("\t");
						w.print(rgid);
						w.print("\t");
						w.print(result.numTotal);
						w.print("\t");
						final OptionalDouble size = this.calcTelLength(result);
						if(size.isPresent()) {
							w.print(size.getAsDouble());
							}
						else
							{
							w.print("N/A");
							}
						w.println();
						}
					w.flush();
					}
				}
			}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new TelSeq().instanceMainWithExit(args);
	}
}
