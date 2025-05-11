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
package com.github.lindenb.jvarkit.tools.phased;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Example

```
$ find src/test -type f -name "*.bam" > jeter.list
$ java -jar dist/bam2mnv.jar --input jeter.list src/test/resources/rotavirus_rf.vcf.gz --distance 3000

#CHROM1	POS1	REF1	ALT1	CHROM2	POS2	REF2	ALT2	distance	S1	S2	S3	S4	S5
RF03	1221	C	G	RF03	1242	C	A	22	ref	v1	v1	v2	ref
RF03	1221	C	G	RF03	1688	T	G	468	ref	v1	v1	ref	v2
RF03	1221	C	G	RF03	1708	G	T	488	ref	v1	v1	ref	v2
RF03	1242	C	A	RF03	1688	T	G	447	ref	ref	ref	v1	v2
RF03	1242	C	A	RF03	1708	G	T	467	ref	ref	ref	v1	v2
RF03	1688	T	G	RF03	1708	G	T	21	ref	ref	ref	ref	cis
RF03	1688	T	G	RF03	2150	T	A	463	ref	ref	ref	ref	v1
RF03	1688	T	G	RF03	2201	G	C	514	ref	v2	v2	.	.
RF03	1708	G	T	RF03	2150	T	A	443	.	ref	ref	ref	v1
RF03	1708	G	T	RF03	2201	G	C	494	ref	v2	v2	.	.
RF03	2150	T	A	RF03	2201	G	C	52	ref	v2	v2	ref	.
RF03	2150	T	A	RF03	2573	A	G	424	v1	v2	v2	ref	.
RF03	2201	G	C	RF03	2573	A	G	373	.	cis	cis	.	.
RF04	1900	A	C	RF04	1920	A	T	21	v2	ref	ref	ref	v1
RF05	41	T	C	RF05	499	A	T	459	ref	v2	v2	v1	ref
RF05	879	C	A	RF05	1297	T	G	419	ref	cis	cis	ref	ref
RF05	1297	T	G	RF05	1339	A	C	43	ref	cis	cis	ref	.
RF06	517	C	A	RF06	543	G	C	27	.	v2	v2	.	ref
RF06	668	A	G	RF06	695	T	C	28	ref	ref	ref	.	v1
RF07	225	C	A	RF07	684	T	G	460	ref	.	.	v2	.
RF08	926	A	C	RF08	992	G	C	67	ref	v2	v2	.	v1
RF09	294	T	C	RF09	317	C	A	24	.	ref	ref	ref	v2
```

END_DOC
*/
@Program(name="bam2mnv",
description="MNV haplotypes",
keywords={"vcf","phased","genotypes","bam"},
creationDate="20211208",
modificationDate="20211208",
generate_doc=false,
jvarkit_amalgamion = true,
menu="BAM Manipulation"
)
public class BamToMNV extends Launcher {
	private static final Logger LOG=Logger.of(BamToMNV.class);

	
	private enum Phase {
		cis,trans,ambigous,ref,v1,v2
		}

	private class Mnv {
		private final int idx1;
		private final int idx2;
		private Map<String, Phase> sample2phase = new HashMap<>();
		Mnv(int idx1,int idx2) {
			this.idx1 = idx1;
			this.idx2 = idx2;
			}
		VariantContext get(int i) {
			switch(i) {
				case 0: case 1: return BamToMNV.this.all_variants.get(i==0?idx1:idx2);
				default: throw new IllegalArgumentException();
				}
			}
		
		Phase getPhase(final List<SAMRecord> records) {
			int count_alt_alt=0;
			int count_ref_ref=0;
			int count_alt_ref=0;
			int count_ref_alt=0;
			
			final EqualIterator<SAMRecord> equal_range = new EqualIterator<>(records.iterator(),(A,B)->A.getReadName().compareTo(B.getReadName()));
			while(equal_range.hasNext()) {
				final List<SAMRecord> L = equal_range.next();
				int[] count_REF = new int[]{0,0};
				int[] count_ALT = new int[]{0,0};

				for(int side=0;side<2;++side) {
					final VariantContext ctx = get(side);
					for(SAMRecord rec:L) {
						final Allele allele = findAllele(ctx, rec).orElse(null);
						if(allele==null) continue;
						if(allele.equals(ctx.getAlleles().get(0))) count_REF[side]++;
						else if(allele.equals(ctx.getAlleles().get(1))) count_ALT[side]++;
						}
					}
				if(count_ALT[0]>0 && count_ALT[1]>0 && count_REF[0]==0 && count_REF[1]==0) {
					count_alt_alt++;
					}
				else if(count_ALT[0]==0 && count_ALT[1]==0 && count_REF[0]>0 && count_REF[1]>0) {
					count_ref_ref++;
					}
				else if(count_ALT[0]>0 && count_REF[0]==0 && count_ALT[1]==0 && count_REF[1]>0) {
					count_alt_ref++;
					}
				else if(count_ALT[0]==0 && count_REF[0]>0 && count_ALT[1]>0 && count_REF[1]==0) {
					count_ref_alt++;
					}
				}
			equal_range.close();
			
			if(count_alt_alt>0 && count_alt_ref==0 && count_ref_alt==0 && count_ref_ref==0) {
				return Phase.cis;
				}
			else if(count_alt_alt==0 && count_alt_ref>0 && count_ref_alt>0 && count_ref_ref==0) {
				return Phase.trans;
				}
			else if(count_alt_alt==0 && count_alt_ref==0 && count_ref_alt==0 && count_ref_ref>0) {
				return Phase.ref;
				}
			else if(count_alt_alt==0 && count_alt_ref>0 && count_ref_alt==0 && count_ref_ref==0) {
				return Phase.v1;
				}
			else if(count_alt_alt==0 && count_alt_ref==0 && count_ref_alt>0 && count_ref_ref==0) {
				return Phase.v2;
				}
			return Phase.ambigous;
			}
		
		}

	

	

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	protected Path faidx=null;
	@Parameter(names={"-I","--input","--bam"},description="add this indexed bam")
	protected List<String> input_bams = new ArrayList<>();
	@Parameter(names={"--distance"},description="max distance distance between two snp.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int min_distance_mnv = 100;
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int minmapq=1;
	@Parameter(names={"--pedigree","--ped"},description=PedigreeParser.OPT_DESC)
	private Path pedigreePath = null;
	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	final List<VariantContext> all_variants = new ArrayList<>();


	
	private byte alleleToChar(final Allele a) {
		if(a==null || a.isSymbolic() || a.isNoCall() || a.length()!=1) throw new IllegalArgumentException("length!=1 for allele "+a);
		return a.getBases()[0];
	}
	
	private Optional<Allele> findAllele(final VariantContext ctx,final SAMRecord rec) {
		if(!ctx.overlaps(rec)) return Optional.empty();
		final byte[] bases = rec.getReadBases();
		final byte ref = alleleToChar(ctx.getReference());
		final byte alt = alleleToChar(ctx.getAlleles().get(1));
		for(final AlignmentBlock ab:rec.getAlignmentBlocks()) {
			final int readPos1 = ab.getReadStart();
			final int refPos1 = ab.getReferenceStart();
			for(int x=0;x< ab.getLength();++x) {
				final int refPos1_x = refPos1+x;
				if(refPos1_x<ctx.getStart()) continue;
				if(refPos1_x>ctx.getStart()) break;				
				final byte readBase = bases [ (readPos1-1) + x ];
				if(readBase==ref) return Optional.of(ctx.getReference());
				if(readBase==alt) return Optional.of(ctx.getAlleles().get(1));
				return Optional.empty();
				}
			}
		return Optional.empty();
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		try {
			final List<Path> bams = IOUtils.unrollPaths(this.input_bams);
			if(bams.isEmpty()) {
				LOG.error("No bam was defined");
				return -1;
				}
			final Pedigree pedigree;
			if(this.pedigreePath!=null) {
				pedigree = new PedigreeParser().parse(this.pedigreePath);
				pedigree.checkUniqIds();
			} else {
				pedigree=null;
			}
			
			try(VCFReader reader=VCFReaderFactory.makeDefault().open(oneAndOnlyOneFile(args),false)) {
				final VCFHeader header = reader.getHeader();
				final OrderChecker<VariantContext> order = new OrderChecker<>(header.getSequenceDictionary(),false);
				try(CloseableIterator<VariantContext> r = reader.iterator()) {
					this.all_variants.addAll(r.stream().
						filter(V->V.isBiallelic() && V.isSNP()).
						map(V->new VariantContextBuilder(V).noGenotypes().make()).
						map(order).
						collect(Collectors.toList()));
					}
				}
			final List<Mnv> all_mnvs = new ArrayList<>();
			for(int i=0; i+1 < this.all_variants.size();i++) {
				final VariantContext v1 = this.all_variants.get(i);
				for(int j=i+1; j < this.all_variants.size();j++) {
					final VariantContext v2 = this.all_variants.get(j);
					if(!v1.withinDistanceOf(v2, min_distance_mnv)) break;
					if(v1.getStart()==v2.getStart()) continue;
					all_mnvs.add(new Mnv(i,j));
					}
				}
			final Set<String> samples = new TreeSet<>();
			final SamReaderFactory srf = super.createSamReaderFactory().referenceSequence(this.faidx);
			for(final Path bam:bams) {
				LOG.info(String.valueOf(bam));
				try(SamReader sr = srf.open(bam)) {
					final SAMFileHeader header = sr.getFileHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					final String sample = header.getReadGroups().stream().
							map(R->R.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(bam));
					if(samples.contains(sample)) {
						LOG.error("duplicate sample "+sample);
						return -1;
						}
					samples.add(sample);
					if(pedigree!=null && pedigree.getSampleById(sample)==null) {
						LOG.warn("sample "+sample+" from "+bam+ " is not in pedigree.");
					}

					if(all_mnvs.isEmpty()) continue;
					
					final QueryInterval[] intervals = QueryInterval.optimizeIntervals(this.all_variants.stream().
						map(V->new QueryInterval(dict.getSequenceIndex(V.getContig()),V.getStart(),V.getEnd())).
						toArray(X->new QueryInterval[X]));
					final List<SAMRecord> sam_reads = new ArrayList<>();
					try(CloseableIterator<SAMRecord> iter = sr.query(intervals, false)) {
						while(iter.hasNext()) {
							final SAMRecord record = iter.next();
							if(!SAMRecordDefaultFilter.accept(record,this.minmapq)) continue;
							if(record.getReadBases()==SAMRecord.NULL_SEQUENCE) continue;
							sam_reads.add(record);
							}
						}
					//sort on query name
					Collections.sort(sam_reads,(A,B)->A.getReadName().compareTo(B.getReadName()));
					
					for(final Mnv mnv:all_mnvs) {
						final Phase phase = mnv.getPhase(sam_reads);
						if(phase.equals(Phase.ambigous)) continue;
						mnv.sample2phase.put(sample, phase);
						}
					}
				}
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				pw.print("#CHROM1\tPOS1\tREF1\tALT1");
				pw.print("\tCHROM2\tPOS2\tREF2\tALT2");
				pw.print("\tdistance");
				for(final String sn:samples) pw.print("\t"+sn);
				if(pedigree!=null) {
					pw.print("\tcase-cis\tcase-trans\tcontrol-cis\tcontrol-trans\tfisher");
					}
				pw.println();

				for(final Mnv mnv:all_mnvs) {
					if(mnv.sample2phase.values().stream().allMatch(V->V.equals(Phase.ambigous) || V.equals(Phase.ref))) continue;
					for(int side=0;side<2;++side) {
						final VariantContext ctx= mnv.get(side);
						if(side>0) pw.print("\t");
						pw.print(ctx.getContig());
						pw.print("\t");
						pw.print(ctx.getStart());
						pw.print("\t");
						pw.print(ctx.getReference().getDisplayString());
						pw.print("\t");
						pw.print(ctx.getAlleles().get(1).getDisplayString());
						}
					pw.print("\t");
					pw.print(CoordMath.getLength(mnv.get(0).getStart(), mnv.get(1).getEnd()));
					int case_cis=0;
					int case_trans=0;
					int ctrl_cis=0;
					int ctrl_trans=0;
					for(final String sn:samples) {
						pw.print("\t");
						final Phase phase = mnv.sample2phase.get(sn);
						if(phase==null) {
							pw.print(".");
							continue;
							}
						pw.print(phase.name());
						if(pedigree!=null) {
							final Sample sample = pedigree.getSampleById(sn);
							if(sample==null) {
								//nothing
								}
							else if(sample.isAffected()) {
								if(phase.equals(Phase.cis)) {
									case_cis++;
									}
								else if(phase.equals(Phase.trans)) {
									case_trans++;
									}
								}
							else if(sample.isUnaffected()) {
								if(phase.equals(Phase.cis)) {
									ctrl_cis++;
									}
								else if(phase.equals(Phase.trans)) {
									ctrl_trans++;
									}
								}
							}
						}
					if(pedigree!=null) {
						pw.print("\t");
						pw.print(case_cis);
						pw.print("\t");
						pw.print(case_trans);
						pw.print("\t");
						pw.print(ctrl_cis);
						pw.print("\t");
						pw.print(ctrl_trans);
						pw.print("\t");
						final FisherExactTest fisher=FisherExactTest.compute(case_cis, case_trans, ctrl_cis, ctrl_trans);
						pw.print(fisher.getAsDouble());
						}
					pw.println();
					}
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
		new BamToMNV().instanceMainWithExit(args);

	}

}
