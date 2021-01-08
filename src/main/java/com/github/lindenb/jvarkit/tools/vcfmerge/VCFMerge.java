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
package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFRecordCodec;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

The motivation for this is to merge a large number of VCF files without opening a bunch of temporary files.

For a regular normal number of files you should use  GATK combineVariants or bcftools merge
 
## Example


```bash
$ java -jar dist/vcfmerge.jar -hr src/test/resources/S*.vcf.gz | more
[INFO][VCFMerge]merging...5 vcfs
##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	970	.	A	C	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	0/0	1/1
RF01	3246	.	A	G	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	0/0	1/1
RF02	578	.	G	A	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	1/1	0/0
RF02	877	.	T	A	.	.	AC=1;AF=0.100;AN=10	GT	0/1	0/0	0/0	0/0	0/0
RF02	1962	.	TACA	TA	.	.	AC=1;AF=0.100;AN=10	GT	0/1	0/0	0/0	0/0	0/0
RF02	2332	.	AT	A	.	.	AC=1;AF=0.100;AN=10	GT	0/0	0/0	0/0	0/1	0/0
RF02	2662	.	G	C	.	.	AC=2;AF=0.200;AN=10	GT	0/0	0/0	0/0	0/0	1/1

```

END_DOC
 */
@Program(name="vcfmerge",
	description="Merge a large number of VCF Files",
	keywords={"vcf","sort","merge"},
	creationDate="20130916",
	modificationDate="20201209"
	)
public class VCFMerge
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFMerge.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-homref","--homref","-hr"},description="Use HomRef 0/0 for unknown variant")
	private boolean useHomRefForUnknown = false;
	@Parameter(names={"-region","--region","-r"},description="Merge in that region: " + IntervalParserFactory.OPT_DESC )
	private String regionStr = "";
	@Parameter(names={"--ploidy"},description="Ploidy")
	private int ploidy=2;
	@Parameter(names={"--fields"},description="print the following INFO/FORMAT Fields.")
	private String fieldsString= String.join(",",VCFConstants.ALLELE_COUNT_KEY,
			VCFConstants.ALLELE_NUMBER_KEY,
			VCFConstants.ALLELE_FREQUENCY_KEY,
			VCFConstants.DEPTH_KEY,
			VCFConstants.GENOTYPE_QUALITY_KEY,
			VCFConstants.GENOTYPE_ALLELE_DEPTHS,
			VCFConstants.GENOTYPE_PL_KEY
			);;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

	
	@Override
	public int doWork(final List<String> args) {
		final List<Path> userVcfFiles=new ArrayList<Path>();
		final Set<String> genotypeSampleNames=new TreeSet<String>();
		SAMSequenceDictionary dict=null;
		VariantContextWriter w=null;
		SortingCollection<VariantContext> array = null;
		CloseableIterator<VariantContext> iter=null;
		try
			{
			userVcfFiles.addAll(IOUtils.unrollPaths(args));
			final Set<String> fieldsSet = Arrays.stream(this.fieldsString.split("[,; ]")).collect(Collectors.toSet());
			final boolean with_ac = fieldsSet.contains(VCFConstants.ALLELE_COUNT_KEY);
			final boolean with_an = fieldsSet.contains(VCFConstants.ALLELE_NUMBER_KEY);
			final boolean with_af = fieldsSet.contains(VCFConstants.ALLELE_FREQUENCY_KEY);
			final boolean with_dp = fieldsSet.contains(VCFConstants.DEPTH_KEY);
			final boolean with_gq = fieldsSet.contains(VCFConstants.GENOTYPE_QUALITY_KEY);
			final boolean with_ad = fieldsSet.contains(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
			final boolean with_pl = fieldsSet.contains(VCFConstants.GENOTYPE_PL_KEY);
			
			
			
			if(userVcfFiles.isEmpty())
				{
				LOG.error("No input");
				return -1;
				}
		
			final boolean requireIndex = !StringUtils.isBlank(this.regionStr);

			for(final Path vcfFile:userVcfFiles) {
				try(VCFReader in= VCFReaderFactory.makeDefault().open(vcfFile,requireIndex)){
					final VCFHeader header= in.getHeader();
					for(final String sn:header.getSampleNamesInOrder()) {
						if(genotypeSampleNames.contains(sn)) {
							LOG.error("duplicate sample name "+sn);
							return -1;
							}
						genotypeSampleNames.add(sn);
						}
					final SAMSequenceDictionary dict1= SequenceDictionaryUtils.extractRequired(header);
					if(dict==null) {
						dict=dict1;
					} else  {
						SequenceUtil.assertSequenceDictionariesEqual(dict, dict1);
					}
				}
			}
			
			if(dict==null) {
				LOG.error("No sequence dictionary defined");
				return -1;
			}
			
			final SAMSequenceDictionary finalDict = dict;
			final Function<String,Integer> contig2tid=C->{
				final int tid = finalDict.getSequenceIndex(C);
				if(tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(C, finalDict);
				return tid;
				};
			
			final Comparator<String> compareContigs = (C1,C2)->{
				if(C1.equals(C2)) return 0;
				return contig2tid.apply(C1) - contig2tid.apply(C2);
				};
			
			final Comparator<VariantContext> compareChromPos = (V1,V2)->{
				int i = compareContigs.compare(V1.getContig(),V2.getContig());
				if( i!=0 ) return i;
				return V1.getStart() - V2.getStart();
				};	
			final Comparator<VariantContext> compareChromPosRef = (V1,V2)->{
				int i = compareChromPos.compare(V1,V2);
				if( i!=0 ) return i;
				return V1.getReference().compareTo(V2.getReference());
				};

			
			final SimpleInterval rgn;
			final Predicate<VariantContext> accept;

			if(!StringUtil.isBlank(VCFMerge.this.regionStr)) {
				rgn = IntervalParserFactory.newInstance().
						dictionary(dict).
						enableWholeContig().
						make().
						apply(VCFMerge.this.regionStr).
						orElseThrow(IntervalParserFactory.exception(VCFMerge.this.regionStr));
				accept = (CTX)->{
					return rgn.overlaps(CTX);
					};
				}
			else
				{
				accept = (VOL) -> true;	
				rgn = null;
				}
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_KEY);
			
	
	
			if(with_dp) VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.DEPTH_KEY);
			if(with_gq) VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_QUALITY_KEY);
			if(with_pl) VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_PL_KEY);
			if(with_ad) VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_ALLELE_DEPTHS);			
			
			if(with_ac) VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.ALLELE_COUNT_KEY);
			if(with_an) VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.ALLELE_NUMBER_KEY);
			if(with_af) VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.ALLELE_FREQUENCY_KEY);
			if(with_dp) VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.DEPTH_KEY);
			
			final 	VCFHeader	mergedHeader=new VCFHeader(
					metaData,
					genotypeSampleNames
					);
			mergedHeader.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, mergedHeader);
			
			array= SortingCollection.newInstance(
					VariantContext.class,
					new VCFRecordCodec(mergedHeader),
					compareChromPosRef,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			array.setDestructiveIteration(true);
			
			for(final Path vcfFile:userVcfFiles) {
				try(VCFReader in= VCFReaderFactory.makeDefault().open(vcfFile,requireIndex)){
					try(CloseableIterator<VariantContext> lit=(in.isQueryable() && rgn!=null ?in.query(rgn):in.iterator())) {
						while(lit.hasNext())
							{					
							final VariantContext  ctx = lit.next();
							if(!accept.test(ctx)) continue;
							array.add(new VariantContextBuilder(ctx).
									unfiltered().
									genotypes(ctx.getGenotypes().stream().filter(G->G.isCalled()).map(G->{
										final GenotypeBuilder gb= new GenotypeBuilder(G);
										gb.noAttributes();
										return gb.make();
										}).collect(Collectors.toList())).
									rmAttributes(new ArrayList<>(ctx.getAttributes().keySet())).make());
							}
						}
		
					}
				}
			array.doneAdding();
			LOG.info("merging..."+userVcfFiles.size()+" vcfs");
	
			
	
			//create the context writer
			w= this.writingVariantsDelegate.open(outputFile);
			w.writeHeader(mergedHeader);
			iter= array.iterator();
			EqualRangeIterator<VariantContext> eqiter = new EqualRangeIterator<>(iter, compareChromPosRef);
			while(eqiter.hasNext())
				{
				final List<VariantContext> row = eqiter.next();
				final VariantContext first = row.get(0);
				final List<Allele> alleles = new ArrayList<>();
				alleles.add(first.getReference());
				alleles.addAll(row.stream().flatMap(VC->VC.getAlternateAlleles().stream()).
						filter(A->!A.isNoCall()).
						collect(Collectors.toSet()));
				final VariantContextBuilder vcb=new VariantContextBuilder(
						null,
						first.getContig(),
						first.getStart(),
						first.getEnd(),
						alleles
						);
				final String id = row.stream().filter(V->V.hasID()).map(ST->ST.getID()).findFirst().orElse(null);
				if(!StringUtils.isBlank(id)) {
					vcb.id(id);
				}
				
				final Map<String,Genotype> sample2genotypes = new HashMap<>(genotypeSampleNames.size());
				final Set<String> remainingSamples=new HashSet<String>(genotypeSampleNames);
				int an=0;
				int dp=-1;
				final Counter<Allele> ac = new Counter<>();
				
				for(final VariantContext ctx:row)
					{
					for(final Genotype gt:ctx.getGenotypes()) {
						if(gt.isNoCall()) continue;
						for(Allele a: gt.getAlleles()) {
							ac.incr(a);
							an++;
							}
						final GenotypeBuilder gb=new GenotypeBuilder(gt.getSampleName(), gt.getAlleles());
						if(with_dp && gt.hasDP()) {
							if(dp<0) dp=0;
							dp+= gt.getDP();
							gb.DP(gt.getDP());
						}
						if(with_gq && gt.hasGQ()) gb.GQ(gt.getGQ());
						if(with_pl && gt.hasPL()) gb.PL(gt.getPL());
						if(with_ad && gt.hasAD()) {
							final int src_ad[]= gt.getAD();
							final int dest_ad[]= new int[alleles.size()];
							Arrays.fill(dest_ad, 0);
							for(int i=0;i< src_ad.length && i< ctx.getAlleles().size();i++) {
								final Allele a1 = ctx.getAlleles().get(i);
								final int dest_idx = alleles.indexOf(a1);
								if(dest_idx>=0 && dest_idx < dest_ad.length) {
									dest_ad[dest_idx] = src_ad[i];
								}
							gb.AD(dest_ad);
							}
						}
						sample2genotypes.put(gt.getSampleName(), gb.make());
						}
					
					}
				remainingSamples.removeAll(sample2genotypes.keySet());
				for(String sampleName:remainingSamples)
					{
					final Genotype gt;
					if(this.useHomRefForUnknown)
						{
						final List<Allele> list = new ArrayList<>(ploidy);
						for(int i=0;i<ploidy;i++) list.add(first.getReference());
						an+=ploidy;
						gt= GenotypeBuilder.create(sampleName, list);
						}
					else
						{
						gt = GenotypeBuilder.createMissing(sampleName,this.ploidy);
						}	
					
					sample2genotypes.put(sampleName,gt);
					}
				if(with_an) vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY, an);
				if(with_ac) vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, 
						alleles.subList(1, alleles.size()).stream().mapToInt(A->(int)ac.count(A)).toArray()
						);

				if(with_af && an>0) {
					final double finalAn = an;
					vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, 
							alleles.subList(1, alleles.size()).stream().mapToDouble(A->ac.count(A)/finalAn).toArray()
							);

					}
				if(with_dp && dp>=0) {
					vcb.attribute(VCFConstants.DEPTH_KEY,dp);
					}
				vcb.genotypes(sample2genotypes.values());
				
				w.add(vcb.make());
				}
			eqiter.close();
			
			CloserUtil.close(w);w=null;
			array.cleanup();array=null;
			CloserUtil.close(iter);iter=null;
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(iter);
			if(array!=null) array.cleanup();
			}
		}
	

	public static void main(final String[] args)
		{
		new VCFMerge().instanceMainWithExit(args);
		}
	
	}
