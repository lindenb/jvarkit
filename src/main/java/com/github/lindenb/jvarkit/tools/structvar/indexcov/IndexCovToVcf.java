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
package com.github.lindenb.jvarkit.tools.structvar.indexcov;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.primitive.FloatArray;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Input

input is a tab-delimited file created by e.g: indexcov (https://github.com/brentp/goleft/tree/master/indexcov)

```
#chrom  start  end     SampleBB  SampleBC  SampleBD  SampleBE  SampleBF  SampleBG  SampleBH
chr1    23778  40778   1.59      1.31      1.67      1.61      1.83      1.52      1.48
chr1    29106  46106   1.9       1.54      1.72      1.97      1.88      1.53      1.95
chr1    84581  101581  0.764     0.841     1.2       1.16      1.18      1.13      1.23
chr1    15220  32220   0.355     0.704     1.09      0.784     0.81      1.37      0.954
chr1    58553  75553   0.353     0.436     0.912     0.836     1.16      1.09      0.611
chr1    19347  36347   0.381     0.411     0.811     0.795     1.16      1.22      0.495
chr1    81062  98062   1.09      0.972     1.35      1.22      1.66      1.76      1.1
chr1    17353  34353   1.06      1.06      1.23      1.26      1.44      1.43      1.03
chr1    48498  65498   1.08      0.996     1.28      1.44      1.52      1.57      1.05
```

output:

```
##fileformat=VCFv4.2
##FILTER=<ID=ALL_DEL,Description="number of samples >1 and all are deletions">
##FILTER=<ID=ALL_DUP,Description="number of samples >1 and all are duplication">
##FILTER=<ID=NO_SV,Description="There is no DUP or DEL in this variant">
##FORMAT=<ID=DEL,Number=1,Type=Integer,Description="set to 1 if relative number of copy <= 0.6">
##FORMAT=<ID=DUP,Number=1,Type=Integer,Description="set to 1 if relative number of copy >= 1.9">
##FORMAT=<ID=F,Number=1,Type=Float,Description="Relative number of copy: 0.5 deletion 1 normal 2.0 duplication">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=NDEL,Number=1,Type=Integer,Description="Number of samples being deleted">
##INFO=<ID=NDUP,Number=1,Type=Integer,Description="Number of samples being duplicated">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SampleBB	SampleBC	SampleBD	SampleBE	SampleBF	(...)
chr1	0	.	N	<DUP>	.	.	END=16384;NDEL=0;NDUP=8	GT:DUP:F	0:0:1.59	0:0:1.31	0:0:1.67	0:0:1.61	0:0:1.83 (...)
```


END_DOC
 */
@Program(
		name="indexcov2vcf",
		description="convert indexcov data to vcf",
		keywords={"cnv","duplication","deletion","sv","indexcov"},
		creationDate = "20200528",
		modificationDate="20400313",
		jvarkit_amalgamion = true,
		menu="CNV/SV"
		)
public class IndexCovToVcf extends Launcher {
	private static final Logger LOG = Logger.of(IndexCovToVcf.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-t","--treshold"},description=IndexCovUtils.TRESHOLD_OPT_DESC+". " +FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double treshold = IndexCovUtils.DEFAULT_TRESHOLD;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path refFile = null;
	@Parameter(names={"--limit"},description="limit to x lines for debugging",hidden = true)
	private long max_line_count=-1;
	@Parameter(names={"--no-merge"},description="disable adjacent block merging for the same variant (keep the original bed structure)")
	private boolean disable_block_merge=false;

	@ParametersDelegate
	private CasesControls caseControls=new CasesControls();
	@ParametersDelegate
	private WritingVariantsDelegate writingDelegate = new WritingVariantsDelegate();
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection= new WritingSortingCollection();
	
	
	
	private static class GenotypeV {
		int tid;
		int chromStart;
		int chromEnd;
		int sampleIdx;
		float normDepth;
		}
	
	
	
	private static class GenotypeVCodec extends AbstractDataCodec<GenotypeV> {
		@Override
		public AbstractDataCodec<GenotypeV> clone() {
			return new GenotypeVCodec();
		}
		@Override
		public GenotypeV decode(DataInputStream dis) throws IOException {
			GenotypeV g = new GenotypeV();
			try {
				g.tid = dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			g.chromStart = dis.readInt();
			g.chromEnd = dis.readInt();
			g.sampleIdx = dis.readInt();
			g.normDepth = dis.readFloat();
			return g;
		}
		@Override
		public void encode(DataOutputStream dos, GenotypeV g) throws IOException {
			dos.writeInt(g.tid);
			dos.writeInt(g.chromStart);
			dos.writeInt(g.chromEnd);
			dos.writeInt(g.sampleIdx);
			dos.writeFloat(g.normDepth);
		}
	}
	
	public IndexCovToVcf() {
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		final IndexCovUtils indexCovUtils = new IndexCovUtils(this.treshold);
		BufferedReader br = null;
		VariantContextWriter vcw  = null;
		TabixReader tabix = null;
		SortingCollection<GenotypeV> sorter1=null;
		SortingCollection<GenotypeV> sorter2=null;
		try {

			final SAMSequenceDictionary dictionary =  SequenceDictionaryUtils.extractRequired(this.refFile);
			
			final String input = oneFileOrNull(args);
			if(input!=null && input.endsWith(".bed.gz")) {
				Path tbi = Paths.get(input+FileExtensions.TABIX_INDEX);
				if(Files.exists(tbi)) {
					tabix = new TabixReader(input);
					}
				}
			if(tabix==null) {
				LOG.warn("Cannot instantiate a tabix reader= using stdin or non-indexed file");
				}
			
			br = super.openBufferedReader(input);
			String line = br.readLine();
			if(line==null) {		
				LOG.error( "Cannot read first line of input");
				return -1;
				}
			 final FileHeader tabixFileHeader = new FileHeader(line,S->Arrays.asList(CharSplitter.TAB.split(S)));
			 tabixFileHeader.
			 	assertColumn("#chrom",0).
			 	assertColumn("start",1).
			 	assertColumn("end",2)
			 	;
							 
			
			final Function<String,List<GenotypeV>> linetoGenotype= LINE-> {
				final List<GenotypeV> genotypes=new ArrayList<>(tabixFileHeader.size());
				if(StringUtil.isBlank(LINE)) return genotypes;
				final List<String> tokens =  tabixFileHeader.split(LINE);
				final SAMSequenceRecord ssr = dictionary.getSequence(tokens.get(0));
				if(ssr==null) {
					throw new JvarkitException.ContigNotFoundInDictionary(tokens.get(0),dictionary);
					}
				
				for(int i=3;i< tokens.size();i++) {
					final float normDepth = Float.parseFloat(tokens.get(i));
					
					 if(!IndexCovUtils.isValidNumber(normDepth)) {
							LOG.error("Bad fold "+normDepth+" for sample "+ tabixFileHeader.get(i)+" in "+LINE);
							continue;
						 	}
					
					
					
					final GenotypeV g = new GenotypeV();
					g.chromStart = Integer.parseInt(tokens.get(1));
					g.chromEnd = Integer.parseInt(tokens.get(2));
					g.normDepth = normDepth ;
					g.sampleIdx = i;
					
					
					if(g.chromEnd>ssr.getSequenceLength()) {
						LOG.warn("WARNING sequence length in "+LINE+" is greater than in dictionary ");
						}
					
					g.tid = ssr.getSequenceIndex();
					genotypes.add(g);
					}
				return genotypes;
				};
			
			final CloseableIterator<List<GenotypeV>> iterator_by_loc;
			
			if(!this.disable_block_merge) {
				long nLine=0L;
				
				final Comparator<GenotypeV> sortByPosition = (A,B)->{
					int i= Integer.compare(A.tid, B.tid);
					if(i!=0) return i;
					i= Integer.compare(A.chromStart, B.chromStart);
					if(i!=0) return i;
					return Integer.compare(A.chromEnd, B.chromEnd);
					};
					
				final Comparator<GenotypeV> sortBySamplePos = (A,B)->{
					int i= Integer.compare(A.sampleIdx, B.sampleIdx);
					if(i!=0) return i;
					return sortByPosition.compare(A,B);
					};	
				
				sorter1 =  SortingCollection.newInstance(
						GenotypeV.class,
						new GenotypeVCodec(),
						sortBySamplePos,
						writingSortingCollection.getMaxRecordsInRam(),
						writingSortingCollection.getTmpPaths()
						);
				sorter1.setDestructiveIteration(true);

				
				
				while((line=br.readLine())!=null) {
					++nLine;
					if(this.max_line_count>=0 && nLine>this.max_line_count) break;
					if(StringUtil.isBlank(line)) continue;
					for(GenotypeV g: linetoGenotype.apply(line)) {
						if(!indexCovUtils.isVariant(g.normDepth)) continue;
						sorter1.add(g);
						}
					}
				sorter1.doneAdding();
				
				
				sorter2 =  SortingCollection.newInstance(
						GenotypeV.class,
						new GenotypeVCodec(),
						sortByPosition,
						writingSortingCollection.getMaxRecordsInRam(),
						writingSortingCollection.getTmpPaths()
						);
				sorter2.setDestructiveIteration(true);
				
				try(CloseableIterator<GenotypeV> iter0=sorter1.iterator()) {
					try(final EqualRangeIterator<GenotypeV> eq = new EqualRangeIterator<>(iter0,(A,B)->{
						int i= Integer.compare(A.sampleIdx, B.sampleIdx);
						if(i!=0) return i;
						return  Integer.compare(A.tid, B.tid);
						})) {
							while(eq.hasNext()) {
								final List<GenotypeV> array = new ArrayList<>(eq.next());
								while(!array.isEmpty()) {
									final GenotypeV first = array.remove(0);
									int i=0;
									while(!array.isEmpty() && !disable_block_merge) {
										final GenotypeV o = array.get(i);
										if(o.tid!=first.tid) throw new IllegalStateException();
										if(o.sampleIdx!=first.sampleIdx) throw new IllegalStateException();
										if(o.chromStart<first.chromStart) throw new IllegalStateException();
										if(first.chromEnd!=o.chromStart) break;
										array.remove(0);
										first.chromEnd = o.chromEnd;
										}
									sorter2.add(first);
									}
							}
						}
					}
				
				sorter2.doneAdding();
				sorter1.cleanup();
				sorter1=null;
				
				final CloseableIterator<GenotypeV> iter0=sorter2.iterator();
				iterator_by_loc = new EqualRangeIterator<>(iter0,sortByPosition);
				}
			else
				{
				final BufferedReader fbr=br;
				iterator_by_loc = new AbstractCloseableIterator<List<GenotypeV>>() {
					int nLine=0;
					@Override
					protected List<GenotypeV> advance() {
						try {
							final String line = fbr.readLine();
							if(line==null) return null;
							if(max_line_count>=0 && nLine>max_line_count) return null;
							nLine++;
							return linetoGenotype.apply(line);
							}
						catch(IOException err) {
							throw new RuntimeIOException(err);
							}
						}
					public void close() {
						
						}
					};
				}
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true,VCFConstants.GENOTYPE_KEY,VCFConstants.GENOTYPE_QUALITY_KEY);
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.END_KEY);
			
			/** raw value in indexcov */
			final VCFFormatHeaderLine foldHeader = new VCFFormatHeaderLine("F", 1, VCFHeaderLineType.Float,"Relative number of copy: 0.5 deletion 1 normal 2.0 duplication");
			metaData.add(foldHeader);

			final VCFFilterHeaderLine filterAllDel = new VCFFilterHeaderLine("ALL_DEL", "number of samples greater than 1 and all are deletions");
			metaData.add(filterAllDel);
			final VCFFilterHeaderLine filterAllDup = new VCFFilterHeaderLine("ALL_DUP", "number of samples  greater than  1 and all are duplication");
			metaData.add(filterAllDup);
			final VCFFilterHeaderLine filterNoSV= new VCFFilterHeaderLine("NO_SV", "There is no DUP or DEL in this variant");
			metaData.add(filterNoSV);
			final VCFFilterHeaderLine filterHomDel = new VCFFilterHeaderLine("HOM_DEL", "There is one Homozygous deletion.");
			metaData.add(filterHomDel);
			final VCFFilterHeaderLine filterHomDup = new VCFFilterHeaderLine("HOM_DUP", "There is one Homozygous duplication.");
			metaData.add(filterHomDup);

			final VCFInfoHeaderLine infoSVLEN = new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer, "SV length");
			metaData.add(infoSVLEN);
			final VCFInfoHeaderLine infoNumDup = new VCFInfoHeaderLine("NDUP", 1, VCFHeaderLineType.Integer,"Number of samples being duplicated");
			metaData.add(infoNumDup);
			final VCFInfoHeaderLine infoNumDel = new VCFInfoHeaderLine("NDEL", 1, VCFHeaderLineType.Integer,"Number of samples being deleted");
			metaData.add(infoNumDel);
			final VCFInfoHeaderLine infoSingleton = new VCFInfoHeaderLine("SINGLETON", 1, VCFHeaderLineType.Flag,"Singleton candidate");
			metaData.add(infoSingleton);
			final VCFInfoHeaderLine infoAllAffected = new VCFInfoHeaderLine("ALL_CASES", 1, VCFHeaderLineType.Flag,"All cases are affected");
			metaData.add(infoAllAffected);
			final VCFInfoHeaderLine infoFisher = new VCFInfoHeaderLine("FISHER", 1, VCFHeaderLineType.Float,"Fisher test case/control");
			metaData.add(infoFisher);
			metaData.add(new VCFInfoHeaderLine("MEDIAN_ALL", 1, VCFHeaderLineType.Float, "Median normalized depth"));
			metaData.add(new VCFInfoHeaderLine("MEAN_ALL", 1, VCFHeaderLineType.Float,"Mean normalized depth"));
			metaData.add(new VCFInfoHeaderLine("STDDEV_ALL", 1, VCFHeaderLineType.Float,"Stddev normalized depth"));
			metaData.add(new VCFInfoHeaderLine("MEDIAN_VAR", 1, VCFHeaderLineType.Float, "Median normalized depth for Variants"));
			metaData.add(new VCFInfoHeaderLine("MEAN_VAR", 1, VCFHeaderLineType.Float,"Mean normalized depth for Variants"));
			metaData.add(new VCFInfoHeaderLine("STDDEV_VAR", 1, VCFHeaderLineType.Float,"Stddev normalized depth for Variants"));
			metaData.add(new VCFInfoHeaderLine("MEDIAN_REF", 1, VCFHeaderLineType.Float, "Median normalized depth for REF"));
			metaData.add(new VCFInfoHeaderLine("MEAN_REF", 1, VCFHeaderLineType.Float,"Mean normalized depth for REF"));
			metaData.add(new VCFInfoHeaderLine("STDDEV_REF", 1, VCFHeaderLineType.Float,"Stddev normalized depth for REF"));
			metaData.add(new VCFInfoHeaderLine("NBLOCKS", 1, VCFHeaderLineType.Integer,"number of indexcov BLOCKS (size="));

			
			final List<String> samples =tabixFileHeader.subList(3,tabixFileHeader.size());
			caseControls.load().retain(samples);
			
			
			final VCFHeader vcfHeader = new VCFHeader(metaData, samples);
	  		JVarkitVersion.getInstance().addMetaData(this, vcfHeader);

			
			vcfHeader.setSequenceDictionary(dictionary);
			
			
			vcw = this.writingDelegate.dictionary(dictionary).open(outputFile);
			vcw.writeHeader(vcfHeader);
			
			//final List<Allele> NO_CALL_NO_CALL = Arrays.asList(Allele.NO_CALL,Allele.NO_CALL);
			final Allele DUP_ALLELE =Allele.create("<DUP>",false);
			final Allele DEL_ALLELE =Allele.create("<DEL>",false);
			final Allele REF_ALLELE =Allele.create("N",true);

			while(iterator_by_loc.hasNext()) {
				final List<GenotypeV> array = new ArrayList<>(iterator_by_loc.next());
				final GenotypeV first = array.get(0);
				
				
				final Set<Allele> alleles =  new HashSet<>();
				alleles.add(REF_ALLELE);
				final String contig= dictionary.getSequence(first.tid).getSequenceName();
				final VariantContextBuilder vcb = new VariantContextBuilder();
				vcb.chr(contig);
				vcb.start(first.chromStart+1);
				vcb.stop( first.chromEnd);
				vcb.attribute(VCFConstants.END_KEY,  first.chromEnd);
				vcb.attribute(infoSVLEN.getID(), first.chromEnd-first.chromStart);
				vcb.attribute("NBLOCKS", (first.chromEnd-first.chromStart)/IndexCovUtils.BAI_BLOCK_SIZE);
				
				
				final Map<String,Float> sample2fold = new HashMap<>(samples.size());
				for(GenotypeV g:array) {
					final String sampleName = tabixFileHeader.get(g.sampleIdx);
					sample2fold.put(sampleName, g.normDepth);
					}

				if(tabix!=null) {
					final AutoMap<String, Sum,Sum> sample2sum = AutoMap.make(()->new Sum());
					final FloatArray floatArrayAll =new FloatArray();
					final FloatArray floatArrayWithVariant =new FloatArray();
					final FloatArray floatArrayWithNoVariant =new FloatArray();
					TabixReader.Iterator tbxr =tabix.query(contig,first.chromStart+1, first.chromEnd);
					for(;;) {
						line = tbxr.next();
						if(line==null) break;
						final List<String> tokens = tabixFileHeader.split(line);
						for(int i=3;i<tokens.size();i++) {
							final float normDepth = Float.parseFloat(tokens.get(i));
							if(!IndexCovUtils.isValidNumber(normDepth)) continue;
							floatArrayAll.add(normDepth);
							final String sn =  tabixFileHeader.get(i);
							if(!sample2fold.containsKey(sn)) {
								sample2sum.insert(sn).increment(normDepth);
								}
							
							if(indexCovUtils.isVariant(normDepth)) {
								floatArrayWithVariant.add(normDepth);
								}
							else {
								floatArrayWithNoVariant.add(normDepth);
								}
							}
						}
					double[] values;
					if(!floatArrayAll.isEmpty()) {
						values= floatArrayAll.stream().mapToDouble(X->X).toArray();
						vcb.attribute("MEDIAN_ALL",new Median().evaluate(values));
						vcb.attribute("MEAN_ALL",new Mean().evaluate(values));
						vcb.attribute("STDDEV_ALL",new StandardDeviation().evaluate(values));
						}
					if(!floatArrayWithVariant.isEmpty()) {
						values= floatArrayWithVariant.stream().mapToDouble(X->X).toArray();
						vcb.attribute("MEDIAN_VAR",new Median().evaluate(values));
						vcb.attribute("MEAN_VAR",new Mean().evaluate(values));
						vcb.attribute("STDDEV_VAR",new StandardDeviation().evaluate(values));
						}
					if(!floatArrayWithNoVariant.isEmpty()) {
						values= floatArrayWithNoVariant.stream().mapToDouble(X->X).toArray();
						vcb.attribute("MEDIAN_REF",new Median().evaluate(values));
						vcb.attribute("MEAN_REF",new Mean().evaluate(values));
						vcb.attribute("STDDEV_REF",new StandardDeviation().evaluate(values));
						}
					for(final String sn:sample2sum.keySet()) {
						final Sum sum = sample2sum.get(sn);
						sample2fold.put(sn,(float)(sum.getResult()/sum.getN()));
						}
					}
				
				
				int count_dup=0;
				int count_del=0;
				
				
				
				final List<Genotype> genotypes = new ArrayList<>(samples.size());
				
				int count_sv_cases = 0;
				int count_sv_controls = 0;
				int count_ref_cases = 0;
				int count_ref_controls = 0;

				for(final String sampleName:sample2fold.keySet())
					{
					final float normDepth = sample2fold.get(sampleName);
					final IndexCovUtils.SvType type = indexCovUtils.getType(normDepth);
					
					final GenotypeBuilder gb;

					switch(type) {
						case REF: {
							gb = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,REF_ALLELE));
							break;
							}
						case HET_DEL:{
							gb = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,DEL_ALLELE));						
							alleles.add(DEL_ALLELE);
							count_del++;
							break;
							}
						case HOM_DEL: {
							gb = new GenotypeBuilder(sampleName,Arrays.asList(DEL_ALLELE,DEL_ALLELE));						
							alleles.add(DEL_ALLELE);
							count_del++;
							vcb.filter(filterHomDel.getID());
							break;
							}
						case HET_DUP: {
							gb = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,DUP_ALLELE));
							alleles.add(DUP_ALLELE);
							count_dup++;
							break;
							}
						case HOM_DUP: {
							gb = new GenotypeBuilder(sampleName,Arrays.asList(DUP_ALLELE,DUP_ALLELE));
							alleles.add(DUP_ALLELE);
							vcb.filter(filterHomDup.getID());
							count_dup++;
							break;
							}
						default:
							{
							gb = new GenotypeBuilder(sampleName,Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
							break;
							}
						}
					gb.attribute(foldHeader.getID(),normDepth);
					gb.GQ(type.getGenotypeQuality(normDepth));

					if(type.isVariant()) {
						if(caseControls.isCase(sampleName))
							{
							count_sv_cases++;
							}
						else if(caseControls.isControl((sampleName)))
							{
							count_sv_controls++;
							}
						}
					else
						{
						if(caseControls.isCase(sampleName))
							{
							count_ref_cases++;
							}
						else if(caseControls.isControl((sampleName)))
							{
							count_ref_controls++;
							}
						}
										
					genotypes.add(gb.make());
					}
				
				
				vcb.alleles(alleles);
				
				if(!caseControls.isEmpty()) {
					if(		count_sv_cases==1 && 
							count_sv_controls==0
						) {
						vcb.attribute(infoSingleton.getID(),Boolean.TRUE);
					} else if( count_sv_cases>0 && 
							count_sv_cases == caseControls.getCases().size() &&
							count_sv_controls==0
						) {
						vcb.attribute(infoAllAffected.getID(),Boolean.TRUE);
						}

					vcb.attribute(infoFisher.getID(),FisherExactTest.compute(
							count_sv_cases, count_sv_controls,
							count_ref_cases, count_ref_controls
							).getAsDouble());
					}
					
				
				vcb.genotypes(genotypes);
				
				if(count_dup == samples.size() && samples.size()!=1) {
					vcb.filter(filterAllDup.getID());
				}
				if(count_del == samples.size() && samples.size()!=1) {
					vcb.filter(filterAllDel.getID());
				}
				
				
				
				vcb.attribute(infoNumDel.getID(), count_del);
				vcb.attribute(infoNumDup.getID(), count_dup);
				
				
				if(alleles.size()==1) continue;
				
				vcw.add(vcb.make());
				}
			iterator_by_loc.close();
			vcw.close();
			vcw=null;
			br.close();
			br=null;
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(br!=null) try {br.close();} catch(Throwable err) {}
			if(vcw!=null) try {vcw.close();} catch(Throwable err) {}
			if(tabix!=null) try {tabix.close();} catch(Throwable err) {}
			if(sorter1!=null) try {sorter1.cleanup();} catch(Throwable err) {}
			if(sorter2!=null) try {sorter2.cleanup();} catch(Throwable err) {}
			}
		}

		public static void main(final String[] args) {
			new IndexCovToVcf().instanceMainWithExit(args);
			}

}
