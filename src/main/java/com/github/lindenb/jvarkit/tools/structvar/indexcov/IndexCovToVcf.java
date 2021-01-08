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
package com.github.lindenb.jvarkit.tools.structvar.indexcov;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.IndexCovUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
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

## history

  * 20191112 : add pedigree

END_DOC
 */
@Program(
		name="indexcov2vcf",
		description="convert indexcov data to vcf",
		keywords={"cnv","duplication","deletion","sv","indexcov"},
		modificationDate="20200227"
		)
public class IndexCovToVcf extends Launcher {
	private static final Logger LOG = Logger.build(IndexCovToVcf.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-t","--treshold"},description=IndexCovUtils.TRESHOLD_OPT_DESC+". " +FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double treshold = IndexCovUtils.DEFAULT_TRESHOLD;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path refFile = null;
	@Parameter(names={"-p","--pedigree"},description="Optional Pedigree. "+PedigreeParser.OPT_DESC)
	private Path pedFile = null;

	@ParametersDelegate
	private WritingVariantsDelegate writingDelegate = new WritingVariantsDelegate();
	
	
	public IndexCovToVcf() {
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		final IndexCovUtils indexCovUtils = new IndexCovUtils(this.treshold);
		
		final CharSplitter tab = CharSplitter.TAB;
		BufferedReader r = null;
		VariantContextWriter vcw  = null;
		try {

			
			final SAMSequenceDictionary dict;
			if(this.refFile==null) {
				dict = null;
			} else
			{
				dict= SequenceDictionaryUtils.extractRequired(this.refFile);
			}
			
			
			
			r = super.openBufferedReader(oneFileOrNull(args));
			String line = r.readLine();
			if(line==null) {		
				LOG.error( "Cannot read first line of input");
				return -1;
				}
			String tokens[] = tab.split(line);
			if(tokens.length<4 ||
				!tokens[0].equals("#chrom") ||
				!tokens[1].equals("start") ||
				!tokens[2].equals("end")) {
				LOG.error( "bad first line "+line );
				return -1;
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

			
			final VCFInfoHeaderLine infoNumDup = new VCFInfoHeaderLine("NDUP", 1, VCFHeaderLineType.Integer,"Number of samples being duplicated");
			metaData.add(infoNumDup);
			final VCFInfoHeaderLine infoNumDel = new VCFInfoHeaderLine("NDEL", 1, VCFHeaderLineType.Integer,"Number of samples being deleted");
			metaData.add(infoNumDel);
			final VCFInfoHeaderLine infoSingleton = new VCFInfoHeaderLine("SINGLETON", 1, VCFHeaderLineType.Flag,"Singleton candidate");
			metaData.add(infoSingleton);
			final VCFInfoHeaderLine infoAllAffected = new VCFInfoHeaderLine("ALL_CASES", 1, VCFHeaderLineType.Flag,"All cases are affected");
			metaData.add(infoAllAffected);
			
			
			final List<String> samples = Arrays.asList(tokens). subList(3,tokens.length);
			final Pedigree pedigree;
			
			final int count_cases_in_pedigree;
			if(this.pedFile==null) {
				pedigree = PedigreeParser.empty();
				count_cases_in_pedigree = 0;
			} else
			{
				pedigree = new PedigreeParser().parse(this.pedFile);
				final Set<String> set = new HashSet<>(samples);
				count_cases_in_pedigree = (int)pedigree.getAffectedSamples().stream().filter(S->set.contains(S.getId())).count();
			}

			
			final VCFHeader vcfHeader = new VCFHeader(metaData, samples);
	  		JVarkitVersion.getInstance().addMetaData(this, vcfHeader);

			
			if(dict!=null) {
				vcfHeader.setSequenceDictionary(dict);
			}
			
			vcw = this.writingDelegate.dictionary(dict).open(outputFile);
			vcw.writeHeader(vcfHeader);
			
			//final List<Allele> NO_CALL_NO_CALL = Arrays.asList(Allele.NO_CALL,Allele.NO_CALL);
			final Allele DUP_ALLELE =Allele.create("<DUP>",false);
			final Allele DEL_ALLELE =Allele.create("<DEL>",false);
			final Allele REF_ALLELE =Allele.create("N",true);

			while((line=r.readLine())!=null) {
				if(StringUtil.isBlank(line)) continue;
				tokens =  tab.split(line);
				if(tokens.length!=3+samples.size()) {
					r.close();
					vcw.close();
					throw new JvarkitException.TokenErrors("expected "+(samples.size()+3)+ "columns.", tokens);
					}
				
				final Set<Allele> alleles =  new HashSet<>();
				alleles.add(REF_ALLELE);
				
				final VariantContextBuilder vcb = new VariantContextBuilder();
				vcb.chr(tokens[0]);
				vcb.start(Integer.parseInt(tokens[1]));
				final int chromEnd = Integer.parseInt(tokens[2]);
				vcb.stop(chromEnd);
				vcb.attribute(VCFConstants.END_KEY, chromEnd);
				
				if(dict!=null) {
					final SAMSequenceRecord ssr = dict.getSequence(tokens[0]);
					if(ssr==null) {
						LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(tokens[0],dict));
						return -1;
					}
					if(chromEnd>ssr.getSequenceLength()) {
						LOG.warn("WARNING sequence length in "+line+" is greater than in dictionary ");
					}
				}
				
				int count_dup=0;
				int count_del=0;
				final Map<String,Float> sample2fold = new HashMap<>(samples.size());
				for(int i=3;i<tokens.length;i++) {
					final String sampleName = samples.get(i-3);
					final float f = Float.parseFloat(tokens[i]);
					 if(f<0 || Float.isNaN(f) ||! Float.isFinite(f)) {
						 LOG.error("Bad fold "+f+" for sample "+sampleName+" in "+line);
					 	}
					sample2fold.put(sampleName, f);
					}
				
				
				final List<Genotype> genotypes = new ArrayList<>(samples.size());
				
				int count_sv_cases = 0;
				int count_sv_controls = 0;
				int count_ref_cases = 0;
				int count_ref_controls = 0;

				boolean got_sv = false;
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
					if(type.isVariant()) got_sv=true;
					gb.attribute(foldHeader.getID(),normDepth);
					gb.GQ(type.getGenotypeQuality(normDepth));

					final Sample sn = pedigree.getSampleById(sampleName);
					if(sn!=null) {
						if(type.isVariant()) {
							if(sn.isAffected())
								{
								count_sv_cases++;
								}
							else if(sn.isUnaffected())
								{
								count_sv_controls++;
								}
							}
						else
							{
							if(sn.isAffected())
								{
								count_ref_cases++;
								}
							else if(sn.isUnaffected())
								{
								count_ref_controls++;
								}
							}
						}
					
					
					genotypes.add(gb.make());
					}
				
				
				vcb.alleles(alleles);
				
				if(!pedigree.isEmpty() &&
						count_sv_cases==1 && 
						count_ref_cases>0 &&
						count_sv_controls==0 &&
						count_ref_controls>0
					) {
					vcb.attribute(infoSingleton.getID(),Boolean.TRUE);
				} else if(!pedigree.isEmpty() &&
						count_sv_cases>0 && 
						count_sv_cases == count_cases_in_pedigree &&
						count_ref_cases==0 &&
						count_sv_controls==0 &&
						count_ref_controls>0
					) {
					vcb.attribute(infoAllAffected.getID(),Boolean.TRUE);
				}
				
					
				
				vcb.genotypes(genotypes);
				
				if(count_dup == samples.size() && samples.size()!=1) {
					vcb.filter(filterAllDup.getID());
				}
				if(count_del == samples.size() && samples.size()!=1) {
					vcb.filter(filterAllDel.getID());
				}
				
				if(!got_sv) {
					vcb.filter(filterNoSV.getID());
					}
				
				vcb.attribute(infoNumDel.getID(), count_del);
				vcb.attribute(infoNumDup.getID(), count_dup);
				
				vcw.add(vcb.make());
				}
			vcw.close();
			vcw=null;
			r.close();
			r=null;
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(vcw);
			}
		}

		public static void main(final String[] args) {
			new IndexCovToVcf().instanceMainWithExit(args);
			}

}
