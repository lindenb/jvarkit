/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.misc.VcfCadd;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.


### Output


#### INFO column

 *  BurdenF1Fisher : Fisher test

#### FILTER column

 *  BurdenF1Fisher :Fisher test doesn't meet  user's requirements

### see also

 *  VcfBurdenFilter3


END_DOC
*/


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC


### Output

#### INFO column


 *  BurdenF1Fisher : Fisher test

#### FILTER column

 *  BurdenF1Fisher :Fisher test doesn't meet  user's requirements

### see also

 *  VcfBurdenFilter3

### History 

  * [20180910] add NFE https://www.youtube.com/watch?v=Fi5dLGAH8R0
  * [20180831] add CADD values from VcfCadd

END_DOC
*/

@Program(name="vcfburdenrscriptv",
	description="Fisher Case / Controls per Variant (Vertical)",
	keywords={"vcf","burden","fisher"}
	)
public class VcfBurdenRscriptV
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurdenRscriptV.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-if","--ignorefilter"},description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
	private boolean acceptFiltered = false;

	@Parameter(names={"-f","--function"},description="User defined R function to be called after each VCF")
	private String userDefinedFunName = "";

	@Parameter(names={"-t","--title"},description="Try to find ##(TITLE)=abcdefghijk in the VCF header and use it as the name of the VCF chunk")
	private String titleHeaderStr = "";

	@Parameter(names={"-minusnineiszero","--minusnineiszero"},description="No Call is '0' (default is -9)")
	private boolean nocalliszero = false;
	
	@Parameter(names={"--cadd","-cadd"},description="[20180831] Include CADD data, if available (INFO/"+VcfCadd.DEFAULT_CADD_FLAG_PHRED +" INFO/"+VcfCadd.DEFAULT_CADD_FLAG_SCORE+")")
	private boolean include_vcf_cadd = false;
	@Parameter(names={"-cpm","--cadd-phred-missing"},description="[20180831] value for CADD / phred missing data")
	private String cadd_phred_missing_value = "NA";
	@Parameter(names={"-csm","--cadd-score-missing"},description="[20180921] value for CADD / score missing data")
	private String cadd_score_missing_value = "NA";
	@Parameter(names={"--nfe"},description="[20180910] INCLUDE gnomad genome NFE AC")
	private boolean include_gnomad_genome_nfe_ac  = false;
	@Parameter(names={"-p","--pedigree"},description="[20180831] pedigree file (or I will try to extract the pedigree from the vcf header. " +Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	
	
	
	public VcfBurdenRscriptV()
		{
		}
	
	private static class Variant {
		String contig;
		int start;
		int end;
		Allele ref;
		Allele alt;
		Double maf=null;
		Double cadd_phred = null;
		Double cadd_score = null;
		Integer gnomad_genome_nfe_ac = null;
	}
	
	@Override
	public int doWork(final List<String> args) {
		/* 
		 mail matile April 11: requires calling R for SKAT: needs
		 matrix of samples and genotypes */
		PrintWriter pw = null;
		VCFIterator in = null;
		LineIterator lr = null;
		long vcf_id = System.currentTimeMillis();
		try {
			String inputName = oneFileOrNull(args);
			lr =(
				inputName==null?
				IOUtils.openStreamForLineIterator(stdin()):
				IOUtils.openURIForLineIterator(inputName)
				);
			
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			pw.println("# This is the header of the generated R script. ");
			pw.println("# generated by "+ getProgramName());
			pw.println("# version "+getGitHash());
			

			
			if(this.userDefinedFunName!=null && !this.userDefinedFunName.trim().isEmpty()) {
				pw.println("# user defined function "+userDefinedFunName+" should have been defined BEFORE this header. ");
				pw.println("#  something like 'cat user_function.R this_file.R |  R --no-save > result.txt'");
			}
			
			
			
			if(!lr.hasNext()) {
				LOG.warn("No input found. Couldn't read any VCF header file");
			}
			
			while(lr.hasNext()) {
			vcf_id++;
			in = VCFUtils.createVCFIteratorFromLineIterator(lr,true);
			
			final VCFHeader header=in.getHeader();
			final Set<Pedigree.Person> samples;
			
			if(this.pedigreeFile==null) {
				samples = new TreeSet<>( super.getCasesControlsInPedigree(header));
				}
			else
				{
				final Pedigree ped = Pedigree.newParser().parse(this.pedigreeFile);
				samples = new TreeSet<>(new Pedigree.CaseControlExtractor().extract(header, ped));
				}
			
			final List<Variant> variants = new ArrayList<>();

		

			
			boolean first=true;
			pw.println("# BEGIN  VCF ##########################################");
			pw.println("# Title ");
			final VCFHeaderLine vcfTitle= 
					(
					StringUtil.isBlank(this.titleHeaderStr) ?
					null :
					header.getOtherHeaderLine(this.titleHeaderStr.trim())
					);
			if(vcfTitle==null) {
				LOG.warn("No title was found");
				pw.println("# [WARNING] No title was found. ");
				
				pw.println("vcf.title <- \"vcf"+String.format("%04d", vcf_id)+"\"");
			} else {
				pw.println("vcf.title <- \""+vcfTitle.getValue()+"\"");
			}
			
			
			first=true;
			pw.println("# samples ( 0: unaffected 1:affected)");
			
			pw.print("population <- data.frame(family=c(");			
			first=true;for(final Pedigree.Person person : samples) { if(!first) pw.print(","); pw.print("\""+person.getFamily().getId()+"\"");first=false;}
			pw.print("),name=c(");
			first=true;for(final Pedigree.Person person : samples) { if(!first) pw.print(","); pw.print("\""+person.getId()+"\"");first=false;}
			pw.print("),status=c(");
			first=true;for(final Pedigree.Person person : samples) { if(!first) pw.print(","); pw.print(person.isUnaffected()?0:1);first=false;}
			pw.println("))");
						
			
			first=true;
			pw.println();
			pw.println("# genotypes as a list. Should be a multiple of length(samples).");
			pw.println("# 0 is homref (0/0), 1 is het (0/1), 2 is homvar (1/1)");
			pw.println("# if the variant contains another ALT allele: (0/2) and (2/2) are considered 0 (homref)");
			pw.print("genotypes <- c(");
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(in.hasNext()) {
				final VariantContext ctx = progess.watch(in.next());
				if(ctx.isFiltered() && !this.acceptFiltered) continue;
				final int n_alts = ctx.getAlternateAlleles().size();
				if( n_alts == 0) {
					LOG.warn("ignoring variant without ALT allele.");
					continue;
				}
				if( n_alts > 1) {
					LOG.warn("variant with more than one ALT. Using getAltAlleleWithHighestAlleleCount.");
				}
				final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
				final MafCalculator mafCalculator= new MafCalculator(observed_alt,ctx.getContig());
				for (final Pedigree.Person person : samples) {
						final Genotype genotype = ctx.getGenotype(person.getId());
						if(genotype==null) {
							pw.close();pw=null;
							in.close();
							throw new IllegalStateException("Cannot get genotype for "+person.getId());
						}
						
						mafCalculator.add(genotype,  person.isMale());
						
						
						
						if(!first) pw.print(",");
						if (genotype.isHomRef()) {
							pw.print('0');
						} else if (genotype.isHomVar() && genotype.getAlleles().contains(observed_alt)) {
							pw.print('2');
						} else if (genotype.isHet() && 
								genotype.getAlleles().contains(observed_alt) &&
								genotype.getAlleles().contains(ctx.getReference())) {
							pw.print('1');

						}
						/* we treat 0/2 has hom-ref */
						else if (genotype.isHet() && 
								!genotype.getAlleles().contains(observed_alt) &&
								genotype.getAlleles().contains(ctx.getReference())) {
							LOG.warn("Treating "+genotype+" as hom-ref (0) alt="+observed_alt);
							pw.print('0');
							
						} 
						/* we treat 2/2 has hom-ref */
						else if (genotype.isHomVar() && !genotype.getAlleles().contains(observed_alt)) {
							LOG.warn("Treating "+genotype+" as hom-ref (0) alt="+observed_alt);
							pw.print('0');
						}
						else {
							pw.print(this.nocalliszero?"0":"-9");
						}
						first=false;
					}
				
				final Variant variant = new Variant();
				variant.contig = ctx.getContig();
				variant.start = ctx.getStart();
				variant.end = ctx.getEnd();
				variant.ref = ctx.getReference();
				variant.alt = observed_alt;
				
				/** 20181110, include NFE AC */
				if(this.include_gnomad_genome_nfe_ac)
					{
					final Function<Object,Integer> object2int = (O)->{
						if(O==null) return null;
						final String str= String.valueOf(O);
						if(str.isEmpty() || str.equals(".")) return null;
						try {
							return new Integer(str);
							}
						catch (final NumberFormatException e) {
							return null;
							}
						};
					final int allele_index = ctx.getAlleleIndex(observed_alt);
					final List<Object> nfeL = ctx.getAttributeAsList("gnomad_genome_AC_NFE");
					if(!nfeL.isEmpty() && allele_index> 0 /* yes */ && allele_index< nfeL.size()) {
						variant.gnomad_genome_nfe_ac  = object2int.apply(nfeL.get(allele_index-1));/* -1 because allele_index is relative to REF */
						}
					}
				
				/** 2018-08-31 : Matilde wants scores CADD in output */
				if(this.include_vcf_cadd) {
					final Function<Object,Double> object2double = (O)->{
						if(O==null) return null;
						final String str= String.valueOf(O);
						if(str.isEmpty() || str.equals(".")) return null;
						try {
							return Double.valueOf(str);
							}
						catch (final NumberFormatException e) {
							LOG.warn("bad cadd number "+O);
							return null;
							}
						};
					final int allele_index = ctx.getAlleleIndex(observed_alt);
					List<Object> caddL = ctx.getAttributeAsList(VcfCadd.DEFAULT_CADD_FLAG_PHRED);
					if(!caddL.isEmpty() && allele_index> 0 /* yes */ && (allele_index-1)< caddL.size()) {
						variant.cadd_phred  = object2double.apply(caddL.get(allele_index-1));/* -1 because allele_index is relative to REF */
						}
					
					caddL = ctx.getAttributeAsList(VcfCadd.DEFAULT_CADD_FLAG_SCORE);
					if(!caddL.isEmpty() && allele_index> 0 /* yes */ && (allele_index-1)< caddL.size()) {
						variant.cadd_score  = object2double.apply(caddL.get(allele_index-1));/* -1 because allele_index is relative to REF */
						}
					}
				
				if(!mafCalculator.isEmpty())
					{
					variant.maf = mafCalculator.getMaf();
					}
				else
					{
					variant.maf = null;
					}
				variants.add(variant);
				}// end reading vcf
			progess.finish();
			in.close();
			
			pw.println(")");
			first = true;
			
			pw.print("# variants. CONTIG/START/END/REF/ALT/MAF");
			if(this.include_vcf_cadd)
				{
				pw.print("/CADD_SCORE/CADD_PHRED");
				}
			pw.println();
			pw.print("variants <- data.frame(chrom=c(");
			first=true; for(final Variant v: variants) {if(!first) pw.print(",");pw.print("\""+v.contig+"\"");first=false;}
			pw.print("),chromStart=c(");
			first=true; for(final Variant v: variants) {if(!first) pw.print(",");pw.print(v.start);first=false;}
			pw.print("),chromEnd=c(");
			first=true; for(final Variant v: variants) {if(!first) pw.print(",");pw.print(v.end);first=false;}
			pw.print("),refAllele=c(");
			first=true; for(final Variant v: variants) {if(!first) pw.print(",");pw.print("\""+v.ref.getDisplayString()+"\"");first=false;}
			pw.print("),altAllele=c(");
			first=true; for(final Variant v: variants) {if(!first) pw.print(",");pw.print("\""+v.alt.getDisplayString()+"\"");first=false;}
			pw.print("),maf=c(");
			first=true; for(final Variant v: variants) {if(!first) pw.print(",");pw.print(v.maf==null?"NA":String.valueOf(v.maf));first=false;}
			pw.print(")");
			
			if(this.include_gnomad_genome_nfe_ac)
				{
				pw.print(",gnomad_nfe_ac=c(");
				pw.print(variants.stream().
						map(V->V.gnomad_genome_nfe_ac).
						map(D->D==null?"0":String.valueOf(D)).
						collect(Collectors.joining(",")
						));
				pw.print(")");				
				}
			
			if(this.include_vcf_cadd)
				{
				pw.print(",cadd_phred=c(");
				pw.print(variants.stream().
						map(V->V.cadd_phred).
						map(D->D==null?cadd_phred_missing_value:String.valueOf(D)).
						collect(Collectors.joining(",")
						));
				pw.print("),cadd_score=c(");
				pw.print(variants.stream().
						map(V->V.cadd_score).
						map(D->D==null?cadd_score_missing_value:String.valueOf(D)).
						collect(Collectors.joining(",")
						));
				pw.print(")");
				}
			
			pw.println(")");
			
			if(!variants.isEmpty())
				{
				pw.println("# assert sizes");
				pw.println("stopifnot( length(genotypes) %% NROW(population) == 0 )");
				pw.println("stopifnot(NROW(variants) * NROW(population) == length(genotypes) )");
				
				
				if(this.userDefinedFunName==null || this.userDefinedFunName.trim().isEmpty()) {
					pw.println("## WARNING not user-defined R function was defined");
					}
				else
					{
					pw.println("# consumme data with user-defined R function ");
					pw.println(this.userDefinedFunName+"()");
					}
				}
			else
				{
				LOG.warn("No Variant found");
				}
			pw.println("# END VCF ##########################################");
			}
			
			
			
			pw.flush();
			if(pw.checkError()) {
				LOG.error(this.getProgramName()+" : pw.checkError(): I/O error ###### ");
				return -1;
			}
			pw.close();pw=null;
			
					
			
			
			LOG.info("done");
			return RETURN_OK;
			} catch(final Exception err) {
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(pw);
				CloserUtil.close(in);
				CloserUtil.close(lr);
			}
		}
	
	public static void main(String[] args)
		{
		new VcfBurdenRscriptV().instanceMainWithExit(args);
		}
	}
