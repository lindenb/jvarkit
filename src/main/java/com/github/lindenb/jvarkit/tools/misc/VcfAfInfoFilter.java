/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Motivation

I'm often asked to filter out variant that are too frequent in gnomad, but I must keep the data if any ALT allele is NOT in gnomad.

This tool filters VCF containing external allele frequency information (AF or AC/AN). Used as a  complement of VcfGnomad.

## Example

```
$ java -jar dist/vcfafinfofilter.jar -nfe input.vcf
$ java -jar dist/vcfafinfofilter.jar -af 'gnomad_exome_AF_NFE,gnomad_genome_AF_NFE'   input.vcf
$ java -jar dist/vcfafinfofilter.jar -acn 'gnomad_genome_AC_NFE,gnomad_genome_AN_NFE'   input.vcf
$ java -jar dist/vcfafinfofilter.jar -acn 'gnomad_genome_*_NFE'   input.vcf

```



END_DOC
 */
@Program(name="vcfafinfofilter",
	description="Filter VCF annotated with external (AF or AC/AN) frequency information like vcfgnomad",
	keywords={"vcf","annotation","gnomad"})
public class VcfAfInfoFilter extends Launcher{
	
	private static final Logger LOG = Logger.build(VcfAfInfoFilter.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"--filter","-f"},description="set this filter if all ALT fails the treshold. If empty :remove the variant")
	private String filterAllAltInGnomad="";
	@Parameter(names={"--gtfilter","-gtf"},description="set this *GENOTYPE* filter if all ALT for a Genotype fail the treshold. If empty :set genotype to NO_CALL")
	private String genotypeFilter="";
	@Parameter(names={"--treshold","-t"},description="Treshold for allele Frequency. ALT alleles above this AF value will be subject to filtration.")
	private double user_af_treshold = 1E-3;
	@Parameter(names={"-af","--af"},description="A list of AF fields, separated with comma,semicolon or whitespace that will be used to extract a AF field.")
	private String user_af_fields = "";
	@Parameter(names={"-acn","--acn"},description="A list of pairs of AC/AF fields separated with comma,semicolon or whitespace that will be used to calculate the AF. "
			+ "If an attribute contains  '*', it will be replaced by 'AC' and 'AN'. eg: 'gnomad_exome_AC_NFE,gnomad_exome_AN_NFE,my_pop_*'")
	private String user_ac_an_fields = "";
	@Parameter(names={"-i","--no-valid"},description="Ignore INFO Field Validation. (e.g INFO field not declarated in VCF header)")
	private boolean ignore_INFO_field_validation=false;
	@Parameter(names={"-nfe","--nfe"},description="Add INFO fields for the 'NFE' population created by vcfgnomad: gnomad_exome_AC_NFE,gnomad_exome_AF_NFE,gnomad_exome_AN_NFE,gnomad_genome_AC_NFE,gnomad_genome_AF_NFE,gnomad_genome_AN_NF")
	private boolean vcf_gnomad_nfe =false;
	@ParametersDelegate
	private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();

	
	private abstract class AFExtractor
		{
		abstract boolean validateHeader(final VCFHeader header);
		abstract List<Object> parse(final VariantContext ctx);
		}
	
	private class AFFieldExtractor extends AFExtractor
		{
		private String afAttr;
		AFFieldExtractor(final String afAttr) {
			this.afAttr = afAttr;
			if(StringUtil.isBlank(this.afAttr)) throw new IllegalArgumentException("AF is blank");
			}
		@Override
		boolean validateHeader(final VCFHeader header) {
			final VCFInfoHeaderLine hdr = header.getInfoHeaderLine(this.afAttr);
			if(hdr==null) {
				final String msg = "VCF header doesn't contain the INFO field " +this.afAttr;
				LOG.warn(msg);
				return false;
				}
			if(hdr.getType()!=VCFHeaderLineType.Float) {
				final String msg = "INFO=<ID=" +this.afAttr+"> extected Type=Float but got "+hdr.getType();
				LOG.warn(msg);
				return false;
				}
			if(hdr.getCountType()!=VCFHeaderLineCount.A) {
				final String msg = "INFO=<ID=" +this.afAttr+"> extected Number=A but got "+hdr.getCountType();
				LOG.warn(msg);
				return false;
				}
			return true;
			}
		@Override
		List<Object> parse(final VariantContext ctx) {
			if(!ctx.hasAttribute(this.afAttr)) return null;
			return ctx.getAttributeAsList(this.afAttr);
			}
		@Override
		public String toString() {
			return "AF Extractor using INFO/"+this.afAttr;
			}
		}
	private class ACANFieldsExtractor extends AFExtractor
		{
		private String acAttr;
		private String anAttr;
		ACANFieldsExtractor(final String acAttr,final String anAttr) {
			this.acAttr = acAttr;
			this.anAttr = anAttr;
			if(StringUtil.isBlank(this.acAttr)) throw new IllegalArgumentException("AC is blank");
			if(StringUtil.isBlank(this.anAttr)) throw new IllegalArgumentException("AN is blank");
			if(this.acAttr.equals(this.anAttr)) throw new IllegalArgumentException(acAttr+"="+this.anAttr);
			}
		@Override
		boolean validateHeader(final VCFHeader header) {
			VCFInfoHeaderLine hdr = header.getInfoHeaderLine(this.acAttr);
			if(hdr==null) {
				final String msg = "VCF header doesn't contain the INFO field " +this.acAttr;
				LOG.warn(msg);
				return false;
				}
			if(hdr.getType()!=VCFHeaderLineType.Integer) {
				final String msg = "INFO=<ID=" +this.acAttr+"> expected Type=Integer but got "+hdr.getType();
				LOG.warn(msg);
				return false;
				}
			if(hdr.getCountType()!=VCFHeaderLineCount.A) {
				final String msg = "INFO=<ID=" +this.acAttr+"> expected Number=A but got "+hdr.getCountType();
				LOG.warn(msg);
				return false;
				}
			hdr = header.getInfoHeaderLine(this.anAttr);
			if(hdr==null) {
				final String msg = "VCF header doesn't contain the INFO field " +this.anAttr;
				LOG.warn(msg);
				return false;
				}
			if(hdr.getType()!=VCFHeaderLineType.Integer) {
				final String msg = "INFO=<ID=" +this.anAttr+"> expected Type=Integer but got "+hdr.getType();
				LOG.warn(msg);
				return false;
				}
			
			if(!hdr.isFixedCount())
				{
				final String msg = "INFO=<ID=" +this.anAttr+"> expected fixed count "+hdr.getCountType();
				LOG.warn(msg);
				return false;
				}
			
			if(hdr.getCount()!=1) {
				final String msg = "INFO=<ID=" +this.anAttr+"> expected Number=1 but got "+hdr.getCount();
				LOG.warn(msg);
				return false;
				}
			return true;
			}
		
		@Override
		List<Object> parse(final VariantContext ctx) {
			if(!ctx.hasAttribute(this.acAttr)) return null;
			if(!ctx.hasAttribute(this.anAttr)) return null;
			final double an = (double)ctx.getAttributeAsInt(this.anAttr, 0);
			if(an>0) {
				return ctx.getAttributeAsList(this.acAttr).
						stream().
						map(O->objectToFrequency(O)).
						map(V->V==null?null:V/an).
						collect(Collectors.toList());
				}
			return null;
			}
		@Override
		public String toString() {
			return "AF Extractor using INFO/"+this.acAttr+" and INFO/"+this.anAttr;
			}
		}

	public VcfAfInfoFilter() {
		
	}
	
	private Double objectToFrequency(final Object o) {
		if(o==null) return null;
		if(o instanceof Double) return Double.class.cast(o);
		if(o instanceof Float) return Float.class.cast(o).doubleValue();
		if(o instanceof Integer) return Integer.class.cast(o).doubleValue();
		final String s = o.toString();
		if(StringUtil.isBlank(s) || s.equals("null") || s.equals(".")) return null;
		return Double.parseDouble(s);
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final  VcfIterator in,
			final  VariantContextWriter out
			) {
		try
			{
			if(this.vcf_gnomad_nfe) {
				this.user_af_fields += ",gnomad_exome_AF_NFE,gnomad_genome_AF_NFE";
				this.user_ac_an_fields += ",gnomad_exome_AC_NFE,gnomad_exome_AN_NFE,"
						+ "gnomad_genome_AC_NFE,gnomad_genome_AN_NFE";
			
			}
			
			
			final List<AFExtractor> afExtractors = new ArrayList<>();
			final VCFHeader header = in.getHeader();
			
			if(!StringUtil.isBlank(this.user_af_fields)) {
				for(final String s:this.user_af_fields.split("[,; \t\n]+"))
					{
					if(StringUtil.isBlank(s)) continue;
					if(s.contains("*")) {
						LOG.error("INFO/"+s+" contains a '*'");
						return -1;
						}
					final AFFieldExtractor extractor = new AFFieldExtractor(s);
					if(!extractor.validateHeader(header) && !this.ignore_INFO_field_validation)
						{
						LOG.warn("Ignoring "+extractor);
						continue;
						}
					afExtractors.add(extractor);
					}
				}
			if(!StringUtil.isBlank(this.user_ac_an_fields)) {
				final String array[]=this.user_ac_an_fields.split("[,; \t\n]+");
				int i=0;
				while(i<array.length)
					{
					String s= array[i];
					if(StringUtil.isBlank(s)) {i++; continue;}
					final String acf;
					final String anf;
					
					if(s.contains("*"))
						{
						acf = s.replace("*", "AC");
						anf = s.replace("*", "AN");
						}
					else
						{
						acf = s;
						i++;
						while(i<array.length)
							{
							s= array[i];
							if(StringUtil.isBlank(s)) {i++; continue;}
							break;
							}
						if(i==array.length) {
							LOG.error("missing AN for "+acf+ " in "+this.user_ac_an_fields);
							return -1;
							}
						anf = s;
						}
					
					
					final ACANFieldsExtractor extractor = new ACANFieldsExtractor(acf,anf);
					if(!extractor.validateHeader(header) && !this.ignore_INFO_field_validation)
						{
						LOG.warn("Ignoring "+extractor);
						i++;
						continue;
						}
					afExtractors.add(extractor);
					i++;
					}
				}
			
			if(afExtractors.isEmpty()) {
				LOG.warn("No extractor was defined !");
			}
			
			final Set<VCFHeaderLine> headerLines = new HashSet<>();
			final VCFFilterHeaderLine noAltVariantFilter = new VCFFilterHeaderLine(
					filterAllAltInGnomad,
					"All ALT alleles don't pass the gnomad treshold AF < "+this.user_af_treshold
					);
			if(!StringUtil.isBlank(this.filterAllAltInGnomad))
				{
				headerLines.add(noAltVariantFilter);
				}
			VCFStandardHeaderLines.addStandardFormatLines(headerLines, true, VCFConstants.GENOTYPE_FILTER_KEY);
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header);
			this.recalculator.setHeader(header);
			headerLines.stream().forEach(H->header.addMetaDataLine(H));
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			out.writeHeader(header);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.watch(in.next());

				if(!ctx.isVariant())
					{
					out.add(ctx);
					continue;
					}
				final List<Allele> alt_alleles = ctx.getAlternateAlleles();
				final Set<Allele> ok_alleles = new HashSet<>(alt_alleles);
				ok_alleles.remove(Allele.SPAN_DEL);
				for(final AFExtractor afExtractor : afExtractors)
					{
					if(ok_alleles.isEmpty()) break;
					final List<Object> afo_list = afExtractor.parse(ctx);
					if(afo_list==null) continue;
					if(afo_list.size()!=alt_alleles.size()) {
						LOG.warn("in "+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference()+" illegal number of AF values "+afExtractor);
						}
					for(int x=0;x< afo_list.size() && x < alt_alleles.size();++x)
						{
						final Double af = objectToFrequency(afo_list.get(x));
						if(af==null) continue;
						if(af.doubleValue()> this.user_af_treshold)
							{
							
							ok_alleles.remove(alt_alleles.get(x));
							}
						}
					}

				if(ok_alleles.isEmpty() )
					{
					if(!StringUtil.isBlank(this.filterAllAltInGnomad))
						{
						out.add(new VariantContextBuilder(ctx).filter(this.filterAllAltInGnomad).make());
						}
					continue;
					}
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				if(ok_alleles.isEmpty())
					{
					vcb.filter(noAltVariantFilter.getID());
					}
				final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
				for(final Genotype gt:ctx.getGenotypes()) {
					if(gt.isNoCall() || gt.isHomRef())
						{
						genotypes.add(gt);
						continue;
						}
					boolean got_good_allele = false;
					boolean got_bad_alt=false;
					for(final Allele gta :gt.getAlleles())
						{
						if(!gta.isCalled() || gta.equals(Allele.SPAN_DEL)) {
							continue;
							}
						else if(gta.isReference() || ok_alleles.contains(gta)) {
							got_good_allele = true;
							}
						else if(alt_alleles.contains(gta)) {
							got_bad_alt = true;
							}
						}
					if(got_good_allele || !got_bad_alt) {
						genotypes.add(gt);
						continue;
						}
					
					if(StringUtil.isBlank(this.genotypeFilter)) {
						genotypes.add(GenotypeBuilder.createMissing(gt.getSampleName(),gt.getPloidy()));
						}
					else
						{
						genotypes.add(new GenotypeBuilder(gt).filter(this.genotypeFilter).make());
						}
					}
				vcb.genotypes(genotypes);
				out.add(this.recalculator.apply(vcb.make()));
				}
			progress.finish();
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		}

public static void main(final String[] args) {
	new VcfAfInfoFilter().instanceMainWithExit(args);
	}
}
