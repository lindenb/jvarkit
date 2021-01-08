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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/*
BEGIN_DOC

## Example

```
$ wget -O - -q "https://github.com/KennethJHan/Bioinformatics_Programming_101/raw/f40a5daa01cb9b37232d7d7a576100a10f867e80/GATK_BtPractice/SRR000982.filtered.variants.vcf" |  java -jar dist/vcfadfraction.jar | java -jar dist/vcf2table.jar | grep AD_RATIO -C 10

 Genotype Types
 +---------+-------+-----+
 | Type    | Count | %   |
 +---------+-------+-----+
 | HOM_VAR | 1     | 100 |
 +---------+-------+-----+
 Genotypes
 +-----------+---------+-----+----+----------+----+-----+-----+---------------+---------+
 | Sample    | Type    | AD  | DP | FT       | GQ | GT  | PGT | PID           | PL      |
 +-----------+---------+-----+----+----------+----+-----+-----+---------------+---------+
 | SRR000982 | HOM_VAR | 1,4 | 5  | AD_RATIO | 5  | 1/1 | 1|1 | 196625750_A_T | 163,5,0 |
 +-----------+---------+-----+----+----------+----+-----+-----+---------------+---------+
<<GRCh37 chr3:196625765/A (n. 90)
>>GRCh37 chr3:197119835/G (n. 91)
 Variant
 +-------+-----------+
 | Key   | Value     |
 +-------+-----------+
 | CHROM | chr3      |
 | POS   | 197119835 |
 | end   | 197119835 |
--
 Genotype Types
 +------+-------+-----+
 | Type | Count | %   |
 +------+-------+-----+
 | HET  | 1     | 100 |
 +------+-------+-----+
 Genotypes
 +-----------+------+------+----+----------+----+-----+----------+
 | Sample    | Type | AD   | DP | FT       | GQ | GT  | PL       |
 +-----------+------+------+----+----------+----+-----+----------+
 | SRR000982 | HET  | 3,16 | 19 | AD_RATIO | 56 | 0/1 | 472,0,56 |
 +-----------+------+------+----+----------+----+-----+----------+
<<GRCh37 chr10:42385236/A (n. 167)
>>GRCh37 chr10:42385255/G (n. 168)
```


END_DOC
 */
@Program(name="vcfadfraction",
	description="filter VCF for strange FORMAT/AD fraction",
	keywords= {"vcf","allele-depth","AD","depth"},
	creationDate="20190723",
	modificationDate="20190726"
	)
public class VcfAlleleDepthFraction extends Launcher {
	private static final Logger LOG = Logger.build(VcfAlleleDepthFraction.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null; 
	@Parameter(names={"-f","--filtered"},description="ignore FILTER-ed **GENOTYPES**")
	private boolean ignore_filtered_genotypes = false;
	@Parameter(names={"-het","--het"},description="AD ratio for **HET** genotypes. HET genotype should have x <= AD[1]/(AD[0]+AD[1])<= (1-x)")
	private double het_treshold = 0.2;
	@Parameter(names={"-hom","--hom"},description="AD ratio for **HOM_REF** or **HOM_VAR** genotypes. "
			+ "HOM_REF genotype should have x <= AD[1]/(AD[0]+AD[1]). "
			+ "HOM_VAR genotype should have  AD[1]/(AD[0]+AD[1]) >= (1-x). "
			)
	private double hom_treshold = 0.05;
	@Parameter(names={"-gtf","--gtf"},description="Genotype FILTER")
	private String genotype_filter = "AD_RATIO";
	@Parameter(names={"-filter","--filter"},description="Variant FILTER")
	private String variant_filter = "AD_RATIO";
	@Parameter(names={"-maxFilteredGenotypes","--maxFilteredGenotypes"},description="Set Variant FILTER if number of BAD genotype is greater than 'x'. Negative is ignore.")
	private int maxFilteredGenotypes=-1;
	@Parameter(names={"-maxFractionFilteredGenotypes","--maxFractionFilteredGenotypes"},description="Set Variant FILTER if percent of BAD genotype is greater than 'x'. Negative is ignore.")
	private double maxFractionFilteredGenotypes=-1;
	@Parameter(names={"-dp","--dp"},description="Only consider Genotypes having DP> 'x'")
	private int min_depth = -1;

		
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator iterin,
		final VariantContextWriter out
		) {

		final VCFHeader header = iterin.getHeader();
		
		final VCFHeader header2 = new VCFHeader(header);
		if(header.hasGenotypingData()) {
			header2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
			}
		if(!StringUtils.isBlank(this.variant_filter)) {
			header2.addMetaDataLine(new VCFFilterHeaderLine(this.variant_filter, "Variant filter for AD. HET:"+this.het_treshold+" HOM_REF/HOM_VAR:"+this.hom_treshold));
			}
		JVarkitVersion.getInstance().addMetaData(VcfAlleleDepthFraction.class.getSimpleName(), header2);
		out.writeHeader(header2);
		
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.
				newInstance().
				logger(LOG).
				dictionary(header).
				build();
		while(iterin.hasNext())
			{
			final VariantContext ctx = progress.apply(iterin.next());
			if(!ctx.hasGenotypes() || !ctx.isVariant())
				{
				out.add(ctx);
				continue;
				}
			
			int count_bad_gt=0;
			final List<Genotype> newgt = new ArrayList<>(ctx.getNSamples());
			for(final Genotype gt: ctx.getGenotypes())
				{
				if(!gt.hasAD())
					{
					newgt.add(gt);
					continue;
					}
				
				if(this.min_depth>=0 && (!gt.hasDP() || gt.getDP() < this.min_depth)) {
					newgt.add(gt);
					continue;
					}
				
				final int ad[]=gt.getAD();
				if(ad==null  || ad.length==1 || ad.length!=ctx.getNAlleles()) {
					newgt.add(gt);
					continue;
					}
				
				
				
				boolean ok_ad=true;
				
				if(gt.isHet())
					{
					final int altidx = gt.getAlleles().stream().filter(A->!A.isReference()).mapToInt(A->ctx.getAlleleIndex(A)).findFirst().orElse(-1);
					if(altidx>0 && altidx < ad.length && ad[0]+ad[altidx]>0)
						{
						final double ratio = ad[altidx]/(double)(ad[0]+ad[altidx]);
						if(ratio< this.het_treshold || ratio >(1-this.het_treshold)) {
							ok_ad = false;
							}
						}
					}
				else if(gt.isHomRef()) {
					final int ad_not_ref= Arrays.stream(ad).skip(1L).sum();
					if(ad_not_ref>0 && ad[0]+ ad_not_ref >0)
						{
						final double ratio = ad_not_ref/(double)(ad[0]+ad_not_ref);
						if(ratio > this.hom_treshold) {
							ok_ad = false;
							}
						}
					}
				else if(gt.isHomVar())
					{
					final int altidx = gt.getAlleles().stream().filter(A->!A.isReference()).mapToInt(A->ctx.getAlleleIndex(A)).findFirst().orElse(-1);
					if(altidx>0 && altidx < ad.length && ad[0]+ad[altidx]>0)
						{
						final double ratio = ad[altidx]/(double)(ad[0]+ad[altidx]);
						if(ratio < (1.0-this.hom_treshold)) {
							ok_ad = false;
							}
						}
					}
				if(ok_ad)
					{
					newgt.add(gt);
					}
				else
					{
					if(StringUtils.isBlank(this.genotype_filter)) {
						newgt.add(gt);
						}
					else
						{
						newgt.add(new GenotypeBuilder(gt).filter(this.genotype_filter).make());
						}
					count_bad_gt++;
					}
				}
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.genotypes(newgt);
			boolean set_filter=false;
			
			
			if(this.maxFilteredGenotypes>=0 && count_bad_gt> this.maxFilteredGenotypes) {
				set_filter = true;
				}
			if(this.maxFractionFilteredGenotypes>=0 && (count_bad_gt/(double)ctx.getNSamples())> this.maxFractionFilteredGenotypes) {
				set_filter = true;
				}
			
			if(!StringUtils.isBlank(this.variant_filter)) 
				{
				//nothing
				}
			else if(set_filter)
				{
				final Set<String> filters = ctx.getFilters();
				filters.add(this.variant_filter);
				vcb.filters(filters);
				}
			else if(!ctx.isFiltered())
				{
				vcb.passFilters();
				}
			
			out.add(vcb.make());
			}
		progress.close();
		return 0;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.het_treshold<0 || this.het_treshold>1.0) {
			LOG.error("bad het treshold "+ this.het_treshold);
			return -1;
			}
		if(this.het_treshold>0.5) this.het_treshold=1-this.het_treshold;
		
		if(this.hom_treshold<0 || this.hom_treshold>0.5) {
			LOG.error("bad hom treshold "+ this.hom_treshold);
			return -1;
			}		
		return doVcfToVcf(args, this.outputFile);
		}
	public static void main(final String[] args) {
		new VcfAlleleDepthFraction().instanceMainWithExit(args);
	}

}
