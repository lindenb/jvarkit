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
package com.github.lindenb.jvarkit.tools.vcfpar;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/*
BEGIN_DOC

 
## Example


```bash
java -jar dist/vcfpar.jar   in.vcf > out.vcf
```

END_DOC
 */
@Program(name="vcfpar",
	description="Flag human sexual regions excluding PAR.",
	keywords={"vcf","sex","par"},
	creationDate="20200908",
	modificationDate="20200908"
	)
public class VcfPseudoAutosomalRegion extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfPseudoAutosomalRegion.class).make();

	@Parameter(names={"--tag"},description="VCF info TAG")
	private String tag = "SEX";
	@Parameter(names={"--autosome"},description="if not empty, autosomal and PAR regions will be annotated with this flag. Ignoring the chromosome looking like mitochondrial chromosome.")
	private String autosome_tag = null;

	public VcfPseudoAutosomalRegion() {
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private static class ParRegionDefition
		{
		private final String name;
		private final int[] Xrgn;
		private final int[] Yrgn;
		
		ParRegionDefition(final String name,int[] Xrgn, int[] Yrgn) {
			if(Xrgn.length%2!=0 || Xrgn.length != Yrgn.length) throw new IllegalArgumentException();
			this.name = name;
			this.Xrgn = Xrgn;
			this.Yrgn = Yrgn;
		}
		
		private boolean isSex(final int[] array,final int start,final int end) {
			for(int i=0;i+1 <array.length;i+=2) {
				if( CoordMath.overlaps(array[i], array[i+1], start, end)) return false;
				}
			return true;
			}
		
		public boolean isSex(final Locatable l) {
			final String contig = l.getContig();
			if(contig.equals("chrX") || contig.equals("X")) {
				return isSex(this.Xrgn,l.getStart(),l.getEnd());
				}
			else if(contig.equals("chrY") || contig.equals("Y")) {
				return isSex(this.Yrgn,l.getStart(),l.getEnd());
				}
			else if(contig.startsWith("chrX_") || contig.startsWith("chrY_")) {
				return true;
				}
			return false;
			}
		private String intervals(final String ctg,int[] array) {
			String s="";
			for(int i=0;i+1 <array.length;i+=2) {
				s+=ctg+":"+array[i]+"-"+array[i+1]+" ";
				}
			return s;
			}
		String getDescription() {
			return "Variant is mapped on a sexual chromosome on "+getName()+" excluding PAR region(s) on "+intervals("X",Xrgn)+intervals("Y",Yrgn); 
			}
		public String getName() {
			return this.name;
			}
		}
	
	private boolean isMitochondrial(final String s) {
		return s.equals("M") || s.equals("MT") || s.equals("chrM") || s.equals("chrMT");
	}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final ParRegionDefition parDefinition;
		final VCFHeader header= iterin.getHeader();
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		// https://en.wikipedia.org/wiki/Pseudoautosomal_region#Location
		if(SequenceDictionaryUtils.isGRCh38(dict)) {
			parDefinition = new ParRegionDefition(
				"GRCh38",
				new int[] {10_001,2_781_479, 155_701_383,156_030_895},
				new int[] {10_001,2_781_479, 56_887_903, 57_217_415}
				);
			}
		else if(SequenceDictionaryUtils.isGRCh37(dict)) {
			parDefinition = new ParRegionDefition(
				"GRCh37",
				new int[] {60_001,2_699_520, 154_931_044,155_260_560 },
				new int[] {10_001,2_649_520, 59_034_050,59_363_566}
				);
			}
		else
			{
			LOG.error("Vcf dictionary doesn't look like a known human reference.");
			return -1;
			}
		final VCFInfoHeaderLine info = new VCFInfoHeaderLine(
				this.tag,1,
				VCFHeaderLineType.Flag,
				parDefinition.getDescription()
				);
		header.addMetaDataLine(info);
		
		
		final VCFInfoHeaderLine autosome;
        
		if(!StringUtils.isBlank(this.autosome_tag)) {
			if(this.tag.equals(this.autosome_tag)) {
				LOG.error("autosome INFO/TAG is same as sexual INFO/TAG.");
				return -1;
				}
			autosome = new VCFInfoHeaderLine(
					this.autosome_tag,1,
					VCFHeaderLineType.Flag,
					"Variant is on autosomal chromosome or PAR."
					);
			header.addMetaDataLine(autosome);
		} else {
			autosome = null;
		}
		
		
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		while(iterin.hasNext()) {
			final VariantContext ctx = iterin.next();
			if(parDefinition.isSex(ctx)) {
				out.add(new VariantContextBuilder(ctx).attribute(info.getID(), Boolean.TRUE).make());
				}
			else if(autosome!=null && !isMitochondrial(ctx.getContig())) 
				{
				out.add(new VariantContextBuilder(ctx).attribute(autosome.getID(), Boolean.TRUE).make());
				}
			else
				{
				out.add(ctx);
				}
			}
		return 0;
		}
	


	public static void main(final String[] args)
		{
		new VcfPseudoAutosomalRegion().instanceMainWithExit(args);
		}
	
	}
