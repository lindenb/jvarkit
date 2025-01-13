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
package com.github.lindenb.jvarkit.tools.vcfpar;



import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.par.PseudoAutosomalRegion;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
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
	modificationDate="20200908",
	menu="VCF Manipulation"
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
	
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		
		final VCFHeader header= iterin.getHeader();
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		final PseudoAutosomalRegion parDefinition = PseudoAutosomalRegion.getInstance(dict).orElse(null);
		// https://en.wikipedia.org/wiki/Pseudoautosomal_region#Location
		if(parDefinition==null)
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
			switch(parDefinition.getLabel(ctx)) {
				case mixed:
				case sexual: {
					out.add(new VariantContextBuilder(ctx).attribute(info.getID(), Boolean.TRUE).make());
					break;
					}
				case autosomal:
				case pseudoautosomal:{
					if(autosome!=null) 
						{
						out.add(new VariantContextBuilder(ctx).attribute(autosome.getID(), Boolean.TRUE).make());
						}
					else
						{
						out.add(ctx);
						}
					break;
					}
				default: out.add(ctx); break;
				}
		}
		return 0;
		}
	


	public static void main(final String[] args)
		{
		new VcfPseudoAutosomalRegion().instanceMainWithExit(args);
		}
	
	}
