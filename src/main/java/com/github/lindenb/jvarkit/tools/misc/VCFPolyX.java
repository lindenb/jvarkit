/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.nio.file.Path;
import java.util.List;


import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.ChromosomeSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.vcf.VCFIterator;

/*
BEGIN_DOC

## Example

```
$ java  -jar dist/vcfpolyx.jar -R reference.fa input.vcf
(...)
2   1133956 .   A   G   2468.84 .   POLYX=23
2   1133956 .   A   AG  3604.25 .   POLYX=23
2   2981671 .   T   G   47.18   .   POLYX=24
(...)
```

END_DOC
*/
@Program(name="vcfpolyx",
	description="Number of repeated REF bases around POS.",
	keywords={"vcf","repeat"},
	modificationDate="20191129"
	)
public class VCFPolyX extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFPolyX.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-n","--filter"},description="if number of repeated bases is greater or equal to 'n' set a FILTER = (tag)")
	private int filterTrehsold = -1 ;
	@Parameter(names={"-t","--tag"},description="Tag used in INFO and FILTER columns.")
	private String polyXtag = "POLYX";
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faixPath = null;
	@Parameter(names={"--skip-filtered"},description="Don't spend some time to calculate the tag if the variant is FILTERed")
	private boolean skip_filtered=false;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator r,
			final VariantContextWriter w
			) 
		{
		ReferenceSequenceFile referenceSequenceFile = null;
		
		if(StringUtil.isBlank(this.polyXtag)) {
			LOG.error("Empty tag");
			return -1;
			}
		
		try
			{
			referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faixPath);
			
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(referenceSequenceFile));

			
			final VCFHeader h2 = new VCFHeader(r.getHeader());
	
			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
					this.polyXtag.trim(),
					1,
					VCFHeaderLineType.Integer,
					"Number of repeated bases around REF")
					;
			h2.addMetaDataLine(infoHeaderLine);

			final VCFFilterHeaderLine filterHeaderLine = new VCFFilterHeaderLine(
					infoHeaderLine.getID()+"_ge_"+this.filterTrehsold,
					"Number of repeated bases around REF is greater or equal to " +
							this.filterTrehsold)
					;
			
			if( this.filterTrehsold>-1) {
				h2.addMetaDataLine(filterHeaderLine);
				}
			ChromosomeSequence genomicContig=null;
			JVarkitVersion.getInstance().addMetaData(this, h2);
			w.writeHeader(h2);
			final ProgressFactory.Watcher<VariantContext> progress= ProgressFactory.newInstance().dictionary(r.getHeader()).logger(LOG).build();
			while(r.hasNext())
				{
				final VariantContext ctx = progress.apply(r.next());
				if(this.skip_filtered && ctx.isFiltered())
					{
					w.add(ctx);
					continue;
					}
				
				final String normalizedContig = contigNameConverter.apply(ctx.getContig());
				if(normalizedContig==null) {
					w.add(ctx);
					continue;
					}
				
				if(genomicContig==null || !genomicContig.hasName(normalizedContig))
					{
					genomicContig= new GenomicSequence(referenceSequenceFile, normalizedContig);
					}
				
				final VariantContextBuilder b = new VariantContextBuilder(ctx);

				
				int pos0=ctx.getStart()-1;
				int count=1;
				char c0=Character.toUpperCase(genomicContig.charAt(pos0));
				//go left
				pos0--;
				while(pos0>=0 && c0==Character.toUpperCase(genomicContig.charAt(pos0)))
					{
					++count;
					pos0--;
					}
				//go right
				pos0=ctx.getEnd()-1;
				c0=Character.toUpperCase(genomicContig.charAt(pos0));
				pos0++;
				while(pos0< genomicContig.length()
					&& c0==Character.toUpperCase(genomicContig.charAt(pos0)))
					{
					++count;
					++pos0;
					}
				b.attribute(infoHeaderLine.getID(),count);
				
				/* filter */
				if(this.filterTrehsold>-1 && count>=this.filterTrehsold)
					{
					b.filter(filterHeaderLine.getID());
					}
				
				w.add(b.make());			
				}
			progress.close();
			w.close();
			}
		finally
			{
			CloserUtil.close(referenceSequenceFile);
			}
		return RETURN_OK;
		}
	
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcfPath(args,this.writingVariants,this.outputFile);
		}

	public static void main(final String[] args)
		{
		new VCFPolyX().instanceMainWithExit(args);
		}

}
