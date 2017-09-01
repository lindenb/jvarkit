/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/**

BEGIN_DOC

For Matilde K: move the information in FILTER to the INFO column to keep a trace of the FILTERs.

## Example

```
$ cat input.vcf | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16057607	rs201535778	G	GAAAA	.	.	AAC=2;AAF=0.0005;BEACON=T|solvebio,T|bob,T|solvebio-133,T|altruist,T|prism,T|kaviar;NS=3690;RAC=3688;RAF=0.9995;VTYPE=SNV
22	16057608	rs201535778	G	T	.	.	AAC=2;AAF=0.0005;BEACON=T|solvebio,T|bob,T|solvebio-133,T|altruist,T|prism,T|kaviar;NS=3690;RAC=3688;RAF=0.9995;VTYPE=SNV


$ java -jar dist/vcfmovefilterstoinfo.jar input.vcf | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16057608	rs201535778	G	T	.	.	AAC=2;AAF=0.0005;NS=3690;PREVIOUSLY_FILTERED_AS=FT1;RAC=3688;RAF=0.9995;VTYPE=SNV
22	16058492	.	G	A	.	.	AAC=2;AAF=0.0005;NS=3708;PREVIOUSLY_FILTERED_AS=FT2;RAC=3706;RAF=0.9995;VTYPE=SNV


$ java -jar dist/vcfmovefilterstoinfo.jar -f OLDFILTER -t FT2 input.vcf | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16057608	rs201535778	G	T	.	FT1	AAC=2;AAF=0.0005;NS=3690;RAC=3688;RAF=0.9995;VTYPE=SNV
22	16058492	.	G	A	.	.	AAC=2;AAF=0.0005;NS=3708;OLDFILTER=FT2;RAC=3706;RAF=0.9995;VTYPE=SNV

```

END_DOC
*/

@Program(name="vcfmovefilterstoinfo",
		description="Move any FILTER to the INFO column. reset FILTER to PASS",
		keywords={"vcf","burden","format","info"}
		)
public class VcfMoveFiltersToInfo
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfMoveFiltersToInfo.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-f","--filter"},description="INFO name. This tag will be used to store the previous filters")
	private String infoName = "PREVIOUSLY_FILTERED_AS";

	@Parameter(names={"-t","--limitto"},description="If not empty, limit to those FILTERS. Multiple separated by comma/space.")
	private Set<String> onlyThoseFiltersTagStr = new HashSet<>();
	
	
	public VcfMoveFiltersToInfo()
		{
		}
	 
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out)
		{
		if(StringUtil.isBlank(this.infoName)) {
			LOG.error("undefined option for infoName");
			return -1;
			}
		final VCFHeader header = in.getHeader();
		final Set<String> limitToThoseFilters = 
		this.onlyThoseFiltersTagStr.stream().flatMap(
				S->Arrays.asList(S.split("[, ]")).stream()).
				filter(S->!StringUtil.isBlank(S)).
				collect(Collectors.toSet())
				;
		
		try {
			final VCFInfoHeaderLine infoHeaderLine = new VCFInfoHeaderLine(
					this.infoName.trim(),
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Variant was previously FILTERed with the given values."
					);
			if(header.getInfoHeaderLine(infoHeaderLine.getID())!=null)
				{
				LOG.error("INFO["+infoHeaderLine.getID()+"] already exists in input VCF.");
				return -1;
				}
			
			
			final VCFHeader h2= new VCFHeader(header);	
			h2.addMetaDataLine(infoHeaderLine);
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary()).logger(LOG);
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.watch(in.next());
				if(ctx.isNotFiltered())
					{
					out.add(ctx);
					}
				else
					{
					final Set<String> INFOfilters = new HashSet<>();
					final Set<String> FILTERfilters = new HashSet<>();
					for(final String filter : ctx.getFilters()) {
						if( filter.equals(VCFConstants.UNFILTERED) ||
							filter.equals(VCFConstants.PASSES_FILTERS_v3) ||
							filter.equals(VCFConstants.PASSES_FILTERS_v4)
							)
							{
							continue;
							}
						if(!limitToThoseFilters.isEmpty() && !limitToThoseFilters.contains(filter)) {
							FILTERfilters.add(filter);
							}
						else
							{
							INFOfilters.add(filter);
							}
						}
					
					
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					if(FILTERfilters.isEmpty()) {
						vcb.unfiltered();
						}
					else
						{
						vcb.filters(FILTERfilters);
						}
					if(!INFOfilters.isEmpty()) {
						vcb.attribute(infoHeaderLine.getID(), new ArrayList<>(INFOfilters));
						}
					out.add(vcb.make());
					}
				}
			progess.finish();
			LOG.info("done.");
			return RETURN_OK;
			} catch(final Exception err) {
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(in);
				CloserUtil.close(out);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args,outputFile);
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfMoveFiltersToInfo().instanceMainWithExit(args);
		}
	}
