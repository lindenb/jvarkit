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

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;

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
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;
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
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	@XmlType(name="vcfmovefilterstoinfo")
	@XmlRootElement(name="vcfmovefilterstoinfo")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			@XmlElement(name="filter")
			@Parameter(names={"-f","--filter"},description="INFO name. This tag will be used to store the previous filters")
			private String infoName = "PREVIOUSLY_FILTERED_AS";
		
			@XmlElement(name="limit")
			@Parameter(names={"-t","--limitto"},description="If not empty, limit to those FILTERS. Multiple separated by comma/space.")
			private Set<String> onlyThoseFiltersTagStr = new HashSet<>();
	
			private class CtxWriter extends DelegateVariantContextWriter
				{
				private VCFInfoHeaderLine infoHeaderLine = null;
				private Set<String> limitToThoseFilters = null;
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				@Override
				public void writeHeader(final VCFHeader header) {
					this.limitToThoseFilters = 
							CtxWriterFactory.this.onlyThoseFiltersTagStr.stream().flatMap(
									S->Arrays.asList(S.split("[, ]")).stream()).
									filter(S->!StringUtil.isBlank(S)).
									collect(Collectors.toSet())
									;
					this.infoHeaderLine = new VCFInfoHeaderLine(
							CtxWriterFactory.this.infoName.trim(),
							VCFHeaderLineCount.UNBOUNDED,
							VCFHeaderLineType.String,
							"Variant was previously FILTERed with the given values."
							);
					if(header.getInfoHeaderLine(infoHeaderLine.getID())!=null)
						{
						throw new JvarkitException.UserError("INFO["+this.infoHeaderLine.getID()+"] already exists in input VCF.");
						}
					
					
					final VCFHeader h2= new VCFHeader(header);	
					h2.addMetaDataLine(this.infoHeaderLine);
					super.writeHeader(h2);
					}
				@Override
				public void add(final VariantContext ctx) {
					if(ctx.isNotFiltered())
						{
						super.add(ctx);
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
							if(!this.limitToThoseFilters.isEmpty() && !this.limitToThoseFilters.contains(filter)) {
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
							vcb.attribute(this.infoHeaderLine.getID(), new ArrayList<>(INFOfilters));
							}
						super.add(vcb.make());
						}
					}
				}
	
			@Override
			public int initialize() {
				if(StringUtil.isBlank(this.infoName)) {
					LOG.error("undefined option for infoName");
					return -1;
					}
				return 0;
				}
			
			@Override
			public VariantContextWriter open(VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			}
			
	public VcfMoveFiltersToInfo()
		{
		}
	 
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator in,
		final VariantContextWriter delegate)
		{
		final VariantContextWriter  out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(in.getHeader()).logger(LOG);
		out.writeHeader(in.getHeader());
		while(in.hasNext() &&  !out.checkError())
			{
			out.add(progess.watch(in.next()));
			}
		progess.finish();
		out.close();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, outputFile);
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	 	
	
	public static void main(final String[] args)
		{
		new VcfMoveFiltersToInfo().instanceMainWithExit(args);
		}
	}
