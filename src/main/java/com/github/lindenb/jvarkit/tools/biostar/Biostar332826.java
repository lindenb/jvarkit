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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
 BEGIN_DOC
 
 ## Example
 
 ```
 $ wget -O - -q "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" |\
 	gunzip -c |\
 	java -jar dist/biostar332826 --ids ids.txt > out.vcf 
 ```
 
 END_DOC
 */
@Program(name="biostar332826",
description="Fast Extraction of Variants from a list of IDs",
keywords= {"vcf","rs","id"},
modificationDate="20210412",
creationDate="20180817",
biostars={332826,433062}
)
public class Biostar332826 extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(Biostar332826.class).make();
	
	@Parameter(names={"-r","-i","--ids"},description="A list of identifiers, one per line")
	private Path rsFile = null;
	@Parameter(names={"-R","-I"},description="A semicolon/comma/space separated list of identifiers")
	private String rsStr = "";
	@Parameter(names={"-d","--delete"},description="When found , remove the ID from the list of identifiers unless it's a '.'. Should be faster but don't use it if two variants have the same ID.")
	private boolean removeIfFound=false;
	@Parameter(names={"--inverse"},description="Inverse: don't print the variants containing the IDS.")
	private boolean filterOutVariantsInSet =false;
	@Parameter(names={"-f","--filter"},description="if not blank soft filter the variants that are NOT in the list. "
			+ "If '--inverse' is specified then soft-filter the variants IN the list.")
	private String filterName = null;


	@Override
	protected Logger getLogger() {
		return LOG;
	}
	

	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final Set<String> rsSet = new HashSet<>(10_000);
		
		if(this.rsFile!=null)
			{
			try {
				rsSet.addAll(Files.readAllLines(this.rsFile));
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		
		for(final String str:this.rsStr.split("[ ;,]"))
			{
			if(StringUtil.isBlank(str)) continue;
			rsSet.add(str);
			}
		rsSet.remove("");
		if(rsSet.isEmpty()) getLogger().warn("NO IDENTIFIER WAS SPECIFIED");
		
		final Predicate<VariantContext> findID = (CTX)->{
			final String id = (!CTX.hasID()?VCFConstants.EMPTY_ID_FIELD:CTX.getID());
			if(removeIfFound && !VCFConstants.EMPTY_ID_FIELD.equals(id)) {
				return rsSet.remove(id);
				}
			else
				{
				return rsSet.contains(id);
				}
			};
			
		
		final VCFHeader header= iterin.getHeader();
		final VCFFilterHeaderLine filterHeader;
		if(StringUtils.isBlank(this.filterName)) {
			filterHeader= null;
		} else 
			{
			filterHeader = new VCFFilterHeaderLine(this.filterName, "Variant filtered for their IDs (N="+rsSet.size()+")");
			header.addMetaDataLine(filterHeader);
			}
		
		
		JVarkitVersion.getInstance().addMetaData(this, header);
		out.writeHeader(header);
		if(filterHeader!=null) {
			while(iterin.hasNext() && !rsSet.isEmpty()) {
				final VariantContext ctx = iterin.next();
				boolean keep=  findID.test(ctx);
				if(this.filterOutVariantsInSet) keep=!keep;
				if(keep) {
					out.add(ctx);
					}
				else
					{
					out.add(new VariantContextBuilder(ctx).filter(filterHeader.getID()).make());
					}
				}
			}
		else if(this.filterOutVariantsInSet)
			{
			while(iterin.hasNext() && !rsSet.isEmpty()) {
				final VariantContext ctx = iterin.next();
				if(!findID.test(ctx)) {
					out.add(ctx);
					}
				}
			//remaining variants
			while(iterin.hasNext()) out.add(iterin.next());
			}
		else
			{
			while(iterin.hasNext() && !rsSet.isEmpty()) {
				final VariantContext ctx = iterin.next();
				if(findID.test(ctx)) {
					out.add(ctx);
					}
				}
			}
		return 0;
		}

public static void main(final String[] args) throws IOException
	{
	new Biostar332826().instanceMainWithExit(args);
	}
}
