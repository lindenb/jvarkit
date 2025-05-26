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
package com.github.lindenb.jvarkit.tools.vcfsplitvep;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;

import htsjdk.samtools.util.RuntimeIOException;

/*
BEGIN_DOC

## Motivation

like `bcftools +split-vep` with aggregates




END_DOC
*/
@Program(name="vcfsplitvep",
	description="Split CSQ vep annotations",
	keywords={"vcf","vep","annotation"},
	creationDate="20250517",
	modificationDate="20250517",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VCFSplitVEP extends AbstractOnePassVcfAnnotator
	{
	private static final Logger LOG = Logger.of(VCFSplitVEP.class);

	@Parameter(names={"-t","--tag","--tags"},description="VEP tags name{:type{:aggregate}},name2{:type2{:aggregate2}},etc...  where type is a VCF info type Integer,String,Float an aggregate is one of none,min,max,uniq,first,random")
	private List<String> tags = new ArrayList<>();
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		try {
			final VCFSplitVEPAnnotator ann = new VCFSplitVEPAnnotator();
			ann.setTags(String.join(",",this.tags));
			return Collections.singletonList(ann);
			}
		catch(Throwable err) {
			throw new RuntimeIOException(err);
			}
		}
	

	public static void main(final String[] args)
		{
		new VCFSplitVEP().instanceMainWithExit(args);
		}

}
