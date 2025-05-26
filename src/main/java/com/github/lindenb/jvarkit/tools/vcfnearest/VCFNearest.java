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
package com.github.lindenb.jvarkit.tools.vcfnearest;

import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;

import htsjdk.samtools.util.RuntimeIOException;

/*
BEGIN_DOC

## Motivation

find nearest gene near a variant




END_DOC
*/
@Program(name="vcfnearest",
	description="find nearest feature near a variant",
	keywords={"vcf","vep","annotation"},
	creationDate="20250523",
	modificationDate="20250523",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VCFNearest extends AbstractOnePassVcfAnnotator
	{
	private static final Logger LOG = Logger.of(VCFNearest.class);

	@Parameter(names={"-t","--tag"},description="Tag name")
	private String tagName = "NEAREST";
	@Parameter(names={"-B","--bed"},description="Bed path. ",required = true)
	private Path bedPath = null;
	@Parameter(names={"-d","--distance"},description="max distance from variant to feature. "+DistanceParser.OPT_DESCRIPTION,splitter = NoSplitter.class,converter = DistanceParser.StringConverter.class)
	private int max_distance=0;
	@Parameter(names={"-n","--max-features"},description="do not print more than 'x' features; disable if x <=0")
	private int limit_n_features=-1;
	@Parameter(names={"-C","--columns"},description="comma separated names of the columns in the bed file. Empty columns will be skipped. e.g: ',,,NAME,TYPE,,SCORE' expects a bed with 7 columns. The 1,2,3 and 6th will be skipped. Two columns are added: the side (0: overlap, -1 bed is before variant , 1 bed is after variant ) and the distance",required = true)
	private String columnStr="";
	@Parameter(names={"--tabix"},description="use tabix index if possible (otherwise, put the bed in memory)")
	private boolean use_tabix;

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		try {
			final VCFNearestAnnotator ann = new VCFNearestAnnotator();
			ann.setTag(this.tagName);
			ann.setBedPath(this.bedPath);
			ann.setMaxDistance(max_distance);
			ann.setLimitCount(limit_n_features);
			ann.setUseTabixIndex(use_tabix);
			ann.setColumns(Arrays.asList(CharSplitter.COMMA.split(this.columnStr)));
			return Collections.singletonList(ann);
			}
		catch(Throwable err) {
			throw new RuntimeIOException(err);
			}
		}
	

	public static void main(final String[] args)
		{
		new VCFNearest().instanceMainWithExit(args);
		}

}
