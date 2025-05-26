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
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.vcf.AbstractOnePassVcfAnnotator;
import com.github.lindenb.jvarkit.wig.BigWigReader;

import htsjdk.samtools.util.RuntimeIOException;
/*
BEGIN_DOC

## Example

```bash
 java -jar dist/vcfbigwig.jar \
 	-T GERP \
 	-B gerp.bw input.vcf.gz 
	
##INFO=<ID=GERP,Number=1,Type=Float,Description="Values from bigwig file: com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig BIGWIG=gerp.bw TAG=GERP IN=input.vcf.gz    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO(...)
A	33926	.	G	A	182	.	GERP=-6.35(...)
A	45365	.	A	G	222	.	GERP=-3.55(...)
```


END_DOC
*/
@Program(name="vcfbigwig",
	description="Annotate a VCF with values from a bigwig file",
	keywords={"vcf","wig","wiggle","bigwig"},
	creationDate="20200506",
	modificationDate="20230819",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VCFBigWig extends AbstractOnePassVcfAnnotator {
	private static final Logger LOG = Logger.of(VCFBigWig.class);

	
	@Parameter(names={"-B","--bigwig"},description= BigWigReader.OPT_DESC,required=true)
	private String userBigWigFileUri = null;

	@Parameter(names={"-T","--tag","-tag"},description="Name of the INFO tag. default: name of the bigwig")
	private String userVcfTag = null;
		
	

	public VCFBigWig()
		{
		}
	@Override
	protected List<VariantAnnotator> createVariantAnnotators() {
		try {
			final BigWigVariantAnnotator ann= new BigWigVariantAnnotator(this.userBigWigFileUri);
			if(!StringUtils.isBlank(this.userVcfTag)) ann.setTag(userVcfTag);
			
			return Collections.singletonList(ann);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	public static void main(final String[] args) throws IOException
		{
		new VCFBigWig().instanceMainWithExit(args);
		}
}
