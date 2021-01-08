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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.nio.file.Path;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFRecordCodec;

/**

BEGIN_DOC

### Deprecated

Use picard SortVcf  http://broadinstitute.github.io/picard/command-line-overview.html#SortVcf.

Use `bcftools sort`

### Example

```
cat input.vcf |\
   java -jar dist/sortvcfonref2.jar  -R ref.fa |\
   bgzip -c > result.vcf.gz && \
   tabix -p vcf -f result.vcf.gz
```

END_DOC
*/

@Program(name="sortvcfonref2",
	description="Sort a VCF using the internal dictionary or an external reference order (Deprecated: use bcftools sort).",
	deprecatedMsg="use picard sortvcf",
	keywords={"vcf","sort"},
	modificationDate="20200518"
	)
public class SortVcfOnRef2 extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(SortVcfOnRef2.class).make();
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path refdict = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
   @Override
   protected int doVcfToVcf(final String inputName, final VCFIterator iterin, final VariantContextWriter w) {
    	CloseableIterator<VariantContext> iter=null;
    	SortingCollection<VariantContext> array=null;
    	try {
    		final VCFHeader header = iterin.getHeader();
    		final VCFHeader h2=new  VCFHeader(header);

		   final SAMSequenceDictionary dict;
		   if(this.refdict!=null) {
				dict= SequenceDictionaryUtils.extractRequired(this.refdict);
				h2.setSequenceDictionary(dict);
				}
		   else
		   		{
				dict= SequenceDictionaryUtils.extractRequired(header);
		   		}

    		array= SortingCollection.newInstance(
					VariantContext.class,
                    new VCFRecordCodec(h2,true),
                    h2.getVCFRecordComparator(),
                    this.writingSortingCollection.getMaxRecordsInRam(),
                    this.writingSortingCollection.getTmpPaths()
                    );
			array.setDestructiveIteration(true);
			while(iterin.hasNext())
				{
				final VariantContext cpl= iterin.next();
				array.add(cpl);
				}
			array.doneAdding();
			
			w.writeHeader(h2);
			iter=array.iterator();
			while(iter.hasNext())
				{
				w.add(iter.next());
				}
			iter.close();
			return 0;
			}
    	catch (final Throwable e)
    		{
			LOG.error(e);
			return -1;
			}
    	finally
	    	{
	    	if(array!=null) array.cleanup();
	    	}
    	
    	}

	public static void main(final String[] args) {
		new SortVcfOnRef2().instanceMainWithExit(args);
	}

}
