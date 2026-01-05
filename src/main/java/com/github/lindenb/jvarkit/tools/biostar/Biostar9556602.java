/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/** 

BEGIN_DOC

## input

input must be sorted on coordinate, Alleles must be normalized (one ALT per variant)

##Example

```bash


$ grep -v "##" jeter.vcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905593	.	A	T	3887.15	PASS	AF=0.1
1	905594	rs867360694	A	T	3887.15	PASS	AF=0.1
1	905594	rs867360694	A	TT	3887.15	PASS	AF=0.2
1	905594	rs867360694	A	TTT	3887.15	PASS	AF=0.3
1	905595	.	A	T	3887.15	PASS	AF=0.1
1	905595	.	A	G	3887.15	PASS	AF=0.01


$ java -jar dist/jvarkit.jar biostar9556602 --sorter HIGHEST_AF  jeter.vcf --filter ZZZZZZZZZZZZ | grep -v "##"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905593	.	A	T	3887.15	PASS	AF=0.1
1	905594	rs867360694	A	TTT	3887.15	PASS	AF=0.3
1	905594	rs867360694	A	TT	3887.15	ZZZZZZZZZZZZ	AF=0.2
1	905594	rs867360694	A	T	3887.15	ZZZZZZZZZZZZ	AF=0.1
1	905595	.	A	T	3887.15	PASS	AF=0.1
1	905595	.	A	G	3887.15	ZZZZZZZZZZZZ	AF=0.01

$ java -jar dist/jvarkit.jar biostar9556602 --sorter LOWEST_AF  jeter.vcf  | grep -v "##"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905593	.	A	T	3887.15	PASS	AF=0.1
1	905594	rs867360694	A	T	3887.15	PASS	AF=0.1
1	905595	.	A	G	3887.15	PASS	AF=0.01


```

END_DOC

*/

@Program(name="biostar9556602",
description="Filtering of tricky overlapping sites in VCF",
keywords= {"vcf","filter"},
biostars=9556602,
jvarkit_amalgamion =  true,
menu="Biostars"
)
public class Biostar9556602 extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(Biostar9556602.class);
	private enum SORTER_TYPE {LOWEST_AF,HIGHEST_AF,LOWEST_DP,HIGHEST_DP,LOWEST_QUAL,HIGHEST_QUAL};
	@Parameter(names={"-filter","--filter"},description="if not blank, do not remove the variants but set the FILTER with this value")
	private String filterStr = null;

	@Parameter(names={"-sorter","--sorter"},description="Best variant sort type")
	private SORTER_TYPE sorterType = SORTER_TYPE.LOWEST_AF;

	private double getDouble(final VariantContext ctx,final String key) {
		if(!ctx.hasAttribute( key)) throw new IllegalArgumentException("INFO/"+ key+" missing for "+ctx);
		final List<Double> afs = ctx.getAttributeAsDoubleList(key, -1.0);
		if(afs.size()!=1)  throw new IllegalArgumentException("INFO/"+key+": expected one and only value for "+ctx);
		return afs.get(0).doubleValue();
		}
	
	private int getInt(final VariantContext ctx,final String key) {
		if(!ctx.hasAttribute( key)) throw new IllegalArgumentException("INFO/"+ key+" missing for "+ctx);
		final List<Integer> afs = ctx.getAttributeAsIntList(key, -1);
		if(afs.size()!=1)  throw new IllegalArgumentException("INFO/"+key+": expected one and only value for "+ctx);
		return afs.get(0).intValue();
		}
	
	private double getAF(final VariantContext ctx) {
		return getDouble(ctx,VCFConstants.ALLELE_FREQUENCY_KEY);
		}
	private int getDP(final VariantContext ctx) {
		return getInt(ctx,VCFConstants.DEPTH_KEY);
		}
	
	private double getQUAL(final VariantContext ctx) {
		if(!ctx.hasLog10PError()) return -1.0;
		return ctx.getPhredScaledQual();
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iter0, VariantContextWriter out) {
		try {
			final VCFHeader headerin = iter0.getHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(headerin);
			final ContigDictComparator dictCmp = new ContigDictComparator(dict);
			final VariantContextComparator ctxComparator = headerin.getVCFRecordComparator();
			final VCFHeader header = new VCFHeader(headerin);
			final VCFFilterHeaderLine filterLine;
			if(StringUtils.isBlank(this.filterStr)) {
				filterLine=null;
				}
			else
				{
				filterLine = new VCFFilterHeaderLine(this.filterStr, "Variants failing "+getProgramName());
				header.addMetaDataLine(filterLine);
				}
			
			final Comparator<VariantContext> comparator;
			switch(this.sorterType) {
				case LOWEST_AF: comparator = (V1,V2)->Double.compare(getAF(V1),getAF(V2)); break;
				case HIGHEST_AF: comparator = (V1,V2)->Double.compare(getAF(V2),getAF(V1)); break;
				case LOWEST_DP: comparator = (V1,V2)->Integer.compare(getDP(V1),getDP(V2)); break;
				case HIGHEST_DP: comparator = (V1,V2)->Integer.compare(getDP(V2),getDP(V1)); break;
				case LOWEST_QUAL: comparator = (V1,V2)->Double.compare(getQUAL(V1),getQUAL(V2)); break;
				case HIGHEST_QUAL: comparator = (V1,V2)->Double.compare(getQUAL(V2),getQUAL(V1)); break;
				default: throw new IllegalArgumentException(this.sorterType.name());
				}
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			out.writeHeader(header);
			try(EqualIterator<VariantContext> iter = new EqualIterator<>(iter0,(V1,V2)->{
				int i = dictCmp.compare(V1.getContig(),V2.getContig());
				if(i!=0) return i;
				i = Integer.compare(V1.getStart(), V2.getStart());
				if(i!=0) return i;
				return V1.getReference().compareTo(V2.getReference());
				}) ) {
				while(iter.hasNext()) {
					final List<VariantContext> array = iter.next().stream().
							sorted(comparator).
							collect(Collectors.toCollection(ArrayList::new));
					if(filterLine==null ) {
						out.add(array.get(0));
						}
					else
						{
						for(int i=1 /* skip first */;i< array.size();i++) {
							array.set(i,
								new VariantContextBuilder(array.get(i)).
								filter(filterLine.getID()).
								make());
							}
						//restore order
						Collections.sort(array,ctxComparator);
						for(VariantContext ctx:array) {
							out.add(ctx);
							}
						}
					
					}
				}
 			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
		
	public static void main(final String[] args) throws IOException
		{
		new Biostar9556602().instanceMainWithExit(args);
		}
		

	}
