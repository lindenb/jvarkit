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
package com.github.lindenb.jvarkit.tools.misc;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

## Example

```
$ java -jar dist/vcfdistancevariants.jar src/test/resources/toy.vcf.gz
##fileformat=VCFv4.2
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=DIST_NEXT,Number=1,Type=Integer,Description="Distance to next variant">
##INFO=<ID=DIST_PREV,Number=1,Type=Integer,Description="Distance to previous variant">
(..)
##contig=<ID=ref,length=45>
##contig=<ID=ref2,length=40>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
ref	11	.	A	C	0.52	.	AC1=1;AF1=0.498013;BQB=1;DIST_NEXT=3;DP=3;DP4=2,0,1,0;FQ=-6.97257;MQ=30;MQ0F=0;MQB=1;PV4=1,1,1,0.106148;RPB=1;SGB=-0.379885	GT:PL	0/1:21,0,46
ref	14	.	AG	AGGG	0.71	.	AC1=2;AF1=1;DIST_PREV=3;DP=3;DP4=0,0,1,0;FQ=-37.5301;IDV=2;IMF=0.666667;INDEL;MQ=30;MQ0F=0;SGB=-0.379885	GT:PL	0/1:30,3,0
ref2	14	.	C	T	0.09	.	AC1=1;AF1=0.487148;BQB=1;DIST_NEXT=0;DP=6;DP4=3,0,1,0;FQ=-14.1557;MQ=30;MQ0F=0;MQB=1;PV4=1,0,1,0.0285955;RPB=1;SGB=-0.379885	GT:PL	0/0:13,0,64
ref2	14	.	CAA	CAAATAA	32.43	.	AC1=2;AF1=1;DIST_NEXT=1;DIST_PREV=0;DP=6;DP4=0,0,3,0;FQ=-43.5253;IDV=1;IMF=0.166667;INDEL;MQ=30;MQ0F=0;SGB=-0.511536;VDB=0.354794	GT:PL	1/1:72,9,0
ref2	17	.	T	A	0.14	.	AC1=1;AF1=0.491968;BQB=1;DIST_PREV=1;DP=6;DP4=4,0,1,0;FQ=-12.2521;MQ=30;MQ0F=0;MQB=1;PV4=1,1,1,0.201057;RPB=1;SGB=-0.379885	GT:PL	0/0:15,0,72
```

END_DOC

*/
@Program(
		name="vcfdistancevariants",
		description="Annotate variants with the distance between previous and next variant.",
		keywords={"vcf","annotation","distance"},
		creationDate="20190410",
		modificationDate="20230510"	,
		jvarkit_amalgamion = true,
		menu="VCF Manipulation"

		)
public class VcfDistanceBetweenVariants extends OnePassVcfLauncher{
	private static final Logger LOG = Logger.of(VcfDistanceBetweenVariants.class);

	@Parameter(names={"-p","--prefix"},description="INFO Attribute Prefix")
	private String prefix="DIST_";
    
	
	private static class Variant {
		VariantContext ctx;
		Integer dist_to_prev=null;
		Integer dist_to_next=null;
		}
	
    protected static int distance(final VariantContext v1,final VariantContext v2) {
    	if(v1.overlaps(v2)) 
    		{
    		return 0;
    		}
    	else if(v1.getEnd() < v2.getStart())
    		{
    		return v2.getStart() - v1.getEnd();
    		}
    	else if( v2.getEnd() < v1.getStart()) {
    		return v1.getStart() - v2.getEnd();
    	} else
    	{
    		throw new IllegalStateException(""+v1+" "+v2);
    	}
    	}
    
    @Override
    protected Logger getLogger() {
    	return LOG;
    	}
    
    
    @Override
    protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter w) {
    	final VCFHeader header= new VCFHeader(in.getHeader());        	
    	final VCFInfoHeaderLine infoPrev= new VCFInfoHeaderLine(prefix+"PREV", 1,VCFHeaderLineType.Integer ,"Distance to previous variant.");
    	final VCFInfoHeaderLine infoNext= new VCFInfoHeaderLine(prefix+"NEXT", 1,VCFHeaderLineType.Integer ,"Distance to next variant.");
    	
    	header.addMetaDataLine(infoPrev);
    	header.addMetaDataLine(infoNext);
    	
    	JVarkitVersion.getInstance().addMetaData(this, header);
    	Variant prev = null;
    	w.writeHeader(header);
    	for(;;) {
    		if(prev!=null) {
				final VariantContextBuilder vcb=new VariantContextBuilder(prev.ctx);
        		vcb.rmAttribute(infoPrev.getID());
        		vcb.rmAttribute(infoNext.getID());
    			if(prev.dist_to_prev!=null) vcb.attribute(infoPrev.getID(),prev.dist_to_prev);
    			if(prev.dist_to_next!=null) vcb.attribute(infoNext.getID(),prev.dist_to_next);
    			w.add(vcb.make());
    			}
        	
    		final VariantContext ctx =  in.hasNext()?in.next():null;
    	
    		if(ctx==null) break;
    		
    		final Variant variant = new Variant();
    		variant.ctx  = ctx;
    		final VariantContext next = in.hasNext()?in.peek():null;
    		if(next!=null && ctx.contigsMatch(next) ) {
    			variant.dist_to_next =  distance(ctx,next);
    			}

    		if(prev!=null && ctx.contigsMatch(prev.ctx) ) {
    			variant.dist_to_prev =  distance(prev.ctx,ctx);
    			}
    		prev=variant;        		
    		}
    	return 0;
		}
    	
    	
    public static void main(final String[] args) {
		new VcfDistanceBetweenVariants().instanceMainWithExit(args);
	}

}
