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

package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFRecordCodec;

/**
BEGIN_DOC

##Example

```bash
]$ curl  "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
   java -jar dist/sortvcfoninfo.jar -F BaseQRankSum | grep -vE "^#" 

chr10	1142208	.	T	C	3404.30	.	AC=8;AF=1.00;AN=8;
chr10	135336656	.	G	A	38.34	.	AC=4;AF=1.00;AN=4;
chr10	52004315	.	T	C	40.11	.	AC=4;AF=1.00;AN=4;
chr10	52497529	.	G	C	33.61	.	AC=4;AF=1.00;AN=4;
chr10	126678092	.	G	A	89.08	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-3.120;
chr16	72057435	.	C	T	572.98	.	AC=1;AF=0.13;AN=8;BaseQRankSum=-2.270;
chr10	48003992	.	C	T	1047.87	.	AC=4;AF=0.50;AN=8;BaseQRankSum=-0.053;
chr10	135210791	.	T	C	65.41	.	AC=4;AF=0.50;AN=8;BaseQRankSum=2.054;
chr10	135369532	.	T	C	122.62	.	AC=2;AF=0.25;AN=8;BaseQRankSum=2.118;
```

END_DOC
 */
@Program(name="sortvcfoninfo",
description="Sort a VCF a field in the INFO column",
keywords={"vcf","sort","annotation"},
creationDate="20140218",
modificationDate="20201204"
)
public class SortVcfOnInfo extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(SortVcfOnInfo.class).make();


	@Parameter(names={"-T","--tag","-t"},description="INFO tag. Special words are '<ID>' to sort on ID, and <QUAL> to sort on QUAL ",required=true)
    private String infoField=null;
	@Parameter(names={"-r","--reverse"},description="reverse order")
    private boolean reverse_it = false;

    private VCFInfoHeaderLine infoDecl;
    
    @ParametersDelegate
    private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();

    
    
    public SortVcfOnInfo()
	    {
	    }
    
    @Override
    protected Logger getLogger() {
    	return LOG;
    	}
    
    
    public String getObject(final VariantContext ctx)
		{
		Object o= ctx.getAttribute(SortVcfOnInfo.this.infoField);
		if(o==null) return null;
		if(o.getClass().isArray())
			{
			Object array[]=(Object[])o;
			o = Arrays.asList(array);
			}
		if(o instanceof List)
			{
			@SuppressWarnings("rawtypes")
			final List L=(List)o;
			if(L.isEmpty()) return null;
			for(final Object o2:L)
				{
				if(o2==null || o2.equals(".")) continue;
				return o2.toString();
				}
			return null;
			}
		return o.toString();
		}
    
    private int defaultCompare(final VariantContext vc1,final VariantContext vc2) {
    	int i = vc1.getContig().compareTo(vc2.getContig());
		if(i!=0) return i;
		i = Integer.compare(vc1.getStart(), vc2.getStart());
		if(i!=0) return i;
		i = Integer.compare(vc1.getEnd(), vc2.getEnd());
		if(i!=0) return i;
		i = vc1.getReference().compareTo(vc2.getReference());
		return i;
    }
    
    private int compareID(final VariantContext vc1,final VariantContext vc2) 
    	{
    	final String id1 = vc1.hasID()?null:vc1.getID();
    	final String id2 = vc2.hasID()?null:vc2.getID();
    	if(id1!=null)
    		{
    		if(id2!=null) {
    			int i= id1.compareTo(id2);
    			if(i!=0) return i;
    			return defaultCompare(vc1,vc2);
    			}
    		else
    			{
    			return 1;
    			}
    		}
    	else {
    		 if(id2!=null) {
    			 return -1;
    		 	}
    		 else
	    		 {
	    		 return defaultCompare(vc1,vc2);
	    		 }
    		}
    	
    	}
    private int compareQUAL(final VariantContext vc1,final VariantContext vc2) 
		{
    	final Double id1 = vc1.hasLog10PError()?null:vc1.getPhredScaledQual();
    	final Double id2 = vc2.hasLog10PError()?null:vc2.getPhredScaledQual();
    	if(id1!=null)
    		{
    		if(id2!=null) {
    			int i= id1.compareTo(id2);
    			if(i!=0) return i;
    			return defaultCompare(vc1,vc2);
    			}
    		else
    			{
    			return 1;
    			}
    		}
    	else {
    		 if(id2!=null) {
    			 return -1;
    		 	}
    		 else
	    		 {
	    		 return defaultCompare(vc1,vc2);
	    		 }
    		}
    	
    	}

    
    	
    private int compareVariants(final VariantContext vc1,final VariantContext vc2) {
		final String o1=getObject(vc1);
		final String o2=getObject(vc2);
		if(o1==null)
			{
			if(o2==null) {
				return defaultCompare(vc1,vc2);
				}
			return -1;
			}
		else if(o2==null)
			{
			return 1;
			}
		switch(this.infoDecl.getType())
			{	
			case Float:
				{	
				return new BigDecimal(o1).compareTo(new BigDecimal(o2));
				}
			case Integer:
				{	
    			return new BigInteger(o1).compareTo(new BigInteger(o2));
				}
			default:
				{
				return o1.compareTo(o2);
				}
			}
    }
    
    @Override
    protected int doVcfToVcf(String inputName, VCFIterator r, VariantContextWriter w) {
    	SortingCollection<VariantContext> sorted=null;
		try {				
			final Comparator<VariantContext> cmp2;
			final VCFHeader header=r.getHeader();
			if(this.infoField!=null && this.infoField.equals("<ID>")) 
				{
				cmp2 = (A,B)->compareID(A, B);
				}
			else if(this.infoField!=null && this.infoField.equals("<QUAL>")) 
				{
				cmp2 = (A,B)->compareQUAL(A, B);
				}
			else
				{
				this.infoDecl=header.getInfoHeaderLine(this.infoField);
				if(this.infoDecl==null)
					{
					LOG.error("VCF doesn't contain the INFO field :"+infoField+". Available:" +
							header.getInfoHeaderLines().stream().map(H->H.getID()).collect(Collectors.joining(" ")));
					return -1;
					}
				cmp2 = (V1,V2)->compareVariants(V1,V2);
				}	
				
			final Comparator<VariantContext> cmp = (reverse_it?cmp2.reversed():cmp2);
			
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header);
			
			sorted=SortingCollection.newInstance(
					VariantContext.class,
					new VCFRecordCodec(header),
					cmp,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorted.setDestructiveIteration(true);
			while(r.hasNext())
				{
				sorted.add(r.next());
				}
			CloserUtil.close(r);r=null;
			
			sorted.doneAdding();			
			w.writeHeader(header);
			
			try(CloseableIterator<VariantContext> iter =sorted.iterator()) {
				while(iter.hasNext())
					{
					w.add(iter.next());
					}
				}
			return 0;
			} 
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			try {
				if(sorted!=null) sorted.cleanup();
				} catch(Exception err){}
			CloserUtil.close(w);
			}
		}
	

	public static void main(final String[] args) {
		new SortVcfOnInfo().instanceMainWithExit(args);
	}

}
