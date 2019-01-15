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
* 2014 creation

*/

package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.List;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

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
@Program(name="sortvcfoninfo"
,description="Sort a VCF a field in the INFO column",
keywords={"vcf","sort","annotation"}
)
public class SortVcfOnInfo extends Launcher
	{
	private static final Logger LOG = Logger.build(SortVcfOnInfo.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-T","--tag","-t"},description="INFO tag",required=true)
    private String infoField=null;
    private VCFInfoHeaderLine infoDecl;
    private AbstractVCFCodec codec;
    
    @ParametersDelegate
    private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
    


    
    public SortVcfOnInfo()
	    {
	    }
    
    private class VcfLine
    	implements Comparable<VcfLine>
    	{
    	VariantContext ctx=null;
    	String line;
    	VcfLine()
    		{
    		}
    	public VcfLine(final String line)
    		{
    		this.line=line;
    		}
    	
    	public String getObject()
    		{
    		Object o=getContext().getAttribute(SortVcfOnInfo.this.infoField);
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
    			for(Object o2:L)
    				{
    				if(o2==null || o2.equals(".")) continue;
    				return o2.toString();
    				}
    			return null;
    			}
    		return o.toString();
    		}
    	@Override
    	public int compareTo(final VcfLine other)
    		{
    		final String o1=getObject();
    		final String o2=other.getObject();
    		if(o1==null)
    			{
    			if(o2==null) return line.compareTo(other.line);
    			return -1;
    			}
    		else if(o2==null)
    			{
    			return 1;
    			}
    		switch(infoDecl.getType())
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
    	public int hashCode() {
    		return line.hashCode();
    		}
    	@Override
    	public boolean equals(Object obj) {
    		return line.equals(VcfLine.class.cast(obj).line);
    		}
    	@Override
    	public String toString() {
    		return line;
    		}
    	VariantContext getContext()
    		{
    		if(this.ctx==null) this.ctx=SortVcfOnInfo.this.codec.decode(this.line);
    		return this.ctx;
    		}
    	}
    
    
	private class VariantCodec extends AbstractDataCodec<VcfLine>
		{
		@Override
		public VcfLine decode(final DataInputStream dis) throws IOException
			{
			VcfLine cpl=new VcfLine();
			try
				{
				cpl.line=readString(dis);
				}
			catch(final IOException err)
				{
				return null;
				}
			return cpl;
			}
		@Override
		public void encode(final DataOutputStream dos, final VcfLine s)
				throws IOException {
			writeString(dos,s.line);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	
	
	@Override
	public int doWork(final List<String> args)
		{
		CloseableIterator<VcfLine> iter=null;
		VariantContextWriter w=null;
		SortingCollection<VcfLine> sorted=null;
		LineIterator r=null;
		try {
			if(args.isEmpty())
			{
			LOG.info("reading from stdin");
			r=IOUtils.openStreamForLineIterator(stdin());
			}
		else if(args.size()==1)
			{
			String filename=args.get(0);
			LOG.info("Reading "+filename);
			r=IOUtils.openURIForLineIterator(filename);
			}
		else
			{
			LOG.error("Illegal number of arguments.");
			return -1;
			}

			
			
			final VCFUtils.CodecAndHeader ch=VCFUtils.parseHeader(r);
			VCFHeader header=ch.header;
			this.codec=ch.codec;
			this.infoDecl=header.getInfoHeaderLine(this.infoField);
			if(this.infoDecl==null)
				{
				final StringBuilder msg=new StringBuilder("VCF doesn't contain the INFO field :"+infoField+". Available:");
				for(VCFInfoHeaderLine vil:header.getInfoHeaderLines()) msg.append(" ").append(vil.getID());
				LOG.error(msg.toString());
				return -1;
				}
			
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
	
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header);
			w=super.openVariantContextWriter(this.outputFile);
			
			w.writeHeader(header);
			
			
			sorted=SortingCollection.newInstance(
					VcfLine.class,
					new VariantCodec(),
					(V1,V2)->V1.compareTo(V2),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorted.setDestructiveIteration(true);
			while(r.hasNext())
				{
				final VcfLine vc=new VcfLine(r.next());
				progress.watch(vc.getContext());
				sorted.add(vc);
				}
			CloserUtil.close(r);r=null;
			
			sorted.doneAdding();
			progress.finish();
			LOG.info("now writing...");
			iter =sorted.iterator();
			while(iter.hasNext())
				{
				w.add(iter.next().getContext());
				}
			iter.close();
			iter=null;
			w.close();
			w=null;
			return 0;
			} 
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(iter);
			try {
				if(sorted!=null) sorted.cleanup();
				} catch(Exception err){}
			CloserUtil.close(w);
			}
		}
	
    
 
	/**
	 * @param args
	 */
	public static void main(final String[] args) {
		new SortVcfOnInfo().instanceMainWithExit(args);
	}

}
