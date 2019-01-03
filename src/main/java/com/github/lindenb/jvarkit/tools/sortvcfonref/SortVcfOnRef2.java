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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC

### Deprecated

Use picard SortVcf  http://broadinstitute.github.io/picard/command-line-overview.html#SortVcf.

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
	description="Sort a VCF using the internal dictionary or an external reference order (Deprecated: use picard SortVcf).",
	deprecatedMsg="use picard sortvcf",
	keywords={"vcf","sort"}
	)
public class SortVcfOnRef2 extends Launcher
	{
	
	private static final Logger LOG = Logger.build(SortVcfOnRef2.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File refdict = null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	

	private SAMSequenceDictionary dict=null;
    private final Pattern tab=Pattern.compile("[\t]");
	    
    
    
    private class ChromPosLine
    	implements Comparable<ChromPosLine>
    	{
    	int tid;
    	Integer pos;
    	Allele ref;
    	String line;
    	ChromPosLine()
    		{
    		}
    	public ChromPosLine(final String line)
    		{
    		String tokens[]= tab.split(line, 5);
    		if(tokens.length<5) throw new IllegalArgumentException("Bad VCF line in "+line); 
    		String chrom=tokens[0];
			this.tid=dict.getSequenceIndex(chrom);
			if(this.tid==-1) throw new RuntimeException("unknown chromosome "+ chrom+" in "+line);
			
    		try
    			{
    			this.pos=new Integer(tokens[1]);
    			}
    		catch(NumberFormatException err)
    			{
    			throw new IllegalArgumentException("Bad POS in VCF line in "+line);
    			}
    		this.ref = Allele.create(tokens[3],true);
    		this.line=line;
    		}
    	@Override
    	public int compareTo(ChromPosLine o)
    		{
    		int i=this.tid-o.tid;
    		if(i!=0) return i;
    		i=this.pos.compareTo(o.pos);
    		if(i!=0) return i;
    		i=this.ref.compareTo(o.ref);
    		if(i!=0) return i;
    		return this.line.compareTo(o.line);
    		}
    	@Override
    	public int hashCode() {
    		return line.hashCode();
    		}
    	@Override
    	public boolean equals(Object obj) {
    		return line.equals(ChromPosLine.class.cast(obj).line);
    		}
    	@Override
    	public String toString() {
    		return line;
    		}
    	}
	private class VariantCodec extends AbstractDataCodec<ChromPosLine>
		{
		@Override
		public ChromPosLine decode(DataInputStream dis) throws IOException
			{
			ChromPosLine cpl=new ChromPosLine();
			try
				{
				cpl.tid=dis.readInt();
				}
			catch(IOException err)
				{
				return null;
				}
			cpl.pos=dis.readInt();
			cpl.ref = Allele.create(readString(dis), true);
			cpl.line= readString(dis);
			return cpl;
			}
		@Override
		public void encode(DataOutputStream dos, ChromPosLine s)
				throws IOException {
			dos.writeInt(s.tid);
			dos.writeInt(s.pos);
			writeString(dos, s.ref.getBaseString());
			writeString(dos,s.line);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	
	private class VariantComparator implements Comparator<ChromPosLine>
		{
		@Override
		public int compare(ChromPosLine o1, ChromPosLine o2)
			{
			return o1.compareTo(o2);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		BufferedReader in =null;
		try
			{
			final String inputName=oneFileOrNull(args);
			if(inputName==null)
				{
				LOG.info("reading from stdin");
				in = IOUtils.openStreamForBufferedReader(stdin());
				return sortvcf(in);
				}
			else 
				{
				LOG.info("Reading "+inputName);
				in=IOUtils.openURIForBufferedReading(inputName);
				}
	
			return sortvcf(in);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
    	}
    
   
	protected int sortvcf(BufferedReader in) throws IOException 
    	{
		if(this.refdict!=null) {
			LOG.info("load dict from "+this.refdict);
			this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.refdict);
			if(this.dict==null) {
				LOG.error("cannot find sam sequence dictionary from "+refdict);
			}
		}
		
		final VCFUtils.CodecAndHeader cah =VCFUtils.parseHeader(in);
    	final VCFHeader h2=new  VCFHeader(cah.header);
    	if(this.dict!=null)
    		{
    		h2.setSequenceDictionary(this.dict);
    		}
    	else
    		{
    		this.dict= h2.getSequenceDictionary();
    		if(this.dict==null)
    			{
    			LOG.error("No internal sequence dictionay found in input");
    			return -1;
    			}
    		}
    	
    	addMetaData(h2);
		
		if(this.dict.isEmpty())
			{
			LOG.warn("SEQUENCE DICTIONARY IS EMPTY/NULL");
			}
		
    	CloseableIterator<ChromPosLine> iter=null;
    	SortingCollection<ChromPosLine> array=null;
    	VariantContextWriter w =null;
    	try {
			array= SortingCollection.newInstance(
					ChromPosLine.class,
					new VariantCodec(),
					new VariantComparator(),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			array.setDestructiveIteration(true);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.dict);
			String line;
			while((line=in.readLine())!=null)
				{
				final ChromPosLine cpl=new ChromPosLine(line);
				progress.watch(cpl.tid,cpl.pos);
				array.add(cpl);
				}
			array.doneAdding();
			progress.finish();
			
			w = super.openVariantContextWriter(outputFile);
			w.writeHeader(h2);
			
			iter=array.iterator();
			while(iter.hasNext())
				{
				w.add(cah.codec.decode(iter.next().line));
				if(w.checkError()) break;
				}
			return RETURN_OK;
			}
    	catch (Exception e)
    		{
			LOG.error(e);
			return -1;
			}
    	finally
	    	{
    		CloserUtil.close(w);
	    	CloserUtil.close(iter);
	    	if(array!=null) array.cleanup();
	    	}
    	
    	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SortVcfOnRef2().instanceMainWithExit(args);
	}

}
