/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.util.Comparator;
import java.util.List;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

/**
 * Sort a VCF on the INFO field
 *
 */
public class SortVcfOnInfo extends AbstractKnimeApplication
	{
    private String infoField=null;
    private  VCFInfoHeaderLine infoDecl;
    private AbstractVCFCodec codec;
    private int maxRecordsInRam=5000;
    private int countFilteredVariants=0;
    
    public SortVcfOnInfo()
    {	
    	
    }
    
    @Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX +"SortVCFOnInfo";
    	}
	    
    @Override
    public String getProgramDescription() {
    	return "Sort a VCF a field in the INFO column";
    	}
    
    public void setMaxRecordsInRam(int maxRecordsInRam) {
		this.maxRecordsInRam = maxRecordsInRam;
		}
    public void setInfoField(String infoField) {
		this.infoField = infoField;
		}
    
    private class VcfLine
    	implements Comparable<VcfLine>
    	{
    	VariantContext ctx=null;
    	String line;
    	VcfLine()
    		{
    		}
    	public VcfLine(String line)
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
    			return array==null || array.length==0 ?null:array[0].toString();
    			}
    		if(o instanceof List)
    			{
    			@SuppressWarnings("rawtypes")
				List L=(List)o;
    			if(L.isEmpty()) return null;
    			return L.get(0).toString();
    			}
    		return o.toString();
    		}
    	@Override
    	public int compareTo(VcfLine other)
    		{
    		String o1=getObject();
    		String o2=other.getObject();
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
		public VcfLine decode(DataInputStream dis) throws IOException
			{
			VcfLine cpl=new VcfLine();
			try
				{
				cpl.line=readString(dis);
				}
			catch(IOException err)
				{
				return null;
				}
			return cpl;
			}
		@Override
		public void encode(DataOutputStream dos, VcfLine s)
				throws IOException {
			writeString(dos,s.line);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	
	private class VariantComparator implements Comparator<VcfLine>
		{
		@Override
		public int compare(VcfLine o1, VcfLine o2)
			{
			return o1.compareTo(o2);
			}
		}
	
	private void doWork(LineIterator r) throws IOException
		{
		VCFUtils.CodecAndHeader ch=VCFUtils.parseHeader(r);
		VCFHeader header=ch.header;
		this.codec=ch.codec;
		this.infoDecl=header.getInfoHeaderLine(this.infoField);
		if(this.infoDecl==null)
			{
			StringBuilder msg=new StringBuilder("VCF doesn't contain the INFO field :"+infoField+". Available:");
			for(VCFInfoHeaderLine vil:header.getInfoHeaderLines()) msg.append(" ").append(vil.getID());
			throw new IOException(msg.toString());
			}
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);

		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		VariantContextWriter w=null;
		if(getOutputFile()==null)
			{
			w=VCFUtils.createVariantContextWriterToStdout();
			}
		else
			{
			w=VCFUtils.createVariantContextWriter(getOutputFile());
			}
		w.writeHeader(header);
		
		
		SortingCollection<VcfLine> sorted=SortingCollection.newInstance(
				VcfLine.class,
				new VariantCodec(),
				new VariantComparator(),
				maxRecordsInRam,
				getTmpDirectories().get(0)
				);
		sorted.setDestructiveIteration(true);
		while(r.hasNext())
			{
			VcfLine vc=new VcfLine(r.next());
			progress.watch(vc.getContext());
			sorted.add(vc);
			}
		sorted.doneAdding();
		progress.finish();
		
		countFilteredVariants=0;
		CloseableIterator<VcfLine> iter=sorted.iterator();
		while(iter.hasNext())
			{
			w.add(iter.next().getContext());
			countFilteredVariants++;
			}
		iter.close();
		w.close();
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -F (field) INFO field for sorting. REQUIRED.");
		out.println(" -T (dir) "+getMessageBundle("add.tmp.dir")+" (optional)");
		out.println(" -N (int) "+getMessageBundle("max.records.in.ram")+" default: "+maxRecordsInRam);
		super.printOptions(out);
		}
    
	@Override
	public int initializeKnime() {
    	if(infoField==null || infoField.trim().isEmpty())
			{
			error("undefined or empty INFO field.");
			return -1;
			}

		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		
		super.disposeKnime();
		}
	
	
	public int getVariantCount()
		{
		return this.countFilteredVariants;
		}

	
    @Override
    public int doWork(String[] args)
    	{
    
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"F:T:N:o:"))!=-1)
			{
			switch(c)
				{
				case 'o': setOutputFile(opt.getOptArg());break;
				case 'F': infoField=opt.getOptArg();break;
				case 'N': setMaxRecordsInRam(Integer.parseInt(opt.getOptArg()));break;
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
    	}
		
    @Override
    public int executeKnime(List<String> args)
    	{
    	LineIterator r=null;
		try
			{
			int ret=0;
			if(args.isEmpty())
				{
				info("reading from stdin");
				r=IOUtils.openStdinForLineIterator();
				}
			else if(args.size()==1)
				{
				String filename=args.get(0);
				info("Reading "+filename);
				r=IOUtils.openURIForLineIterator(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			doWork(r);
			return ret;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			}
    	}
 
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SortVcfOnInfo().instanceMainWithExit(args);
	}

}
