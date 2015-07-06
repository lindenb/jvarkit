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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Comparator;
import java.util.regex.Pattern;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

/**
 * Sort a VCF on the REFERENCE
 *
 */
public class SortVcfOnRef2 extends AbstractCommandLineProgram
	{
	private int max_records_in_ram=10000;
    private SAMSequenceDictionary dict=null;
    private Pattern tab=Pattern.compile("[\t]");
    @Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"SortVCFOnRef";
    	}
	    
    @Override
    public String getProgramDescription() {
    	return "Sort a VCF using the internal dictionary or an external reference order. ";
    	}
    
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
    	public ChromPosLine(String line)
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
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -R (fasta) indexed reference. Optional. The order of this reference will be used for sorting");
		out.println(" -T (dir) add tmp directory (optional)");
		out.println(" -N (int) max records in ram. default: "+ this.max_records_in_ram);
		super.printOptions(out);
		}
    
    @Override
    public int doWork(String[] args)
    	{    
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:T:N:"))!=-1)
			{
			switch(c)
				{
				case 'N': this.max_records_in_ram=Math.max(10,Integer.parseInt(opt.getOptArg()));break;
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;
				case 'R':
					{
					try
						{
						this.dict=new SAMSequenceDictionaryFactory().load(new File(opt.getOptArg()));
						}
					catch(IOException err)
						{
						error(err);
						return -1;
						}
					break;
					}
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
		
    	
		
		try
			{
			int ret=0;
			if(opt.getOptInd()==args.length)
				{
				info("reading from stdin");
				ret=sortvcf(new BufferedReader(new InputStreamReader(System.in)));
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading "+filename);
				BufferedReader in=IOUtils.openURIForBufferedReading(filename);
				ret= sortvcf(in);
				CloserUtil.close(in);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			return ret;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
    	}
    
   
	protected int sortvcf(BufferedReader in) throws IOException 
    	{
    	VCFUtils.CodecAndHeader cah =VCFUtils.parseHeader(in);
    	
    	VCFHeader h2=new  VCFHeader(cah.header);
    	if(this.dict!=null)
    		{
    		h2.setSequenceDictionary(this.dict);
    		}
    	else
    		{
    		this.dict= h2.getSequenceDictionary();
    		if(this.dict==null)
    			{
    			throw new IOException("No internal sequence dictionay found in input");
    			}
    		}
    	
    	h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		
		if(this.dict.isEmpty())
			{
			warning("SEQUENCE DICTIONARY IS EMPTY/NULL");
			}
		
    	CloseableIterator<ChromPosLine> iter=null;
    	SortingCollection<ChromPosLine> array=null;
    	VariantContextWriter w =null;
    	try {
			array= SortingCollection.newInstance(
					ChromPosLine.class,
					new VariantCodec(),
					new VariantComparator(),
					this.max_records_in_ram,
					this.getTmpDirectories()
					);
			array.setDestructiveIteration(true);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.dict);
			String line;
			while((line=in.readLine())!=null)
				{
				ChromPosLine cpl=new ChromPosLine(line);
				progress.watch(cpl.tid,cpl.pos);
				array.add(cpl);
				}
			array.doneAdding();
			progress.finish();
			
			 w = VCFUtils.createVariantContextWriterToStdout();
			w.writeHeader(h2);
			
			iter=array.iterator();
			while(iter.hasNext())
				{
				w.add(cah.codec.decode(iter.next().line));
				if(System.out.checkError()) break;
				}
			return 0;
			}
    	catch (Exception e)
    		{
			error(e);
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
