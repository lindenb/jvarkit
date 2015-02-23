package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.vcf.VCFConstants;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

/**
 * Sort a VCF on the REFERENCE
 *
 */
@Deprecated /** use picard */
public class SortVcfOnRef2 extends AbstractCommandLineProgram
	{
    private SAMSequenceDictionary dict=null;
    private SortingCollectionFactory<ChromPosLine> sortingCollectionFactory=new SortingCollectionFactory<ChromPosLine>();
    
    @Override
    protected String getOnlineDocUrl() {
    	return "https://github.com/lindenb/jvarkit/wiki/SortVCFOnRef";
    	}
	    
    @Override
    public String getProgramDescription() {
    	return "Sort a VCF using the internal dictionary or an external reference order. ";
    	}
    private class ChromPosLine
    	implements Comparable<ChromPosLine>
    	{
    	int tid;
    	int pos;
    	String line;
    	ChromPosLine()
    		{
    		}
    	public ChromPosLine(String line)
    		{
    		int tab1=line.indexOf('\t');
    		if(tab1==-1) throw new IllegalArgumentException("Bad VCF line in "+line);
    		String chrom=line.substring(0,tab1);
			this.tid=dict.getSequenceIndex(chrom);
			if(this.tid==-1) throw new RuntimeException("unknown chromosome "+ chrom+" in "+line);
			int tab2=line.indexOf('\t',tab1+1);
    		if(tab2==-1) throw new IllegalArgumentException("Bad VCF line in "+line);
    		try
    			{
    			this.pos=Integer.parseInt(line.substring(tab1+1, tab2));
    			}
    		catch(NumberFormatException err)
    			{
    			throw new IllegalArgumentException("Bad VCF line in "+line);
    			}
    		this.line=line;
    		}
    	@Override
    	public int compareTo(ChromPosLine o)
    		{
    		int i=this.tid-o.tid;
    		if(i!=0) return i;
    		i=this.pos-o.pos;
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
			cpl.line= readString(dis);
			return cpl;
			}
		@Override
		public void encode(DataOutputStream dos, ChromPosLine s)
				throws IOException {
			dos.writeInt(s.tid);
			dos.writeInt(s.pos);
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
	
    private static int containsvcfformat(String s)
    	{
    	if( s.startsWith("##fileformat=")) return 0;
    	if (s.startsWith("##format=")) return 1;
    	return 2;
    	}
    
    private int contig(String line)
		{
		String tokens1[]=line.split("[<=,]");
		for(int i=0;i+1< tokens1.length;++i)
			{
			if(tokens1[i].equals("ID"))
				{
				return dict.getSequenceIndex(tokens1[i+1]);
				}
			}
		return -1;
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -R (fasta) indexed reference. Optional. The order of this reference will be used for sorting");
		out.println(" -T (dir) add tmp directory (optional)");
		out.println(" -N (int) max records in ram. default: "+sortingCollectionFactory.getMaxRecordsInRAM());
		super.printOptions(out);
		}
    
    @Override
    public int doWork(String[] args)
    	{
    	sortingCollectionFactory.setMaxRecordsInRAM(50000);
    	sortingCollectionFactory.setCodec(new VariantCodec());
    	sortingCollectionFactory.setComparator(new VariantComparator());
    	sortingCollectionFactory.setComponentType(ChromPosLine.class);
    
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:T:N:"))!=-1)
			{
			switch(c)
				{
				case 'N': sortingCollectionFactory.setMaxRecordsInRAM(Integer.parseInt(opt.getOptArg()));break;
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
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		sortingCollectionFactory.setTmpDirs(this.getTmpDirectories());
    	
		
		try
			{
			int ret=0;
			if(opt.getOptInd()==args.length)
				{
				info("reading from stdin");
				ret=doWork(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading "+filename);
				InputStream in=IOUtils.openURIForReading(filename);
				ret=doWork(in);
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
    
   
    @SuppressWarnings("resource")
	protected int doWork(InputStream in) throws IOException 
    	{
		final String contigStart=VCFConstants.CONTIG_HEADER_START+"=";
		String chromHeader=null;
    	LineReader lr=LineReaderUtil.fromBufferedStream(in);
    	LineIterator liter=new LineIteratorImpl(lr);
    	
    	info("reading header");
    	List<String> headerLines=new ArrayList<String>();
		while(liter.hasNext())
			{
			String line=liter.peek();
			if(line.startsWith("#CHROM"))
				{
				chromHeader=liter.next();
				continue;
				}
			if(!line.startsWith("##"))
				{
				break;
				}
			line=liter.next();
			if(containsvcfformat(line)<2)
				{
				if(!headerLines.isEmpty() &&
					containsvcfformat(headerLines.get(0))<containsvcfformat(line))
					{
					//continue
					}
				else if(!headerLines.isEmpty() &&
						containsvcfformat(headerLines.get(0))< 2 &&
						containsvcfformat(headerLines.get(0))>=containsvcfformat(line))
					{
					headerLines.set(0, line);//replace
					}
				else
					{
					headerLines.add(0, line);//insert front
					}
				}
			else if(this.dict!=null && line.startsWith(contigStart))
				{
				//using external dict, skip those lines
				}
			else
				{
				headerLines.add(line);
				}
			}
		if(headerLines.isEmpty())
			{
			error("no header found in input VCF");
			return -1;
			}
		if(chromHeader==null)
			{
			error("The #CHROM line was not found in header");
			return -1;
			}
		
		//insert the external dict
		if(dict!=null)
			{
			for(SAMSequenceRecord ssr:this.dict.getSequences())
				{
				headerLines.add(VCFUtils.samSequenceRecordToVcfContigLine(ssr));
				}
			}
		
		//#CHROM line below
		Collections.sort(headerLines,new Comparator<String>()
				{
				//sort, we want the format as the first line.
				@Override
				public int compare(String o1, String o2)
					{ 
					
					int i=containsvcfformat(o1)-containsvcfformat(o2);
					if(i!=0) return i;
					if(dict!=null)//use external dict
						{
						if( o1.startsWith(contigStart) &&
							o2.startsWith(contigStart)
							)
							{
							i=contig(o1)-contig(o2);
							if(i!=0) return i;
							}
						}
					else
						{
						/*  This sort is guaranteed to be stable: equal elements will not be reordered as a result of the sort. */
						if( o1.startsWith(contigStart) &&
							o2.startsWith(contigStart)
							)
							{
							return 0;
							}
						}
					return o1.compareTo(o2);
					}
				});
		
		
		if(dict==null)//use internal dict
			{
			info("Using internal dictionary");
			List<SAMSequenceRecord> ssl=new ArrayList<SAMSequenceRecord>();
			for(String contig:headerLines)
				{
				if(!contig.startsWith(contigStart)) continue;
				ssl.add(VCFUtils.contigLineToSamSequenceRecord(contig));
				}
			if(ssl.isEmpty())
				{
				warning("There is no internal DICTIONARY '##contig=' in the VCF.");
				}
			this.dict=new SAMSequenceDictionary(ssl);
			}
		
		if(this.dict.isEmpty())
			{
			warning("SEQUENCE DICTIONARY IS EMPTY/NULL");
			}
		
    	CloseableIterator<ChromPosLine> iter=null;
    	SortingCollection<ChromPosLine> array=null;
    	try {
    		
    		for(String line:headerLines)
    			{
    			System.out.println(line);
    			}
    		//write the #CHROM header
    		System.out.println(chromHeader);
    			
    		
			
			array=this.sortingCollectionFactory.make();
			array.setDestructiveIteration(true);
			info("Reading body");
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.dict);
			while(liter.hasNext())
				{
				String line=liter.next();
				if(line.startsWith("#"))
					{
					throw new IOException("bad VCF line in "+line);
					}
				ChromPosLine cpl=new ChromPosLine(line);
				progress.watch(cpl.tid,cpl.pos);
				array.add(cpl);
				}
			array.doneAdding();
			lr.close();
			progress.finish();
			
			iter=array.iterator();
			while(iter.hasNext())
				{
				System.out.println(iter.next().line);
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
    		CloserUtil.close(liter);
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
