package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Comparator;

import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.IOUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

public class SortVcfOnRef extends CommandLineProgram
	{
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+ "Sort a VCF using the reference order. ";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF file (or stdin).",optional=true)
	public File IN=null;
    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="output file (or stdout).",optional=true)
	public File OUT=null;

    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference file.",optional=false)
	public File REF=null;
   
    private SAMSequenceDictionary dict=null;
    
	@Override
	public String getProgramVersion() {
		return "1.0";
	}
	    
	private static class VariantCodec extends AbstractDataCodec<String>
		{
		@Override
		public String decode(DataInputStream dis) throws IOException
			{
			try
				{
				return dis.readUTF();
				}
			catch(IOException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, String s)
				throws IOException {
			dos.writeUTF(s);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	
	private class VariantComparator implements Comparator<String>
		{
		private int ref(String chrom,String line)
			{
			int refId=dict.getSequenceIndex(chrom);
			if(refId==-1) throw new RuntimeException("unknown chromosome "+ chrom+" in "+line);
			return refId;
			}
		@Override
		public int compare(String o1, String o2)
			{
			String tokens1[]=o1.split("\t",5);
			String tokens2[]=o2.split("\t",5);
			
			int i=ref(tokens1[0],o1) - ref(tokens2[0],o2)  ;
			if(i!=0) return i;
			i= Integer.parseInt(tokens1[1])-Integer.parseInt(tokens2[1]);
			if(i!=0) return i;
			i=tokens1[3].compareTo(tokens2[3]);
			return i;
			}
		}
	
    
    @Override
    protected int doWork() 
    	{
    	
    	CloseableIterator<String> iter=null;
		BufferedReader in=null;
		PrintWriter out=null;
    	try {
        	IndexedFastaSequenceFile ref=new IndexedFastaSequenceFile(REF);
        	this.dict=ref.getSequenceDictionary();

			if(IN==null)
				{
				in=new BufferedReader(new InputStreamReader(System.in));
				}
			else
				{
				in=IOUtils.openFileForBufferedReading(IN);
				}
			if(OUT==null)
				{
				out=new PrintWriter(out);
				}
			else
				{
				out=new PrintWriter(IOUtils.openFileForBufferedWriting(OUT));
				}
			
			SortingCollection<String> array=SortingCollection.newInstance(
				String.class,
				new VariantCodec(),
				new VariantComparator(),
				super.MAX_RECORDS_IN_RAM
				);
			array.setDestructiveIteration(true);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.startsWith("#"))
					{
					out.println(line);
					continue;
					}
				array.add(line);
				}
			array.doneAdding();
			iter=array.iterator();
			while(iter.hasNext())
				{
				out.println(iter.next());
				}
			
		} catch (Exception e) {
			e.printStackTrace();
			return -1;
			}
    	finally
	    	{
	    	if(in!=null) try { in.close();} catch(Exception err){}
	    	if(out!=null) try { out.flush();} catch(Exception err){}
	    	if(out!=null) try { out.close();} catch(Exception err){}
	    	if(iter!=null) try { iter.close();} catch(Exception err){}
	    	}
    	return 0;
    	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	new SortVcfOnRef().instanceMainWithExit(args);

	}

}
