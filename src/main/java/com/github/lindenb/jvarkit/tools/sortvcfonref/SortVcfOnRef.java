package com.github.lindenb.jvarkit.tools.sortvcfonref;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

/**
 * Sort a VCF on the REFERENCE
 *
 */
public class SortVcfOnRef extends AbstractCommandLineProgram
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
	
    private static int containsvcfformat(String s)
    	{
    	if( s.startsWith("##fileformat=")) return 0;
    	if (s.startsWith("##format=")) return 1;
    	return 2;
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
				out=new PrintWriter(System.out);
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
			List<String> heads=new ArrayList<String>();
			while((line=in.readLine())!=null)
				{
				if(line.startsWith("#"))
					{
					if(line.startsWith("##"))
						{
						heads.add(line);
						continue;
						}
					//#CHROM line below
					Collections.sort(heads,new Comparator<String>()
							{
							//sort, we want the format as the first line.
							@Override
							public int compare(String o1, String o2)
								{ 
								int i=containsvcfformat(o1)-containsvcfformat(o2);
								if(i!=0) return i;
								return o1.compareTo(o2);
								}
							});
					for(String s:heads)
						{
						out.println(s);
						}
					heads.clear();
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
