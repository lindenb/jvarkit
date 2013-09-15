package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.List;

import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFHeader;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Log;
import net.sf.picard.vcf.VcfIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class VCFMerge extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(VCFMerge.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Merge VCFs";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
    		doc="VCF files to process.",
    		optional=false)
	public List<File> IN=new ArrayList<File>();

    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME,
    		doc="Reference file.",
    		optional=false
    		)
	public File REF=null;

    
	private SAMSequenceDictionary dictionary;
	
	private static class VCFHandler
		{
		File vcfFile;
		StringBuilder headerStr=new StringBuilder();
		}
	
	private static class VariantOfFile
		{
		int fileIndex;
		String line;
		boolean same(VariantOfFile var)
			{
			return false;
			}
		}
	
	private static class VariantCodec
		extends AbstractDataCodec<VariantOfFile>
		{
		@Override
		public VariantOfFile decode(DataInputStream dis) throws IOException
			{
			try
				{
				VariantOfFile o=new VariantOfFile();
				o.fileIndex=dis.readInt();
				o.line=dis.readUTF();
				return o;
				}
			catch(IOException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, VariantOfFile s)
				throws IOException {
			dos.writeInt(s.fileIndex);
			dos.writeUTF(s.line);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	
	private class VariantComparator implements Comparator<VariantOfFile>
		{
		private int ref(String chrom,VariantOfFile line)
			{
			int refId=dictionary.getSequenceIndex(chrom);
			if(refId==-1) throw new RuntimeException("unknown chromosome "+ chrom+" in "+line);
			return refId;
			}
		@Override
		public int compare(VariantOfFile o1, VariantOfFile o2)
			{
			String tokens1[]=o1.line.split("\t",5);
			String tokens2[]=o2.line.split("\t",5);
			
			int i=ref(tokens1[0],o1) - ref(tokens2[0],o2)  ;
			if(i!=0) return i;
			i= Integer.parseInt(tokens1[1])-Integer.parseInt(tokens2[1]);
			if(i!=0) return i;
			i=tokens1[3].compareTo(tokens2[3]);
			if(i!=0) return i;
			i=tokens1[4].compareTo(tokens2[4]);
			if(i!=0) return i;
			return o1.fileIndex - o2.fileIndex;
			}
		}

	
	
	@Override
	protected int doWork()
		{
		List<VCFHandler> vcfHandler=new ArrayList<VCFHandler>();
		VCFHeader h;
		try {
	    	this.dictionary=new SAMSequenceDictionaryFactory().load(REF);
	
			
			SortingCollection<VariantOfFile> array=SortingCollection.newInstance(
					VariantOfFile.class,
					new VariantCodec(),
					new VariantComparator(),
					super.MAX_RECORDS_IN_RAM
					);
			array.setDestructiveIteration(true);
			for(int fileIndex=0;fileIndex< this.IN.size();++fileIndex)
				{
				File vcfFile= this.IN.get(fileIndex);
				LOG.info("reading from "+vcfFile);
				VCFHandler handler=new VCFHandler();
				handler.vcfFile=vcfFile;
				vcfHandler.add(handler);
				
				BufferedReader in=IOUtils.openFileForBufferedReading(vcfFile);
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.startsWith("#"))
						{
						handler.headerStr.append(line).append('\n');
						continue;
						}
					VariantOfFile vof=new VariantOfFile();
					vof.fileIndex=fileIndex;
					vof.line=line;
					array.add(vof);
					}
				
				in.close();
				}
			array.doneAdding();
			
			VariantContextWriter w= VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
			CloseableIterator<VariantOfFile> iter= array.iterator();
			List<VariantOfFile> row=new ArrayList<VariantOfFile>();
			for(;;)
				{
				VariantOfFile var=null;
				if(iter.hasNext())
					{
					var=iter.next();
					}
				if(var==null || (var!=null && !row.isEmpty() && row.get(0).same(var)))
					{
					
					var=null;
					break;
					}
				else
					{
					row.add(var);
					}
				}
			w.close();
			}
		catch(Exception err)
			{
			return -1;
			}
		finally
			{
			
			}		
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		// TODO Auto-generated method stub

		}

	}
