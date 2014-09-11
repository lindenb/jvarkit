package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public abstract class AbstractVCFCompare extends AbstractCommandLineProgram
	{
	protected List<Input> inputs=new ArrayList<Input>();
	protected SortingCollectionFactory<LineAndFile> factory=new SortingCollectionFactory<LineAndFile>();

	protected class Input
		{
		String filename;
		VCFHeader header;
		AbstractVCFCodec codec;
		int file_id=-1;
		long count=0L;
		}
	
	protected class LineAndFile
		{
		int fileIdx;
		String line;
		
		private VariantContext _ctx=null;
		VariantContext getContext()
			{
			if(this._ctx==null)
				{
				this._ctx=inputs.get(this.fileIdx).codec.decode(this.line);
				}
			return this._ctx;
			}
		private String token(int index)
			{
			int i=0;
			int prev=0;
			while(prev< line.length())
				{
				int tab=this.line.indexOf('\t',prev);
				if(tab==-1)
					{
					tab=this.line.length();
					}
				if(i==index) return this.line.substring(prev, tab);
				prev=tab+1;
				++i;
				}
			throw new IllegalStateException("Bad VCF line :"+line+" cannot get tokens["+index+"]");
			}
		
		public String getChrom()
			{	
			return token(0);
			}
		public int getStart()
			{	
			return Integer.parseInt(token(1));
			}
		public String getReference()
			{	
			return token(3);
			}
		}
	
	protected class LineAndFileCodec extends AbstractDataCodec<LineAndFile>
		{
		@Override
		public LineAndFile decode(DataInputStream dis) throws IOException
			{
			LineAndFile v=new LineAndFile();
			try {
				v.fileIdx=dis.readInt();
			} catch (Exception e) {
				return null;
				}
			v.line= readString(dis);
			return v;
			}
		@Override
		public void encode(DataOutputStream dos, LineAndFile v)
				throws IOException
			{
			dos.writeInt(v.fileIdx);
			writeString(dos,v.line);
			}
		@Override
		public AbstractDataCodec<LineAndFile> clone() {
			return new LineAndFileCodec();
			}
		}

	
	
	
	
	
	protected class LineAndFileComparator implements Comparator<LineAndFile>
		{
		@Override
		public int compare(LineAndFile v1, LineAndFile v2)
			{
			int i=v1.getChrom().compareTo(v2.getChrom());
			if(i!=0) return i;
			i=v1.getStart()-v2.getStart();
			if(i!=0) return i;
			i=v1.getReference().compareToIgnoreCase(v2.getReference());
			if(i!=0) return i;
			return 0;
			}
		}
	
	
	
	protected AbstractVCFCompare()
		{
		}
	
	
	@Override
	public String getProgramDescription()
		{
		return "Compares two VCF files.";
		}
	
	@Override
	public void printOptions(PrintStream out) {
		out.println(" -M (int) Max records in RAM. Optional.");
		out.println(" -T (dir) add temporary directory. Optional");
		super.printOptions(out);
		}
	
	@Override
	protected String getGetOptDefault() {
		return super.getGetOptDefault()+"M:T:";
		}
	
	protected Comparator<LineAndFile> createLineAndFileComparator()
		{
		return new LineAndFileComparator();
		}

	private static class ListStringLineReader implements LineReader
		{
		private List<String> L;
		ListStringLineReader(List<String> L)
			{
			this.L=L;
			}
		@Override
		public void close()
			{
			this.L=null;
			}
		@Override
		public String readLine() throws IOException
			{
			return L==null || L.isEmpty()?null:L.remove(0);
			}
		}
	
	protected Input createInput()
		{
		return new Input();
		}
	
	private Input createInput(String filename,List<String> headerLines)
		{
		Input input=createInput();
		input.filename=filename;
		StringWriter sw=new StringWriter();
		for(String s:headerLines) sw.append(s).append("\n");
		LineIteratorImpl li=new LineIteratorImpl(new ListStringLineReader(headerLines));
		VcfIterator iter=new VcfIterator(li);
		input.header=iter.getHeader();
		input.codec=iter.getCodec();
		CloserUtil.close(iter);
		return input;
		}
	
	protected String simplify(String line,int input_idx)
		{
		return line;
		}
	
	protected Input put(SortingCollection<LineAndFile> variants, String vcfUri)
		throws IOException
		{
		info("begin inserting "+vcfUri);
		LineIterator iter=null;
		if(vcfUri==null)
			{
			vcfUri="stdin";
			iter=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(System.in));
			}
		else
			{
			iter=IOUtils.openFileForLineIterator(new File(vcfUri));
			}
		
		List<String> headerLines=new ArrayList<String>();
		while(iter.hasNext() && iter.peek().startsWith("#"))
			{
			headerLines.add(iter.next());
			}
		if(headerLines.isEmpty()) throw new IOException("Not header found in "+vcfUri);
		Input input=createInput(vcfUri, headerLines);
		input.file_id=this.inputs.size();
		this.inputs.add(input);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(input.header.getSequenceDictionary());
		while(iter.hasNext())
			{
			String line=iter.next();
			LineAndFile laf=new LineAndFile();
			laf.fileIdx=input.file_id;
			laf.line=simplify(line,laf.fileIdx);
			progress.watch(laf.getChrom(), laf.getStart());
			variants.add(laf);
			input.count++;
			}
		progress.finish();
		info("end inserting "+vcfUri+" N="+input.count);
		CloserUtil.close(iter);
		return input;
		}
	
	@Override
	protected GetOptStatus
		handleOtherOptions(int c, GetOpt opt, String[] args)
		{
		switch(c)
			{
			case 'M': factory.setMaxRecordsInRAM(Math.max(1,Integer.parseInt(opt.getOptArg()))); return GetOptStatus.OK;
			case 'T': super.addTmpDirectory(new File(opt.getOptArg())) ; return GetOptStatus.OK;
			default: return super.handleOtherOptions(c, opt, args);
			}
		}
	
	}
