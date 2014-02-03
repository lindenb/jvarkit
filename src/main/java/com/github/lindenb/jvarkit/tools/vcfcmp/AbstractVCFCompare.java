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

import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SortingCollection;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

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
		VCFCodec codec;
		int file_id=-1;
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
			v.line=dis.readUTF();
			return v;
			}
		@Override
		public void encode(DataOutputStream dos, LineAndFile v)
				throws IOException
			{
			dos.writeInt(v.fileIdx);
			dos.writeUTF(v.line);
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
			VariantContext ctx1=v1.getContext();
			VariantContext ctx2=v2.getContext();
			int i=ctx1.getChr().compareTo(ctx2.getChr());
			if(i!=0) return i;
			i=ctx1.getStart()-ctx2.getStart();
			if(i!=0) return i;
			i=ctx1.getReference().compareTo(ctx2.getReference());
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
	
	protected Input put(SortingCollection<LineAndFile> variants, String vcfUri)
		throws IOException
		{
		info("begin inserting "+vcfUri);
		LineIterator iter=IOUtils.openFileForLineIterator(new File(vcfUri));
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
			laf.line=line;
			VariantContext ctx=laf.getContext();
			progress.watch(ctx.getChr(), ctx.getStart());
			variants.add(laf);
			}
		progress.finish();
		info("end inserting "+vcfUri);
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
