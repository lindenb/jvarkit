package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * VcfRegulomeDB
 * @author lindenb
 *
 */
public class VcfRegulomeDB extends AbstractVCFFilter2
	{
	private String bedFile=null;
	private String infoTag="REGULOMEDB";
	private RegDataTabixFileReader regDataTabixFileReader;
	private int extend=5;
	private Pattern acceptRegex=null;
	
	private static class RegData
		{
		@SuppressWarnings("unused")
		String chrom;
		int chromSart;
		@SuppressWarnings("unused")
		int chromEnd;
		String rank;
		RegData(String tokens[])
			{
			this.chrom=tokens[0];
			this.chromSart=Integer.parseInt(tokens[1]);
			this.chromEnd=Integer.parseInt(tokens[2]);
			this.rank=tokens[3];
			}
		}
	
	private static class RegDataTabixFileReader
		extends AbstractTabixObjectReader<RegData>
		{
		RegDataTabixFileReader(String uri) throws IOException
			{
			super(uri);
			}
		@Override
		protected Iterator<RegData> iterator(Iterator<String> delegate) {
			return new MyIterator(delegate);
			}
		private class MyIterator
	    	extends AbstractMyIterator
	    	{
	    	private Pattern tab=Pattern.compile("[\t]");
	    	MyIterator(Iterator<String> delegate)
	    		{
	    		super(delegate);
	    		}
	    	
	    	@Override
	    	public RegData next() {
	    		return new RegData(this.tab.split(delegate.next(),5));
	    		}
	    	}

		}
	
	private VcfRegulomeDB()
		{
		
		}
	
	@Override
	public String getProgramDescription() {
		return "Annotate a VCF with the Regulome data (http://regulome.stanford.edu/)";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfRegulomeDB";
		}
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		header.addMetaDataLine(new VCFInfoHeaderLine(
				this.infoTag,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"Format: Position|Distance|Rank"
				));
		
		
		out.writeHeader(header);
		
		while(in.hasNext())
			{
			List<String>  regDataList=new ArrayList<String>();
			VariantContext ctx=in.next();
			
			progress.watch(ctx.getContig(),ctx.getStart());
			
			int start=Math.max(0,ctx.getStart()-this.extend);
			int end=ctx.getEnd()+this.extend;
			
			for(Iterator<RegData> iter=this.regDataTabixFileReader.iterator(ctx.getContig(), start, end);
					iter.hasNext();
					)
				{
				RegData curr=iter.next();
				if(this.acceptRegex!=null && 
				   !this.acceptRegex.matcher(curr.rank).matches()
				   )
					{
					continue;
					}
				String str=
						String.valueOf(curr.chromSart)+"|"+
						String.valueOf(Math.abs(curr.chromSart-(ctx.getStart()-1)))+"|"+
						curr.rank;
				regDataList.add(str);
				}
			if(regDataList.isEmpty())
				{
				out.add(ctx);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.attribute(this.infoTag, regDataList.toArray());
			out.add(vcb.make());
			}
		progress.finish();
		}
		@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -b (file) bed indexed with tabix. Format: chrom(tab)start(tab)end(tab)rank Required.");
		out.println(" -T (string) tag in vcf INFO. Default:"+this.infoTag);
		out.println(" -x (int) base pairs. look.for data around the variation +/- 'x' Default:"+this.extend);
		out.println(" -r (regex-pattern) if defined, only accept the rank matching the regular expression. Optional.");
		super.printOptions(out);
		}
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"b:T:x:r:"))!=-1)
			{
			switch(c)
				{
				case 'b': bedFile= opt.getOptArg();break;
				case 'T': infoTag= opt.getOptArg();break;
				case 'x': this.extend= Integer.parseInt(opt.getOptArg());break;
				case 'r': this.acceptRegex=Pattern.compile(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(bedFile==null)
			{
			error("Bed file indexed with tabix is missing");
			return -1;
			}
		try
			{
			info("Opening "+bedFile);
			this.regDataTabixFileReader=new RegDataTabixFileReader(bedFile);
			return super.doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.regDataTabixFileReader);
			}
		}
	public static void main(String[] args)
		{
		new VcfRegulomeDB().instanceMainWithExit(args);
		}
	}
