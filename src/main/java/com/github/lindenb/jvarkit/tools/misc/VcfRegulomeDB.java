package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.tabix.AbstractTabixObjectReader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * VcfRegulomeDB
 * @author lindenb
 *
 */
public class VcfRegulomeDB extends AbstractVcfRegulomeDB
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfRegulomeDB.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractVcfRegulomeDB.AbstractVcfRegulomeDBCommand
		{		
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
		@Override
	protected Collection<Throwable> doVcfToVcf(String inputName,
			VcfIterator in, VariantContextWriter out) throws IOException
		{
		VCFHeader header=in.getHeader();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		addMetaData(header);
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
		return RETURN_OK;
		}
	@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			
			if(this.bedFile==null)
				{
				return wrapException("Bed file indexed with tabix is missing");
				}
			if(this.acceptRegexStr!=null)
				{
				this.acceptRegex = Pattern.compile(this.acceptRegexStr);
				}
			try
				{
				LOG.info("Opening "+bedFile);
				this.regDataTabixFileReader=new RegDataTabixFileReader(bedFile);
				return doVcfToVcf(inputName);
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(this.regDataTabixFileReader);
				this.acceptRegex=null;
				this.regDataTabixFileReader=null;
				}
			}
		}
	public static void main(String[] args)
		{
		new VcfRegulomeDB().instanceMainWithExit(args);
		}
	}
