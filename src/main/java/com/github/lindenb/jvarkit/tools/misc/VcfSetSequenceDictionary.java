package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfSetSequenceDictionary extends AbstractVCFFilter2
	{
	private SAMSequenceDictionary dict=null;
	private LinkedHashMap<String, Integer> newdict=null;
	
	private VcfSetSequenceDictionary()
		{
		}
	

	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		Set<VCFHeaderLine> meta2=new LinkedHashSet<VCFHeaderLine>();
		for(VCFHeaderLine L:header.getMetaDataInInputOrder())
			{
			if(!L.getKey().equals(VCFHeader.CONTIG_KEY))
				{
				meta2.add(L);
				}
			}
		
		meta2.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		meta2.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		
		if(dict!=null)
			{	
			meta2.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(dict));
			}
		else
			{
			warning("No sequence dictionary was defined");
			}
		VCFHeader header2=new VCFHeader(meta2, header.getSampleNamesInOrder());

		out.writeHeader(header2);
		
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			if(dict!=null && dict.getSequenceIndex(ctx.getContig())==-1)
				{
				warning("Unknown chromosome "+ctx.getContig());
				}
			if(newdict!=null)
				{
				Integer length=this.newdict.get(ctx.getContig());
				if(length==null) length=0;
				if(ctx.getEnd()>length)
					{
					this.newdict.put(ctx.getContig(),ctx.getEnd());
					}
				}
			
			
			out.add(ctx);
			}
			
		}
	
	@Override
	public String getProgramDescription() {
		return "Set the ##contig lines in a VCF header";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfSetSequenceDictionary";
		}

	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -r (fasta) indexed reference. Optional.");
		out.println(" -d (file.dict) at the end, save an alternate dict in that file. Optional.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File newDictOut=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "r:d:"))!=-1)
			{
			switch(c)
				{
				case 'd':
					{
					newDictOut=new File(opt.getOptArg());
					if(!newDictOut.getName().endsWith(".dict"))
						{
						error("dictionary should end with .dict :"+newDictOut);
						return -1;
						}
					this.newdict=new LinkedHashMap<String, Integer>();
					break;
					}
				case 'r':
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
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		int err= doWork(opt.getOptInd(), args);
		
		if(newDictOut!=null)
			{
			info("Saving alt dictionary "+newDictOut);
			FileWriter out=null;
			try
				{
				
				List<SAMSequenceRecord> list=new ArrayList<SAMSequenceRecord>(newdict.size());
				for(String k:this.newdict.keySet())
					{
					list.add(new SAMSequenceRecord(k, this.newdict.get(k)));
					}
				SAMFileHeader sfh=new SAMFileHeader();
				sfh.setSequenceDictionary(
						new SAMSequenceDictionary(list)
						);
		        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
		        codec.setValidationStringency(htsjdk.samtools.ValidationStringency.SILENT);
				out=new FileWriter(newDictOut);
				codec.encode(out, sfh);
				out.flush();
				}
			catch(Exception err2)
				{
				error(err2);
				return -1;
				}
			finally
				{
				CloserUtil.close(out);
				}
			}
		
		return err;
		}

	public static void main(String[] args)
		{
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
		}
	}
