package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;


public class VcfIn extends AbstractVCFCompare
	{
	
	private static final int DATABASE_VCF=0;
	private static final int USER_VCF=1;
	private VcfIn()
		{
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfIn";
		}

		@Override
	public String getProgramDescription() {
		return "Only prints variants that are contained/not contained into another VCF.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -v : inverse. Print variant that are not part of the VCF-database.");
		super.printOptions(out);
		}
	
	
	private final Pattern tab=Pattern.compile("[\t]");
	/* for file0 we only have to keep chrom/start/end */
	@Override
	protected String simplify(String line, int input_idx)
		{
		if(input_idx==DATABASE_VCF)
			{
			String tokens[]=tab.split(line,5);
			return tokens[0]+"\t"+tokens[1]+"\t.\t"+tokens[3]+"\t.";
			}
		return line;
		}
	
	@Override
	public int doWork(String[] args)
		{
		boolean print_in_both_file=true;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"v"))!=-1)
			{
			switch(c)
				{
				case 'v': print_in_both_file=false;break;
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
		
		SortingCollection<LineAndFile> variants=null;
		VariantContextWriter vcw=null;
		try
			{
			final LineAndFileComparator posCompare=new LineAndFileComparator();

			
			factory.setComponentType(LineAndFile.class);
			factory.setComparator(posCompare);
			factory.setTmpDirs(this.getTmpDirectories());
			factory.setCodec(new LineAndFileCodec());
			variants=this.factory.make();
			variants.setDestructiveIteration(true);
			Input inputs[]=new Input[2];
			if(opt.getOptInd()+1 == args.length)
				{
				inputs[DATABASE_VCF]=super.put(variants, args[opt.getOptInd()]);
				inputs[USER_VCF]=super.put(variants, null);
				}
			else if(opt.getOptInd()+2 == args.length)
				{
				inputs[DATABASE_VCF]=super.put(variants, args[opt.getOptInd()+0]);
				inputs[USER_VCF]=super.put(variants, args[opt.getOptInd()+1]);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			VCFHeader header=super.inputs.get(1).header;
			vcw = VCFUtils.createVariantContextWriterToStdout();
			
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			vcw.writeHeader(header);
			
			List<LineAndFile> row=new ArrayList<LineAndFile>(super.inputs.size());
			
			
			CloseableIterator<LineAndFile> iter=variants.iterator();
			int n_printed=0;
			for(;;)
				{
				LineAndFile rec=null;
				if(iter.hasNext())
					{
					rec=iter.next();
					}
			
				if(rec==null || (!row.isEmpty() && posCompare.compare(row.get(0),rec)!=0))
					{
					if(!row.isEmpty())
						{
						boolean found[]=new boolean[]{false,false};
						for(LineAndFile laf:row)
							{
							found[laf.fileIdx]=true;
							}
						
					
						
						if(found[USER_VCF] && ((print_in_both_file && found[DATABASE_VCF]) || ((!print_in_both_file && !found[DATABASE_VCF]))))
							{
							for(LineAndFile laf:row)
								{
								if(laf.fileIdx==USER_VCF)
									{
									vcw.add(laf.getContext());
									n_printed++;
									}
								}
							}
						row.clear();
						}
					
					if(rec==null) break;
					}
				row.add(rec);
				}
			CloserUtil.close(iter);
			info("Done. N="+n_printed);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcw);
			try
				{
				if(variants!=null) variants.cleanup();
				}
			catch(Exception err)
				{
				}
			}
		}
	public static void main(String[] args) {
		new VcfIn().instanceMainWithExit(args);
	}
	}
