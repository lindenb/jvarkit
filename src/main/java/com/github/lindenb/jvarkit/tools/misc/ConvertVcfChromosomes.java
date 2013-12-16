package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class ConvertVcfChromosomes extends AbstractVCFFilter2 {
	private boolean use_original_chrom_name_if_no_mapping=false;
	private Map<String,String> customMapping=new HashMap<String,String>();
	private Set<String> unmappedChromosomes=new HashSet<String>();
	private boolean ignore_if_no_mapping=false;
	private ConvertVcfChromosomes()
		{
		
		}
	
	private String convertName(String chrom)throws IOException
		{
		if(chrom==null) throw new NullPointerException();
		String newname=customMapping.get(chrom);
		if(newname==null)
			{
			if(!unmappedChromosomes.contains(chrom))
				{
				warning("unmapped chromosome "+chrom);
				unmappedChromosomes.add(chrom);
				}
			if(ignore_if_no_mapping) return null;
			
			if(use_original_chrom_name_if_no_mapping)
				{	
				return chrom;
				}
			throw new IOException("No mapping found to convert name of chromosome \""+chrom+"\"");
			}
		return newname;
		}
	
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException {
		VCFHeader header1=in.getHeader();
	
		Set<VCFHeaderLine> meta2=new LinkedHashSet<VCFHeaderLine>();
		for(VCFHeaderLine L:header1.getMetaDataInInputOrder())
			{
			if(!L.getKey().equals(VCFHeader.CONTIG_KEY))
				{
				meta2.add(L);
				}
			}
		VCFHeader header2=new VCFHeader(meta2,header1.getSampleNamesInOrder());

		if(header1.getSequenceDictionary()!=null)
			{
			List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>();

			for(int i=0;i< header1.getSequenceDictionary().size();++i)
				{
				SAMSequenceRecord ssr=header1.getSequenceDictionary().getSequence(i);
				String newName=convertName(ssr.getSequenceName());
				if(newName==null)
					{
					//skip unknown chromosomes
					continue;
					}
				ssr=new SAMSequenceRecord(newName, ssr.getSequenceLength());
				ssrs.add(ssr);
				}
			header2.setSequenceDictionary(new SAMSequenceDictionary(ssrs));
			}
		
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		out.writeHeader(header2);
		
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			String newName=convertName(ctx.getChr());
			if(newName==null)
				{
				//skip unknown chromosomes
				continue;
				}
			vcb.chr(newName);
			ctx=vcb.make();

			out.add(ctx);
			}
		
		if(!unmappedChromosomes.isEmpty())
			{
			warning("Unmapped chromosomes: "+unmappedChromosomes);
			}
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfRenameChromosomes";
		}
	
	@Override
	public String getProgramDescription() {
		return "Convert the names of the chromosomes in a VCF file.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -f (file) load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+");
		out.println(" -i if no mapping found, skip that record.");
		out.println(" -C if no mapping found, use the original name instead of throwing an error. ");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"f:Ci"))!=-1)
			{
			switch(c)
				{
				case 'i': ignore_if_no_mapping=true;break;
				case 'C': use_original_chrom_name_if_no_mapping=true;break;
				case 'f':
					{
					File f=new File(opt.getOptArg());
					BufferedReader in=null;
					try
						{
						info("Loading custom mapping "+f);
						in=IOUtils.openFileForBufferedReading(f);
						String line;
						while((line=in.readLine())!=null)
							{
							if(line.isEmpty() || line.startsWith("#")) continue;
							String tokens[]=line.split("[\t]");
							if(tokens.length!=2
									|| tokens[0].trim().isEmpty()
									|| tokens[1].trim().isEmpty()
									) throw new IOException("Bad mapping line: \""+line+"\"");
							tokens[0]=tokens[0].trim();
							tokens[1]=tokens[1].trim();
							if(customMapping.containsKey(tokens[0]))
								{
								throw new IOException("Mapping defined twice for: \""+tokens[0]+"\"");
								}
							customMapping.put(tokens[0], tokens[1]);
							}
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					finally
						{
						CloserUtil.close(in);
						}
					break;
					}
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		
		return super.doWork(opt.getOptInd(), args);
		}
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new ConvertVcfChromosomes().instanceMainWithExit(args);
		}
	}
