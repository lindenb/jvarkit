
package com.github.lindenb.jvarkit.tools.vcfconcat;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfConcat extends AbstractKnimeApplication
	{
	public static final String VARIANTSOURCE="VARIANTSOURCE";
	private Set<String> inputFiles=new HashSet<String>();
	private File outputfile=null;
	
	
	public VcfConcat()
		{
		}
	
	@Override
	public void setOutputFile(File out)
		{
		this.outputfile=out;
		}
	public File getOutputFile()
		{
		return outputfile;
		}
	
	public Set<String> getInputFiles()
		{
		return inputFiles;
		}

	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"VcfConcat";
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Concatenate VCF files.";
		}
	
	private int compareHeader(VCFHeader h1,VCFHeader h2)
		{
		if(!h1.getSampleNameToOffset().equals(h2.getSampleNameToOffset()))
			{
			error("Not the same sample list to concatenate:"+
					h1.getSampleNameToOffset()+" vs "+
					h2.getSampleNameToOffset()
					);
			return -1;
			}
		
		return 0;
		}
	
	protected VariantContextWriter createVariantContextWriter()
			throws IOException
			{
			if(getOutputFile()==null)
				{
				return VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				info("opening vcf writer to "+getOutputFile());
				return VCFUtils.createVariantContextWriter(getOutputFile());
				}

			}

	
	private int fromFiles() throws IOException
		{
		VariantContextWriter out=null;
		List<VcfIterator> inputs=new ArrayList<VcfIterator>(getInputFiles().size());
		List<String> inputFiles=new ArrayList<>(getInputFiles().size());
		List<String> samples=new ArrayList<>();
		SAMSequenceDictionary dict=null;
		try
			{
			Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
			
			for(String vcfFile:this.getInputFiles())
				{
				info("Opening "+vcfFile);
				VcfIterator r=VCFUtils.createVcfIterator(vcfFile);
				
				
				VCFHeader header = r.getHeader();
				if(dict==null && inputs.isEmpty())
					{
					dict = header.getSequenceDictionary();
					}
				else if(header.getSequenceDictionary()==null && dict!=null)
					{
					error(getMessageBundle("no.dict.in.vcf"));
					return -1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, header.getSequenceDictionary()))
					{
					error(getMessageBundle("not.the.same.sequence.dictionaries"));
					return -1;
					}
				if(inputs.isEmpty())
					{
					samples = header.getSampleNamesInOrder();
					}
				if(inputs.size()>0 && !header.getSampleNamesInOrder().equals(samples))
					{
					error("No same samples");
					return -1;
					}
				
				metaData.addAll(header.getMetaDataInInputOrder());
				inputs.add(r);
				inputFiles.add(VCFUtils.escapeInfoField(vcfFile));
				}
			final Comparator<VariantContext> comparator = (dict==null?
					VCFUtils.createChromPosRefComparator():
					VCFUtils.createTidPosRefComparator(dict)
					);
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			metaData.add(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
			metaData.add(new VCFInfoHeaderLine(
					VARIANTSOURCE,1,VCFHeaderLineType.String,
					"Origin File of Varant"
					));
			VCFHeader h2 = new VCFHeader(
					metaData,
					samples
					);
			out= createVariantContextWriter();
			out.writeHeader(h2);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			for(;;)
				{
				/* get 'smallest' variant */
				VariantContext smallest = null;
				int idx=0;
				int best_idx = -1;
				
				while(idx < inputs.size())
					{
					VcfIterator in= inputs.get(idx);
					if(!in.hasNext())
						{
						CloserUtil.close(in);
						inputs.remove(idx);
						inputFiles.remove(idx);
						}
					else
						{
						VariantContext ctx = in.peek();
						if( smallest==null ||
							comparator.compare(smallest,ctx)>0)
							{
							smallest = ctx;
							best_idx = idx;
							}
						++idx;
						}
					}
				
				if(smallest==null) break;
				
				VariantContext ctx = progress.watch(inputs.get(best_idx).next());
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(VARIANTSOURCE, inputFiles.get(best_idx));
				out.add(vcb.make());
				}
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(inputs);
			CloserUtil.close(out);
			}

		}
	
	@Override
	public int doWork(String[] args)
		{
		GetOpt getopt=new GetOpt();
		int c;
		while((c=getopt.getopt(args, super.getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default: switch(super.handleOtherOptions(c, getopt, args))
					{
					case EXIT_FAILURE: return -1;
					case EXIT_SUCCESS: return 0;
					default:break;
					}
				}
			}
				
		VariantContextWriter w=null;
		try
			{
			w=VCFUtils.createVariantContextWriter(null);
			if(getopt.getOptInd()==args.length)
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			else
				{
				for(int i=0;i< getopt.getOptInd();++i)
					{
					String filename=args[i];
					getInputFiles().addAll(IOUtils.unrollFiles(Arrays.asList(filename)));	
					}
				return fromFiles();
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfConcat().instanceMainWithExit(args);
		}

	}
