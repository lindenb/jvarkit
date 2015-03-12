/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfcmp;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfCompareCallersOneSample
	extends AbstractVCFFilter3
	{
	private Set<File> challengerVcf=new HashSet<File>();
	private int minCountInclusive=0;
	private int maxCountInclusive=Integer.MAX_VALUE-1;
	
	
	public VcfCompareCallersOneSample()
		{
		}
	
	public Set<File> getChallengerVcf() {
		return challengerVcf;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfCompareCallersOneSample";
		}

		@Override
	public String getProgramDescription() {
		return "For my colleague Julien: VCF with one sample called using different callers. *"
				+ "Only keep variant if it was found in min<x=other-files<max";
		}
	
	public void setMinCountInclusive(int minCountInclusive) {
		this.minCountInclusive = minCountInclusive;
		}
		
	public void setMaxCountInclusive(int maxCountInclusive) {
		this.maxCountInclusive = maxCountInclusive;
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -f (file.vcf|vcf.list) can be specified multiple time. List of VCF to be challenged. "
				+ "Must be sorted on dict."
				+ "Must contain a dict. "
				+ "Only VCF containing the same sample will be considered.");
		out.println(" -m (int) min number of challengers found, inclusive. default:0");
		out.println(" -M (int) max number of challengers found, inclusive. default:unbounded");
		super.printOptions(out);
		}
	
	

	@Override
	public int executeKnime(List<String> args)
		{
		File inputFile=null;
		List<VcfIterator> listChallengers = new ArrayList<>();
		VariantContextWriter vcw=null;
		VcfIterator in=null;
		try {
			if(args.isEmpty())
				{
				in = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				inputFile=new File(args.get(0));
				in = VCFUtils.createVcfIteratorFromFile(inputFile);
				}
			else
				{
				System.err.println("#####");
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			
			
			VCFHeader header=in.getHeader();
			if(header.getNGenotypeSamples()!=1)
				{
				error(getMessageBundle("vcf.must.have.only.one.sample"));
				return -1;
				}
			
			VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

			
			
			SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict==null)
				{
				error(getMessageBundle("no.dict.in.vcf"));
				return -1;
				}
			
			/* load files to be challenged */
			for(File cf :this.challengerVcf)
				{
				//do not challenge vs itself
				if(inputFile!=null && inputFile.equals(cf))
					{
					info("Ignoring challenger (self): "+cf);
					continue;
					}
				VcfIterator cin = VCFUtils.createVcfIteratorFromFile(cf);
				VCFHeader ch=cin.getHeader();
				if(ch.getNGenotypeSamples()!=1)
					{
					warning(getMessageBundle("vcf.must.have.only.one.sample"));
					cin.close();
					continue;
					}
				if(!header.getSampleNamesInOrder().get(0).equals(
						ch.getSampleNamesInOrder().get(0)))
					{
					warning("Ignoring "+cf+" because not the same sample.");
					cin.close();
					continue;
					}
				SAMSequenceDictionary hdict = ch.getSequenceDictionary();
				if(hdict==null ||
					!SequenceUtil.areSequenceDictionariesEqual(dict, hdict))
					{
					error(getMessageBundle("not.the.same.sequence.dictionaries"));
					return -1;
					}
				listChallengers.add(cin);
				}
			
			vcw= super.createVariantContextWriter();
			vcw.writeHeader(h2);
			
			Comparator<VariantContext> ctxComparator = VCFUtils.createTidPosRefComparator(dict);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			while(in.hasNext() && !checkOutputError())
				{
				VariantContext ctx = progress.watch(in.next());
				List<VariantContext> ctxChallenging= new ArrayList<>();
				
				int countInOtherFiles=0;
				int idx=0;
				while( idx < listChallengers.size())
					{
					VcfIterator citer = listChallengers.get(idx);
					if(!citer.hasNext())
						{
						citer.close();
						listChallengers.remove(idx);
						continue;
						}
					boolean foundInThatFile=false;
					int diff = 0;
					while(diff<=0 && citer.hasNext())
						{
						diff=ctxComparator.compare(
							citer.peek(),
							ctx
							);
						if(diff < 0)
							{
							//consumme & ignore
							citer.next(); 
							}
						else if(diff==0)
							{
							foundInThatFile=true;
							ctxChallenging.add(citer.next());
							}
						else
							{
							break;
							}
						}
					countInOtherFiles+=(foundInThatFile?1:0);
					++idx;
					}
				
				if(countInOtherFiles >= minCountInclusive &&
					countInOtherFiles <= maxCountInclusive)
					{
					this.incrVariantCount();
					vcw.add(ctx);
					}
				
				}
			progress.finish();
			return 0;
			} 
		catch (Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcw);
			CloserUtil.close(listChallengers);
			CloserUtil.close(in);
			}
		
		}
	
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"f:o:m:M:"))!=-1)
			{
			switch(c)
				{
				case 'm': setMinCountInclusive(Integer.parseInt(opt.getOptArg()));break;
				case 'M': setMaxCountInclusive(Integer.parseInt(opt.getOptArg()));break;
				case 'f':
					File f=new File(opt.getOptArg());
					if(f.getName().endsWith(".list"))
						{
						try {
							for(String L: IOUtil.readLines(f))
								{
								if(L.trim().isEmpty() || L.startsWith("#")) continue;
								this.getChallengerVcf().add(new File(L.trim()));
								}
						} catch (Exception e) {
							error(e);
							return -1;
							}
						}
					else
						{
						this.getChallengerVcf().add(f);
						}
					break;
				case 'o': this.setOutputFile(opt.getOptArg()); break;
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
		return mainWork(opt.getOptInd(), args);
		}
	public static void main(String[] args) {
		new VcfCompareCallersOneSample().instanceMainWithExit(args);
	}
	}
