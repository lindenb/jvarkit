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
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfCompareCallers
	extends AbstractKnimeApplication
	{
	private enum Category
		{
		unique_to_file_1,
		unique_to_file_1_snp,
		unique_to_file_1_indel,
		unique_to_file_2,
		unique_to_file_2_snp,
		unique_to_file_2_indel,
		common_context,
		common_context_snp,
		common_context_indel,
		common_context_discordant_id,
		called_and_same,
		called_and_same_hom_ref,
		called_and_same_hom_var,
		called_and_same_het,
		called_but_discordant,
		called_but_discordant_hom1_het2,
		called_but_discordant_het1_hom2,
		called_but_discordant_hom1_hom2,
		called_but_discordant_het1_het2,
		
		}
	
	public VcfCompareCallers()
		{
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"VcfCompareCallers";
		}

		@Override
	public String getProgramDescription() {
		return "Compare two VCFs; Prints a list table f(Y=sample, X=category)";
		}
	
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (filename) output. Default:stdout.");
		super.printOptions(out);
		}
	

	@SuppressWarnings("resource")
	@Override
	public int executeKnime(List<String> args)
		{
		PrintStream pw=null;
		VcfIterator vcfInputs[]=new VcfIterator[]{null,null};
		VCFHeader headers[]=new VCFHeader[]{null,null};
		try {
			if(args.size()==1)
				{
				info("Reading from stdin and "+ args.get(0));
				vcfInputs[0] = VCFUtils.createVcfIteratorStdin();
				vcfInputs[1] = VCFUtils.createVcfIterator( args.get(0));
				}
			else if(args.size()==2)
				{
				info("Reading from stdin and "+ args.get(0)+" and "+ args.get(1));
				vcfInputs[0] = VCFUtils.createVcfIterator( args.get(0));
				vcfInputs[1] = VCFUtils.createVcfIterator( args.get(1));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
		
			for(int i=0;i< vcfInputs.length;++i)
				{
				headers[i] = vcfInputs[i].getHeader();
				}
			/* dicts */
			SAMSequenceDictionary dict0 = headers[0].getSequenceDictionary();
			SAMSequenceDictionary dict1 = headers[1].getSequenceDictionary();
			Comparator<VariantContext> ctxComparator=null;
			if(dict0==null && dict1==null)
				{
				ctxComparator = VCFUtils.createChromPosRefComparator();
				}
			else if(dict0!=null && dict1!=null)
				{
				if( !SequenceUtil.areSequenceDictionariesEqual(dict0, dict1))
					{
					error(getMessageBundle("not.the.same.sequence.dictionaries"));
					return -1;
					}
				ctxComparator = VCFUtils.createTidPosRefComparator(dict0);
				}
			else
				{
				error(getMessageBundle("not.the.same.sequence.dictionaries"));
				return -1;
				}
			/* samples */
			Set<String> samples0=new HashSet<>(headers[0].getSampleNamesInOrder());
			Set<String> samples1=new HashSet<>(headers[1].getSampleNamesInOrder());
			Set<String> samples= new TreeSet<>(samples0);
			samples.retainAll(samples1);
			
			if(samples.size()!=samples0.size() || samples.size()!=samples1.size())
				{
				warning("Warning: Not the same samples set. Using intersection of both lists.");
				}
			Map<String, Counter<Category>> sample2info=new HashMap<String, Counter<Category>>(samples.size());
			for(String sampleName:samples)
				{
				sample2info.put(sampleName, new  Counter<Category>());
				}
			
			SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(dict0);
			VariantContext buffer[]=new VariantContext[vcfInputs.length];
			VariantContext prev[]=new VariantContext[vcfInputs.length];
			for(;;)
				{
				VariantContext smallest=null;
				//refill buffer
				for(int i=0;i< vcfInputs.length;++i)
					{
					if(buffer[i]==null && vcfInputs[i]!=null )
						{
						if(vcfInputs[i].hasNext())
							{
							buffer[i]=vcfInputs[i].peek();
							/* check data are sorted */
							if(prev[i]!=null && ctxComparator.compare(prev[i], buffer[i])>0)
								{
								error("Input "+(i+1)+"/2 is not sorted"+(
									((i==0 && dict0==null) ||(i==1 && dict1==null))?
									"on chrom/pos/ref":
									"on sequence dictionary"
									)+". got\n"+buffer[i]+"\nafter\n"+prev[i]);
								return -1;
								}
							}
						else
							{
							vcfInputs[i].close();
							vcfInputs[i]=null;
							}							
						}
					
					if(buffer[i]!=null )
						{
						if(smallest==null || ctxComparator.compare(buffer[i],smallest)<0)
							{
							smallest=buffer[i];
							}
						}
					}
				
				if(smallest==null) break;
				
				VariantContext ctx0=null;
				VariantContext ctx1=null;
				
				if(buffer[0]!=null && ctxComparator.compare(buffer[0],smallest)==0)
					{
					prev[0] = progress.watch(vcfInputs[0].next());
					ctx0= prev[0];
					buffer[0]=null;
					}
				if(buffer[1]!=null && ctxComparator.compare(buffer[1],smallest)==0)
					{
					prev[1]= progress.watch(vcfInputs[1].next());
					ctx1= prev[1];
					buffer[1]=null;
					}
	
				
				for(String sampleName: sample2info.keySet())
					{
					Counter<Category> sampleInfo=sample2info.get(sampleName);
					Genotype g0=(ctx0==null?null:ctx0.getGenotype(sampleName));
					Genotype g1=(ctx1==null?null:ctx1.getGenotype(sampleName));
					if(g0!=null && (g0.isNoCall() || !g0.isAvailable())) g0=null;
					if(g1!=null && (g1.isNoCall() || !g1.isAvailable())) g1=null;
					
					if(g0==null && g1==null) continue;
					if(g0!=null && g1==null)
						{
						sampleInfo.incr(Category.unique_to_file_1);
						if(ctx0.isIndel())
							{
							sampleInfo.incr(Category.unique_to_file_1_indel);
							}
						else if(ctx0.isSNP())
							{
							sampleInfo.incr(Category.unique_to_file_1_snp);
							}
						continue;
						}
					else if(g0==null && g1!=null)
						{
						sampleInfo.incr(Category.unique_to_file_2);
						if(ctx1.isIndel())
							{
							sampleInfo.incr(Category.unique_to_file_2_indel);
							}
						else if(ctx1.isSNP())
							{
							sampleInfo.incr(Category.unique_to_file_2_snp);
							}
						continue;
						}
					else
						{	
						sampleInfo.incr(Category.common_context);
						if(ctx0.isIndel() && ctx1.isIndel())
							{
							sampleInfo.incr(Category.common_context_indel);
							}
						else if(ctx0.isSNP() && ctx1.isSNP())
							{
							sampleInfo.incr(Category.common_context_snp);
							}
						
						if( (ctx0.hasID() && !ctx1.hasID()) ||
							(!ctx0.hasID() && ctx1.hasID()) ||
							(ctx0.hasID() && ctx1.hasID() && !ctx0.getID().equals(ctx1.getID()))
							)
							{
							sampleInfo.incr(Category.common_context_discordant_id);
							}
						
						
						if(g0.sameGenotype(g1))
							{
							sampleInfo.incr(Category.called_and_same);
							if(g0.isHomRef())
								{
								sampleInfo.incr(Category.called_and_same_hom_ref);
								}
							if(g0.isHomVar())
								{
								sampleInfo.incr(Category.called_and_same_hom_var);
								}
							else if(g0.isHet())
								{	
								sampleInfo.incr(Category.called_and_same_het);
								}
							
							}
						else
							{
							sampleInfo.incr(Category.called_but_discordant);
							if(g0.isHom() && g1.isHet())
								{
								sampleInfo.incr(Category.called_but_discordant_hom1_het2);
								}
							else if(g0.isHet() && g1.isHom())
								{
								sampleInfo.incr(Category.called_but_discordant_het1_hom2);
								}
							else if(g0.isHom() && g1.isHom())
								{
								sampleInfo.incr(Category.called_but_discordant_hom1_hom2);
								}
							else if(g0.isHet() && g1.isHet())
								{
								sampleInfo.incr(Category.called_but_discordant_het1_het2);
								}
							}
						
						}
					}
				}
			progress.finish();
		
			pw  = (getOutputFile()==null?System.out:new PrintStream(getOutputFile()));
			pw.print("#Sample");
			for(Category c:Category.values())
				{
				pw.print('\t');
				pw.print(c.name());
				}
			pw.println();
			for(String sample: sample2info.keySet())
				{
				Counter<Category> count=sample2info.get(sample);
				pw.print(sample);
				for(Category c:Category.values())
					{
					pw.print('\t');
					pw.print(count.count(c));
					}
				pw.println();
				if(pw.checkError()) break;
				}
			pw.flush();
			
			
			return 0;
			} 
		catch (Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			if(getOutputFile()!=null)  CloserUtil.close(pw);
			}
		
		}

	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
				{
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
		new VcfCompareCallers().instanceMainWithExit(args);
	}
	}
