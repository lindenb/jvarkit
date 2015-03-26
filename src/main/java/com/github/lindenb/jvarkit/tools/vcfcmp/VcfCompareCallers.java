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
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

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
		both_missing,
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
		called_but_discordant_others,
		
		}
	private int numberOfExampleVariants=10;
	private File exampleFile=null;
	
	
	public VcfCompareCallers()
		{
		}
	
	
	public void setExampleFile(File exampleFile) {
		this.exampleFile = exampleFile;
		}
	
	public void setNumberOfExampleVariants(int numberOfExampleVariants) {
		this.numberOfExampleVariants = numberOfExampleVariants;
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
		out.println(" -n number of variants to dump in the example file. Default: "+numberOfExampleVariants);
		out.println(" -e (example file.xml). Print a few representative Variants in this XML file. Optional.");
		super.printOptions(out);
		}
	
	private void watch(
			XMLStreamWriter out,
			VariantContext ctx0,
			VariantContext ctx1,
			Genotype g0,
			Genotype g1,
			String sampleName,Counter<Category> count,Category cat
			) throws XMLStreamException
		{
		long n=count.incr(cat);
		if(out==null || n> this.numberOfExampleVariants) return;
		VariantContext variants[]=new VariantContext[]{ctx0,ctx1};
		Genotype gts[]=new Genotype[]{g0,g1};
		out.writeStartElement("diff");
		out.writeAttribute("type", cat.name());
		out.writeAttribute("sample",sampleName);
		for(int i=0;i< 2;++i)
			{
			if(variants[i]==null) continue;
			out.writeStartElement("variant");
			out.writeAttribute("file",String.valueOf(i+1));
			out.writeAttribute("type",String.valueOf(variants[i].getType()));

			
			
			out.writeStartElement("chrom");
			out.writeCharacters(variants[i].getChr());
			out.writeEndElement();

			out.writeStartElement("pos");
			out.writeCharacters(String.valueOf(variants[i].getStart()));
			out.writeEndElement();

			
			
			out.writeStartElement("id");
			out.writeCharacters(variants[i].hasID()?variants[i].getID():".");
			out.writeEndElement();
				
			
			out.writeStartElement("ref");
			out.writeCharacters(String.valueOf(variants[i].getReference()));
			out.writeEndElement();
			
			out.writeStartElement("alts");
			out.writeCharacters(String.valueOf(variants[i].getAlternateAlleles()));
			out.writeEndElement();

			
			
			if(gts[i]!=null)
				{
				out.writeStartElement("genotype");
				out.writeAttribute("type",String.valueOf(gts[i].getType()));

				for(Allele a:gts[i].getAlleles())
					{
					out.writeStartElement("allele");
					out.writeCharacters(a.toString());
					out.writeEndElement();
					}
				if(gts[i].hasDP())
					{
					out.writeStartElement("dp");
					out.writeCharacters(String.valueOf(gts[i].getDP()));
					out.writeEndElement();
					}
				out.writeEndElement();
				}
			
			out.writeEndElement();
			}
		
	
		
		out.writeEndElement();
		out.writeCharacters("\n");
		}
	

	@SuppressWarnings("resource")
	@Override
	public int executeKnime(List<String> args)
		{
		PrintWriter exampleWriter=null;
		XMLStreamWriter  exampleOut=null;
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
			if(samples.isEmpty())
				{	
				error("No common samples");
				return -1;
				}
			
			Map<String, Counter<Category>> sample2info=new HashMap<String, Counter<Category>>(samples.size());
			for(String sampleName:samples)
				{
				sample2info.put(sampleName, new  Counter<Category>());
				}
			
			if(this.exampleFile!=null)
				{
				exampleWriter=new PrintWriter(exampleFile,"UTF-8");
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				exampleOut=xof.createXMLStreamWriter(exampleWriter);
				exampleOut.writeStartDocument("UTF-8", "1.0");
				exampleOut.writeStartElement("compare-callers");
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
					
					if(g0==null && g1==null)
						{
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.both_missing);
						continue;
						}
					else if(g0!=null && g1==null)
						{
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_1);
						
						if(ctx0.isIndel())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_1_indel);
							}
						else if(ctx0.isSNP())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_1_snp);
							}
						continue;
						}
					else if(g0==null && g1!=null)
						{
						
						if(ctx1.isIndel())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_2_indel);
							}
						else if(ctx1.isSNP())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_2_snp);
							}
						continue;
						}
					else
						{	
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context);
						if(ctx0.isIndel() && ctx1.isIndel())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_indel);
							}
						else if(ctx0.isSNP() && ctx1.isSNP())
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_snp);
							}
						
						if( (ctx0.hasID() && !ctx1.hasID()) ||
							(!ctx0.hasID() && ctx1.hasID()) ||
							(ctx0.hasID() && ctx1.hasID() && !ctx0.getID().equals(ctx1.getID()))
							)
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.common_context_discordant_id);
							}
						
						
						if(g0.sameGenotype(g1))
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same);

							if(g0.isHomRef())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same);
								}
							if(g0.isHomVar())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same_hom_var);
								}
							else if(g0.isHet())
								{	
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same_het);
								}
							}
						else
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant);

							if(g0.isHom() && g1.isHet())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_hom1_het2);
								}
							else if(g0.isHet() && g1.isHom())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_het1_hom2);
								}
							else if(g0.isHom() && g1.isHom())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_hom1_hom2);
								}
							else if(g0.isHet() && g1.isHet())
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_het1_het2);
								}
							else 
								{
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_but_discordant_others);
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
			
			if(exampleOut!=null)
				{
				exampleOut.writeEndElement();
				exampleOut.writeEndDocument();
				exampleOut.flush();
				exampleOut.close();
				}
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
			CloserUtil.close(exampleWriter);
			}
		
		}

	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:n:e:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile(opt.getOptArg()); break;
				case 'n': this.setNumberOfExampleVariants(Integer.parseInt(opt.getOptArg())); break;
				case 'e': this.setExampleFile(new File(opt.getOptArg())); break;
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
