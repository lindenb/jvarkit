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
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.Interval;

public class VcfCompareCallers
	extends AbstractVcfCompareCallers
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfCompareCallers.class);

	private enum Category
		{
		off_target_only_1,
		off_target_only_2,
		off_target_both,
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
	
	public VcfCompareCallers()
		{
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
		if(out==null || n> super.numberOfExampleVariants) return;
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
			out.writeCharacters(variants[i].getContig());
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
	
	private void internalTests()
		{
		Allele ref = Allele.create("A",true);
		List<Allele> alleles = new java.util.ArrayList<Allele>();
		alleles.add(ref);
		alleles.add(ref);
		Genotype g=new GenotypeBuilder("sample",alleles).make();
		if(g.isNoCall()) throw new IllegalStateException();
		if(!g.isCalled()) throw new IllegalStateException();
		if(!g.isAvailable()) throw new IllegalStateException();
		if(!g.isHomRef()) throw new IllegalStateException();
		}
	

	@Override
	public Collection<Throwable> call() throws Exception {
		htsjdk.samtools.util.IntervalTreeMap<Boolean> capture = null;
		PrintWriter exampleWriter=null;
		XMLStreamWriter  exampleOut=null;
		PrintStream pw=null;
		VcfIterator vcfInputs[]=new VcfIterator[]{null,null};
		VCFHeader headers[]=new VCFHeader[]{null,null};
		final List<String> args = getInputFiles();
		
		internalTests();
		
		try {
			if(args.size()==1)
				{
				LOG.info("Reading from stdin and "+ args.get(0));
				vcfInputs[0] = VCFUtils.createVcfIteratorStdin();
				vcfInputs[1] = VCFUtils.createVcfIterator( args.get(0));
				}
			else if(args.size()==2)
				{
				LOG.info("Reading from stdin and "+ args.get(0)+" and "+ args.get(1));
				vcfInputs[0] = VCFUtils.createVcfIterator( args.get(0));
				vcfInputs[1] = VCFUtils.createVcfIterator( args.get(1));
				}
			else
				{
				return wrapException(getMessageBundle("illegal.number.of.arguments"));
				}
		
			if( super.captureFile !=null )
				{
				LOG.info("Reading "+super.captureFile);
				capture = super.readBedFileAsBooleanIntervalTreeMap(super.captureFile);
				}
			
			for(int i=0;i< vcfInputs.length;++i)
				{
				headers[i] = vcfInputs[i].getHeader();
				}
			/* dicts */
			final SAMSequenceDictionary dict0 = headers[0].getSequenceDictionary();
			final SAMSequenceDictionary dict1 = headers[1].getSequenceDictionary();
			final Comparator<VariantContext> ctxComparator;
			if(dict0==null && dict1==null)
				{
				LOG.info("using createChromPosRefComparator");
				ctxComparator = VCFUtils.createChromPosRefComparator();
				}
			else if(dict0!=null && dict1!=null)
				{
				if( !SequenceUtil.areSequenceDictionariesEqual(dict0, dict1))
					{
					return wrapException(getMessageBundle("not.the.same.sequence.dictionaries"));
					}
				LOG.info("using createTidPosRefComparator");
				ctxComparator = VCFUtils.createTidPosRefComparator(dict0);
				}
			else
				{
				return wrapException(getMessageBundle("not.the.same.sequence.dictionaries"));
				}
			/* samples */
			Set<String> samples0=new HashSet<>(headers[0].getSampleNamesInOrder());
			Set<String> samples1=new HashSet<>(headers[1].getSampleNamesInOrder());
			Set<String> samples= new TreeSet<>(samples0);
			samples.retainAll(samples1);
			
			if(samples.size()!=samples0.size() || samples.size()!=samples1.size())
				{
				LOG.warn("Warning: Not the same samples set. Using intersection of both lists.");
				}
			if(samples.isEmpty())
				{	
				return wrapException("No common samples");
				}
			
			final Map<String, Counter<Category>> sample2info=new HashMap<String, Counter<Category>>(samples.size());
			for(final String sampleName:samples)
				{
				sample2info.put(sampleName, new  Counter<Category>());
				}
			
			if(super.exampleFile!=null)
				{
				exampleWriter=new PrintWriter(exampleFile,"UTF-8");
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				exampleOut=xof.createXMLStreamWriter(exampleWriter);
				exampleOut.writeStartDocument("UTF-8", "1.0");
				exampleOut.writeStartElement("compare-callers");
				}
			
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(dict0);
			final VariantContext buffer[]=new VariantContext[vcfInputs.length];
			final VariantContext prev[]=new VariantContext[vcfInputs.length];
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
								return wrapException("Input "+(i+1)+"/2 is not sorted"+(
									((i==0 && dict0==null) ||(i==1 && dict1==null))?
									"on chrom/pos/ref":
									"on sequence dictionary"
									)+". got\n"+buffer[i]+"\nafter\n"+prev[i]);
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
				Interval interval= null;

				if(buffer[0]!=null && ctxComparator.compare(buffer[0],smallest)==0)
					{
					prev[0] = progress.watch(vcfInputs[0].next());
					ctx0= prev[0];
					buffer[0]=null;
					interval = new Interval(ctx0.getContig(),ctx0.getStart(),ctx0.getEnd());
					}
				if(buffer[1]!=null && ctxComparator.compare(buffer[1],smallest)==0)
					{
					prev[1]= progress.watch(vcfInputs[1].next());
					ctx1= prev[1];
					buffer[1]=null;
					interval = new Interval(ctx1.getContig(),ctx1.getStart(),ctx1.getEnd());
					}
				
				//paranoid check
				if( buffer[0]!=null && buffer[1]!=null)
					{
					if( ctxComparator.compare(buffer[0],buffer[1])!=0) {
						exampleWriter.close();
						throw new IllegalStateException();
						}	
					}
				
				boolean in_capture=true;
				if(capture!=null && interval!=null)
					{
					in_capture = capture.containsOverlapping(interval);
					}
				
				for(final String sampleName: sample2info.keySet())
					{
					final Counter<Category> sampleInfo=sample2info.get(sampleName);
					
					
					Genotype g0=(ctx0==null?null:ctx0.getGenotype(sampleName));
					Genotype g1=(ctx1==null?null:ctx1.getGenotype(sampleName));
					
					
					
					
					if(g0!=null && (g0.isNoCall() || !g0.isAvailable())) g0=null;
					if(g1!=null && (g1.isNoCall() || !g1.isAvailable())) g1=null;
					
					/* when merging multiple vcf with GATK. there is no homref
					 * because everything is no call. Problem when comparing 
					 * with multiple called VCF.
					 */
					if(g0!=null && super.homRefIsNoCall && g0.isHomRef()) g0=null;
					if(g1!=null && super.homRefIsNoCall && g1.isHomRef()) g1=null;
					

					if(g0==null && g1==null)
						{
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.both_missing);
						continue;
						}
					else if(g0!=null && g1==null)
						{
						if(!in_capture)  {
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.off_target_only_1);
							continue;
							}
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
						if(!in_capture) {
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.off_target_only_2);
							continue;
							}
						watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.unique_to_file_2);
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
						if(!in_capture)
							{
							watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.off_target_both);
							continue;
							}
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
								watch(exampleOut,ctx0,ctx1,g0,g1,sampleName,sampleInfo,Category.called_and_same_hom_ref);
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
		
			pw  = openFileOrStdoutAsPrintStream();
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
			return RETURN_OK;
			} 
		catch (Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			if(getOutputFile()!=null)  CloserUtil.close(pw);
			CloserUtil.close(exampleWriter);
			}
		
		}

	
	
	public static void main(String[] args) {
		new VcfCompareCallers().instanceMainWithExit(args);
	}
}
