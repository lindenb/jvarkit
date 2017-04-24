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
package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;


/**
BEGIN_DOC

## See also

https://academic.oup.com/bioinformatics/article-abstract/33/7/964/2623048/Improved-VCF-normalization-for-accurate-VCF?redirectedFrom=fulltext

END_DOC
 */
@Program(name="vcffixindels",description="Fix samtools indels (for @SolenaLS)")
public class VCFFixIndels extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFFixIndels.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	public VCFFixIndels()
		{
		}
	
	
	@Override
	protected int doVcfToVcf(String inputName,
			VcfIterator r, VariantContextWriter w)
		{
		long nChanged=0L;
		final String TAG="INDELFIXED";
		final VCFHeader header=r.getHeader();
		
		final VCFHeader h2=new VCFHeader(header);
		addMetaData(h2);
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,1,VCFHeaderLineType.String,"Fix Indels for @SolenaLS (position|alleles...)"));

		w.writeHeader(h2);
		SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
		
		while(r.hasNext())
			{
			boolean somethingChanged=false;
			VariantContext ctx=progress.watch(r.next());
			/* map old allele to new allele */
			Map<Allele, Allele> original2modified=new HashMap<Allele, Allele>();
			/* did we found a strange allele (symbolic, etc ) */
			boolean strange_allele_flag=false;
			for(Allele a:ctx.getAlleles())
				{
				original2modified.put(a, a);
				if(a.isSymbolic() || a.isNoCall())
					{
					strange_allele_flag=true;
					break;
					}
				}
			
			if(strange_allele_flag ||
				original2modified.size()<2 /* at least 2 alleles: REF+ ALT */
				)
				{
				w.add(ctx);
				continue;
				}
			
			/* record chromStart if we need to shift to the right */
			int chromStart= ctx.getStart();
			
			/* trim right then left  */
			for(int side=0;side<2;++side)
				{
				boolean done=false;
				while(!done)
					{
					boolean ok_side=true;
					/* common nucleotide at end/start */
					Character targetChar = null;
					done=true;
					//scan side
					Set<Allele> keys = new HashSet<>(original2modified.keySet());  
					for(Allele k:keys)
						{
						Allele newAllele = original2modified.get(k);
						if(newAllele.isSymbolic()) 
							{
							ok_side = false;
							break;
							}
						String baseString = newAllele.getBaseString().trim().toUpperCase();
						if(baseString.length()<2)
							{
							ok_side = false;
							break;
							}
						/* first or last char or all sequences
						 * side==0 : right
						 * side==1 : left
						 * */
						Character baseChar=
							(side==0?
							baseString.charAt(baseString.length()-1):
							baseString.charAt(0)
							);
						if(targetChar==null)
							{
							targetChar = baseChar;
							}
						else if(!targetChar.equals(baseChar))
							{
							/* doesn't end with same nucleotide */
							ok_side = false;
							break;
							}
						}
					/* ok we can shift all alleles */
					if(ok_side && targetChar!=null)
						{
						done=false;
						somethingChanged=true;
						for(Allele k:keys)
							{
							Allele newAllele = original2modified.get(k);
							String baseString = newAllele.getBaseString().trim().toUpperCase();
							if(side==0)//remove last nucleotide
								{
								newAllele = Allele.create(
									baseString.substring(0,baseString.length()-1),
									newAllele.isReference());
								}
							else
								{
								newAllele = Allele.create(
										baseString.substring(1),
										newAllele.isReference());
								
								}
							original2modified.put(k, newAllele);
							}
						if(side==1) chromStart++;
						}
					}/* end of while done */
				
				}/* end side */
			
			if(!somethingChanged)
				{
				w.add(ctx);
				continue;
				}
			
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			b.start(chromStart);
			Allele newRef=original2modified.get(ctx.getReference());
			b.stop(chromStart+newRef.getBaseString().length()-1);
			b.alleles(original2modified.values());
			List<Genotype> genotypes=new ArrayList<>();
			for(String sample:header.getSampleNamesInOrder())
				{
				Genotype g = ctx.getGenotype(sample);
				if(g.isNoCall()) 
					{
					genotypes.add(g);
					continue;
					}
				GenotypeBuilder gb=new GenotypeBuilder(g);
				List<Allele> aL=new ArrayList<>();
				for(Allele a:g.getAlleles())
					{
					aL.add(original2modified.get(a));
					}
				gb.alleles(aL);
				genotypes.add(gb.make());
				}
			
			StringBuilder tagContent=new StringBuilder();
			tagContent.append(String.valueOf(ctx.getStart()));
			for(Allele a:ctx.getAlleles())
				{
				tagContent.append("|");
				tagContent.append(a.toString());
				}
			
			
			b.attribute(TAG,tagContent.toString());
			b.genotypes(genotypes);
			
			w.add( b.make());
			++nChanged;
			if(w.checkError()) break;
			}
		progress.finish();
		LOG.info("indels changed:"+nChanged);
		return RETURN_OK;
		}

	
	@Override
	public int doWork(List<String> args) {
		try {
			return doVcfToVcf(args,outputFile);
			} 
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		}
	
	
	public static void main(String[] args) throws IOException
		{
		new VCFFixIndels().instanceMainWithExit(args);
		}
}
