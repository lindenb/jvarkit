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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.Pedigree;


public class VCFTrios
	extends AbstractVCFTrios
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VCFTrios.class);

	private Pedigree pedigree=null;
	
    public VCFTrios()
    	{
    	}
	
    private static String allelesToString(Genotype g)
    	{
    	if(!g.isCalled()) return g.getSampleName()+" not called";
    	if(!g.isAvailable()) return g.getSampleName()+" not available";
    	String s=g.getSampleName()+":";
    	for(Allele a:g.getAlleles())
    		{
    		s+= " "+a.getDisplayString();
    		}
    	return s;
    	}
    
	private boolean isChilOf(
			Allele child1,Allele child2,
			Allele parent1,Allele parent2)
		{
		return	(child1.equals(parent1) && child2.equals(parent2)) ||
				(child1.equals(parent2) && child2.equals(parent1))
				;
		}

	
	private boolean trio(
			Allele child1,Allele child2,
			Allele father1,Allele father2,
			Allele mother1,Allele mother2
			)
		{		
		return	isChilOf(child1,child2,father1,mother1) ||
				isChilOf(child1,child2,father2,mother1) ||
				isChilOf(child1,child2,father1,mother2) ||
				isChilOf(child1,child2,father2,mother2)
				;
		}
	
	private boolean trio(Genotype child,Genotype father,Genotype mother)
		{
		return	trio(
				child.getAllele(0),child.getAllele(1),
				father.getAllele(0),father.getAllele(1),
				mother.getAllele(0),mother.getAllele(1)
				);
		}
	
	private boolean duo(
			Allele child1,Allele child2,
			Allele parent1,Allele parent2
			)
		{
		return	 child1.equals(parent1) ||
				 child1.equals(parent2) ||
				 child2.equals(parent1) ||
				 child2.equals(parent2)
				 ;
		}
	
	private boolean duo(Genotype child,Genotype parent)
		{
		return	duo(
				child.getAllele(0),child.getAllele(1),
				parent.getAllele(0),parent.getAllele(1)
				);
		}
	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName,
			VcfIterator r, VariantContextWriter w) throws IOException
		{
		int count_incompats=0;
		VCFHeader header=r.getHeader();
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				"MENDEL",
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"mendelian incompatibilities"
				));
		
		addMetaData(h2);

		if( create_filter) h2.addMetaDataLine(new VCFFilterHeaderLine("MENDEL", "data filtered with VCFTrios"));
		w.writeHeader(h2);
		
		Map<String,Pedigree.Person> samplename2person=new HashMap<String,Pedigree.Person>(h2.getSampleNamesInOrder().size());
		for(String sampleName:h2.getSampleNamesInOrder())
			{
			Pedigree.Person p=null;
			for(Pedigree.Family f:this.pedigree.getFamilies())
				{
				for(Pedigree.Person child:f.getIndividuals())
					{
					if(child.getId().equals(sampleName))
						{
						if(p!=null)
							{
							throw new IllegalArgumentException(sampleName+" found twice in pedigree !");
							}
						p=child;
						}
					}
				}
			if(p==null)
				{
				LOG.info("Cannot find "+sampleName+" in "+pedigreeFile);
				}
			else
				{
				samplename2person.put(sampleName, p);
				}
			}
		
		LOG.info("persons in pedigree: "+samplename2person.size());
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		while(r.hasNext())
			{
			VariantContext ctx= progress.watch(r.next());
			
			Set<String> incompatibilities=new HashSet<String>();
			
			
			for(Pedigree.Person child:samplename2person.values())
				{
				Genotype gChild=ctx.getGenotype(child.getId());
				if(gChild==null)
					{
					LOG.debug("cannot get genotype for child  "+child.getId());
					continue;
					}
				if(gChild.isNoCall())
					{
					continue;
					}
				if(!gChild.isAvailable())
					{
					continue;
					}
				if(gChild.getAlleles().size()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles child:"+ allelesToString(gChild));
					continue;
					}
				
				Pedigree.Person parent=child.getFather();
				Genotype gFather=(parent==null?null:ctx.getGenotype(parent.getId()));
				if(gFather==null && parent!=null)
					{
					LOG.warn("cannot get genotype for father  "+parent.getId());
					}
				if(gFather!=null && gFather.isNoCall()) gFather=null;
				if(gFather!=null && !gFather.isAvailable()) gFather=null;

				if(gFather!=null && gFather.getAlleles().size()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles father: "+ allelesToString(gFather));
					gFather=null;
					}
				parent=child.getMother();
				
				Genotype gMother=(parent==null?null:ctx.getGenotype(parent.getId()));
				
				if(gMother==null && parent!=null)
					{
					LOG.debug("cannot get genotype for mother  "+parent.getId());
					}
				
				if(gMother!=null && gMother.isNoCall()) gMother=null;
				if(gMother!=null && !gMother.isAvailable()) gMother=null;
				if(gMother!=null && gMother.getAlleles().size()!=2)
					{
					LOG.debug(getClass().getSimpleName()+" only handle two alleles mother:"+ allelesToString(gMother));
					gMother=null;
					}
				
				boolean is_ok=true;
				if(gFather!=null && gMother!=null)
					{
					is_ok=trio(gChild,gFather,gMother);
					}
				else if(gFather!=null)
					{
					is_ok=duo(gChild,gFather);
					}
				else if(gMother!=null)
					{
					is_ok=duo(gChild,gMother);
					}
				if(!is_ok)
					{
					incompatibilities.add(child.getId());
					}
				}
				
			if(incompatibilities.isEmpty())
				{
				w.add(ctx);
				continue;
				}
			++count_incompats;
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			if( create_filter) b.filter("MENDEL");
			b.attribute("MENDEL", incompatibilities.toArray());
			w.add(b.make());
			
			if(w.checkError()) break;
			}	
		progress.finish();
		LOG.info("incompatibilities N="+count_incompats);
		return RETURN_OK;
		}
	
	

	
	@Override
	public void disposeKnime() {
		this.pedigree=null;
	}

	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
	}
	
	
	public static void main(String[] args)
		{
		new VCFTrios().instanceMainWithExit(args);
		}

	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.pedigreeFile==null)
			{
			return wrapException("Pedigree undefined.");
			}
		BufferedReader in = null;
		try {
			LOG.info("reading pedigree "+this.pedigreeFile);
			 in=IOUtils.openFileForBufferedReading(this.pedigreeFile);
			this.pedigree=Pedigree.readPedigree(in);
			in.close();
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		return super.initializeKnime();
		}



	}
