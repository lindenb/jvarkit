package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.AbstractVCFFilter;

public class VCFTrios extends AbstractVCFFilter
	{
	private Map<String,Map<String,Person> > pedigree;
	
	private class Person
		{
		String family;
		String id;
		String fatherId=null;
		String motherId=null;
		
		public Map<String,Person> getFamily()
			{
			return VCFTrios.this.pedigree.get(family);
			}
		
		public String getId()
			{
			return this.id;
			}
		
		private Person getParent(String s)
			{
			if(s==null || s.isEmpty() || s.equals("0")) return null;
			Map<String,Person> fam=getFamily();
			if(fam==null) return null;
			return fam.get(s);
			}
		
		public Person getFather()
			{
			return getParent(fatherId);
			}
		
		public Person getMother()
			{
			return getParent(motherId);
			}
		}
	
	private static boolean genotypeHasAllele(final Genotype g,Allele a)
		{
		for(Allele A:g.getAlleles())
			{
			if(A.equals(a)) return true;
			}
		return false;
		}
	
	private boolean childIs(
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
		return	childIs(child1,child2,father1,mother1) ||
				childIs(child1,child2,father2,mother1) ||
				childIs(child1,child2,father1,mother2) ||
				childIs(child1,child2,father2,mother2)
				;
		}
	
	private boolean trio(Genotype child,Genotype father,Genotype mother)
		{
		Set<Allele> childL=new HashSet<Allele>(child.getAlleles());
		Set<Allele> motherL=new HashSet<Allele>(mother.getAlleles());
		Set<Allele> fatherL=new HashSet<Allele>(father.getAlleles());
		
		Set<Allele> fromHisFather=new HashSet<Allele>(childL);
		fromHisFather.retainAll(fatherL);
		if(fromHisFather.isEmpty()) return false;
		
		Set<Allele> fromHisMother=new HashSet<Allele>(childL);
		fromHisMother.retainAll(motherL);
		if(fromHisMother.isEmpty()) return false;
		
		
		
		
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
	protected void doWork(LineReader r, VariantContextWriter w)
			throws IOException
		{
		VCFCodec codeIn=new VCFCodec();		
		VCFHeader header=(VCFHeader)codeIn.readHeader(r);
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"TODO"));
		w.writeHeader(h2);
		String line;
		
		while((line=r.readLine())!=null)
			{
			VariantContext ctx=codeIn.decode(line);
			Set<String> incompatibilities=new HashSet<String>();
			for(Map<String,Person> family:pedigree.values())
				{
				for(Person child:family.values())
					{
					Genotype gChild=ctx.getGenotype(child.getId());
					if(gChild==null) continue;
					Person parent=child.getFather();
					Genotype gFather=(parent==null?null:ctx.getGenotype(parent.getId()));
					if(gFather!=null && !gFather.isCalled()) gFather=null;
					parent=child.getMother();
					Genotype gMother=(parent==null?null:ctx.getGenotype(parent.getId()));
					if(gMother!=null && !gMother.isCalled()) gMother=null;
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
				}
			if(incompatibilities.isEmpty())
				{
				w.add(ctx);
				continue;
				}
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			b.filter("MENDEL");
			b.attribute("MENDEL", incompatibilities.toArray());
			}		
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFTrios().instanceMainWithExit(args);
		}

	}
