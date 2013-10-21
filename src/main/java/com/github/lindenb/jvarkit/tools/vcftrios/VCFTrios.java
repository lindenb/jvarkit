package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFFilterHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;

public class VCFTrios extends AbstractVCFFilter
	{
	private static final Log LOG = Log.getInstance(VCFTrios.class);
    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + "Find mendelian incompatibilitie in VCF (plink format). ";
    @Option(shortName="PED", doc=" Pedigree file",optional=false)
    public File PEDIGREE=null;
    @Option(shortName="PF", doc="Set filter 'MENDEL' if incompatibilities found. ",optional=false)
    public boolean FILTER=false;

	
	
	
	
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
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		LOG.info("reading pedigree");
		Pedigree pedigree=Pedigree.readPedigree(PEDIGREE);
		int count_incompats=0;
		
		VCFHeader header=r.getHeader();
		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine("MENDEL", VCFHeaderLineCount.INTEGER, VCFHeaderLineType.String, "mendelian incompatibilities. PEd file: "+PEDIGREE));
		if( FILTER) h2.addMetaDataLine(new VCFFilterHeaderLine("MENDEL", "data filtered with "+this.getCommandLine()));
		w.writeHeader(h2);
		
		Map<String,Pedigree.Person> samplename2person=new HashMap<String,Pedigree.Person>(h2.getSampleNamesInOrder().size());
		for(String sampleName:h2.getSampleNamesInOrder())
			{
			Pedigree.Person p=null;
			for(Pedigree.Family f:pedigree.getFamilies())
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
				LOG.info("Cannot find "+sampleName+" in "+PEDIGREE);
				}
			else
				{
				samplename2person.put(sampleName, p);
				}
			}
		
		LOG.info("persons in pedigree: "+samplename2person.size());
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			
			
			Set<String> incompatibilities=new HashSet<String>();
			
			
			for(Pedigree.Person child:samplename2person.values())
				{
				Genotype gChild=ctx.getGenotype(child.getId());
				if(gChild==null)
					{
					LOG.debug("cannot get genotype for child  "+child.getId());
					continue;
					}
				if(!gChild.isCalled())
					{
					continue;
					}
				if(gChild.getAlleles().size()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles");
					continue;
					}
				
				Pedigree.Person parent=child.getFather();
				Genotype gFather=(parent==null?null:ctx.getGenotype(parent.getId()));
				if(gFather==null && parent!=null)
					{
					LOG.debug("cannot get genotype for father  "+parent.getId());
					}
				if(gFather!=null && gFather.getAlleles().size()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles");
					gFather=null;
					}
				if(gFather!=null && !gFather.isCalled()) gFather=null;
				parent=child.getMother();
				
				Genotype gMother=(parent==null?null:ctx.getGenotype(parent.getId()));
				
				if(gMother==null && parent!=null)
					{
					LOG.debug("cannot get genotype for mother  "+parent.getId());
					}
				
				if(gMother!=null && !gMother.isCalled()) gMother=null;
				if(gMother!=null && gMother.getAlleles().size()!=2)
					{
					LOG.warn(getClass().getSimpleName()+" only handle two alleles");
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
			if( FILTER) b.filter("MENDEL");
			b.attribute("MENDEL", incompatibilities.toArray());
			w.add(b.make());
			}		
		LOG.info("incompatibilities N="+count_incompats);
		}

	public static void main(String[] args)
		{
		new VCFTrios().instanceMainWithExit(args);
		}

	}
