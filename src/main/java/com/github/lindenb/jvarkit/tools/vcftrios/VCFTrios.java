package com.github.lindenb.jvarkit.tools.vcftrios;

import java.io.IOException;
import java.util.Iterator;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.AbstractVCFFilter;

public class VCFTrios extends AbstractVCFFilter
	{
	public String FATHER;
	public String MOTHER;
	public String CHILD;
	
	private boolean childIs(
			Allele child1,Allele child2,
			Allele parent1,Allele parent2)
		{
		return	child1.equals(parent1) && child2.equals(parent2) ||
				child1.equals(parent2) && child2.equals(parent1)
				;
		}
	
	private boolean trio(
			Allele father1,Allele father2,
			Allele mother1,Allele mother2,
			Allele child1,Allele child2
			)
		{
		return	childIs(child1,child2,father1,mother1) ||
				childIs(child1,child2,father2,mother1) ||
				childIs(child1,child2,father1,mother2) ||
				childIs(child1,child2,father2,mother2)
				;
		}
	
	@Override
	protected void doWork(LineReader r, VariantContextWriter w)
			throws IOException
		{
		VCFCodec codeIn=new VCFCodec();		
		VCFHeader header=(VCFHeader)codeIn.readHeader(r);
		VCFHeader h2=new VCFHeader(header);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName(),"TODO"));
		w.writeHeader(h2);
		String line;
		
		while((line=r.readLine())!=null)
			{
			VariantContext ctx=codeIn.decode(line);
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			Genotype gFather=ctx.getGenotype(FATHER);
			Genotype gMother=ctx.getGenotype(MOTHER);
			Genotype gChild=ctx.getGenotype(CHILD);
			trio(
				gFather.getAllele(0),gFather.getAllele(1),
				gMother.getAllele(0),gMother.getAllele(1),
				gChild.getAllele(0),gChild.getAllele(1)
				);
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
