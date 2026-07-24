package com.github.lindenb.jvarkit.variant.vcf;

import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

public class VCFWriterAnonymizer implements VariantContextWriter {
	final VariantContextWriter delegate;
	final UnaryOperator<String> transform= S->StringUtils.md5(S);
	VCFWriterAnonymizer(VariantContextWriter delegate) {
		this.delegate = delegate;
		}
	
	private String renameSample(String s) {
		return transform.apply(s);
		}
	
	private VCFHeader renameHeader(VCFHeader header) {
		if(!header.hasGenotypingData()) return header;
		return new VCFHeader(
				header.getMetaDataInInputOrder(),
				header.getSampleNamesInOrder().stream().map(S->renameSample(S)).collect(Collectors.toList())
				);
		}
	
	@Override
	public void writeHeader(VCFHeader header) {
		this.delegate.writeHeader(renameHeader(header));
	}

	@Override
	public void close() {
		this.delegate.close();		
	}

	@Override
	public boolean checkError() {
		return this.delegate.checkError();
	}

	@Override
	public void add(VariantContext vc) {
		if(!vc.hasGenotypes()) {
			this.delegate.add(vc);
			}
		else
			{
			this.delegate.add(
				new VariantContextBuilder(vc)
					.genotypes(vc.getGenotypes().stream().map(G->new GenotypeBuilder(G).name(renameSample(G.getSampleName())).make()).collect(Collectors.toList()))
					.make()
				);
			}
		}

	@Override
	public void setHeader(VCFHeader header) {
		this.delegate.setHeader(renameHeader(header));
	}

}
