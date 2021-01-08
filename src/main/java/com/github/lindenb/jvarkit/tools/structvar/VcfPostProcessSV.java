/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.structvar;

import java.util.Arrays;
import java.util.Comparator;
import java.util.ArrayList;
import java.util.List;
import java.util.Collections;
import java.util.stream.Collectors;

import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.Breakend;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFRecordCodec;

/**
BEGIN_DOC

## input

Input is the vcf output of manta/lumpy etc...

## Example:

todo

```
```

END_DOC

 */
@Program(name="vcfpostprocesssv",
description="Postprocess BND pairs in a MANTA/LUMPY-SV.",
keywords= {"cnv","sv","gidss","vcf"},
creationDate="20200612",
modificationDate="20200612"
)
public class VcfPostProcessSV extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfPostProcessSV.class).make();
	@Parameter(names={"-k","--key","--keys"},description="Comma separated list of INFO field to search the ID of mate breakend. Manta: MATEID, lumpy: EVENT ")
	private String keys="MATEID,EVENT";
	@Parameter(names={"-A","--allele"},description="New alternate allele symbol '<x>' after pair of BND junction have been found.")
	private String newAlternateSymbol = "OTHER";
	@Parameter(names={"-no-summary"},description="remove summary count")
	private boolean disable_summary = false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {	
		SortingCollection<VariantContext> sorter1 = null;
		SortingCollection<VariantContext> sorter2 = null;
		try {		
			final Counter<String> counter = new Counter<>();
			if(StringUtils.isBlank(this.keys)) {
				LOG.error("empty --keys");
				return -1;
				}
			if(StringUtils.isBlank(this.newAlternateSymbol)) {
				LOG.error("empty --allele");
				return -1;
				}

			final VCFHeader header= iterin.getHeader();
			
			if(header.getInfoHeaderLine(VCFConstants.SVTYPE)==null) {
				header.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1,VCFHeaderLineType.String,"Type of structural variant"));
				}
			if(header.getInfoHeaderLine("SVLEN")==null) {
				header.addMetaDataLine(new VCFInfoHeaderLine("SVLEN", 1,VCFHeaderLineType.Integer,"Difference in length between REF and ALT alleles"));
				}
			if(header.getInfoHeaderLine(VCFConstants.END_KEY)==null) {
				header.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,VCFHeaderLineType.Integer,"End position of the variant described in this record"));
				}
			
			
			final String mateID = Arrays.stream(CharSplitter.COMMA.split(this.keys)).
					map(S->S.trim()).
					filter(S->!StringUtils.isBlank(S)).
					map(S->header.getInfoHeaderLine(S)).
					filter(H->H!=null && H.getType().equals(VCFHeaderLineType.String) ).
					map(H->H.getID()).findFirst().
					orElse(null);
			if(StringUtils.isBlank(mateID)) {
				LOG.error("Cannot find INFO for mate-id (was "+this.keys+")");
				return -1;
				}
			LOG.info("Mate key is '" + mateID +"'. Reading input.");
						
			final Comparator<VariantContext> comparator = (A,B)->   A.getAttributeAsString(mateID, "").compareTo(B.getAttributeAsString(mateID, ""));
			
			sorter1 = SortingCollection.newInstance(
					VariantContext.class,
                    new VCFRecordCodec(header,true),
                    comparator,
                    this.writingSortingCollection.getMaxRecordsInRam(),
                    this.writingSortingCollection.getTmpPaths()
                    );
			
			sorter2 = SortingCollection.newInstance(
					VariantContext.class,
                    new VCFRecordCodec(header,true),
                    header.getVCFRecordComparator(),
                    this.writingSortingCollection.getMaxRecordsInRam(),
                    this.writingSortingCollection.getTmpPaths()
                    );
			
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				if(ctx.hasAttribute(mateID) && (!ctx.hasAttribute(VCFConstants.SVTYPE) || ctx.getAttributeAsString(VCFConstants.SVTYPE,"").equals("BND"))) {
					sorter1.add(ctx);
					}
				else
					{
					counter.incr("Not a BND variant");
					sorter2.add(ctx);
					}
				}
			sorter1.doneAdding();
			sorter1.setDestructiveIteration(true);
			
			CloseableIterator<VariantContext> iter2=sorter1.iterator();
			PeekIterator<VariantContext> peek = new PeekIterator<>(iter2);
						
			while(peek.hasNext()) {
				final VariantContext first = peek.next();
				final List<VariantContext> variants = new ArrayList<>();
				variants.add(first);
				while(peek.hasNext()) {
					final VariantContext second = peek.peek();
					if(first.hasID() && first.getID().equals(second.getAttributeAsString(mateID,""))) {
						variants.add(peek.next());
					} else if(second.hasID() && second.getID().equals(first.getAttributeAsString(mateID,""))) {
						variants.add(peek.next());
					} else if(first.getAttributeAsString(mateID,"").equals(second.getAttributeAsString(mateID,""))) {
						variants.add(peek.next());
					} else {
					break;
					}
				}


				if(variants.size()!=2) {
					counter.incr("Not a pair of Mate ("+variants.size()+")", variants.size());
					for(final VariantContext ctx:variants) {
						sorter2.add(ctx);
						}
					continue;
					}
				
				Collections.sort(variants,(A,B)->Integer.compare(A.getStart(),B.getStart()));

				final VariantContext ctx1 = variants.get(0);
				final VariantContext ctx2 = variants.get(1);
				if(!(ctx1.getNAlleles()==2 && ctx1.getNAlleles()==2)) {
					counter.incr("expected two alleles.",2L);
					sorter2.add(ctx1);
					sorter2.add(ctx2);
					continue;
					}
					
				final Breakend brk1 = Breakend.parse(ctx1.getAlleles().get(1)).orElse(null);
				final Breakend brk2 = Breakend.parse(ctx2.getAlleles().get(1)).orElse(null);
				if(brk1==null || brk2==null) {
					counter.incr("ALT is not breakend.",2L);
					sorter2.add(ctx1);
					sorter2.add(ctx2);
					}
				/* should be the same breakend 
				difference can be large use CIPOS / CIEND ?
				if( !ctx1.getContig().equals(brk2.getContig()) ||
					!ctx2.getContig().equals(brk1.getContig()) ||
					Math.abs(ctx1.getStart()-brk2.getStart())>1  ||
					Math.abs(ctx2.getStart()-brk1.getStart())>1) {
					counter.incr("mate is not the same position.",2L);
					sorter2.add(ctx1);
					sorter2.add(ctx2);
					continue;
					}
				*/
				final VariantContextBuilder vcb1 = new VariantContextBuilder(ctx1);						
				final VariantContextBuilder vcb2 = new VariantContextBuilder(ctx2);

				if(!ctx1.contigsMatch(ctx2)) {
					vcb1.attribute(VCFConstants.SVTYPE, "TRANSLOC");
					vcb2.attribute(VCFConstants.SVTYPE, "TRANSLOC");
					sorter2.add(vcb1.make());
					sorter2.add(vcb2.make());
					counter.incr("translocation.",2L);
					continue;
					}
					
						
				final Allele ctx1_alt = ctx1.getAlleles().get(1);
				final List<Allele> alleles = Arrays.asList(ctx1.getReference(),Allele.create("<"+this.newAlternateSymbol.trim()+">",false));
				vcb1.attribute(VCFConstants.SVTYPE,this.newAlternateSymbol.trim());
				vcb1.alleles(alleles);
				final int ctx_end = Math.max(ctx1.getEnd(),ctx2.getEnd());
				vcb1.stop(ctx_end);
				vcb1.attribute("END",ctx_end);
				vcb1.attribute("SVLEN",CoordMath.getLength(ctx1.getStart(),ctx_end));
				vcb1.rmAttribute(mateID);
				vcb1.genotypes(ctx1.getGenotypes().stream().map(GT->new GenotypeBuilder(GT).
						alleles(GT.getAlleles().
							stream().
							map(A->A.equals(ctx1_alt)?alleles.get(1):A).
							collect(Collectors.toList())).
							make()
					).collect(Collectors.toList()));
				sorter2.add(vcb1.make());
				counter.incr("Two BND variants converted.",2L);
				}
			iter2.close();
			sorter1.cleanup();
			sorter1= null;
			
			sorter2.doneAdding();
			sorter2.setDestructiveIteration(true);

			iter2=sorter2.iterator();

			
			final VCFHeader header2 = new VCFHeader(header);				
			JVarkitVersion.getInstance().addMetaData(this, header2);			
			out.writeHeader(header2);
			while(iter2.hasNext()) {
				out.add(iter2.next());
				}
		    iter2.close();
		    sorter2.cleanup();
		    sorter2=null;
		    if(!disable_summary) {
		    	LOG.info("SUMMARY COUNT");
			    for(final String key:counter.keySet()) {
			    	LOG.info("\t"+key+" : "+ counter.count(key));
			    }
		    }
		return 0;
	} catch(final Throwable err) {
		err.printStackTrace();
		LOG.error(err);
		return -1;
	} finally {
		if(sorter1!=null) sorter1.cleanup();
		if(sorter2!=null) sorter2.cleanup();
	}
}
	
	public static void main(final String[] args) {
		new VcfPostProcessSV().instanceMainWithExit(args);
	}

}
