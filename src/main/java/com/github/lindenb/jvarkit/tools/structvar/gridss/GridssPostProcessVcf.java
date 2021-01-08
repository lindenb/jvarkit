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
package com.github.lindenb.jvarkit.tools.structvar.gridss;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.Breakend;

import htsjdk.samtools.util.CloseableIterator;
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
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## input

Input is the vcf output of Gridss

## Example:

```
$ wget -O - -q "https://raw.githubusercontent.com/genome/breakdancer/master/test-data/expected_output" | java -jar dist/breakdancer2vcf.jar -R src/test/resources/human_b37.dict 

##fileformat=VCFv4.2
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimate allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CHROM2,Number=1,Type=String,Description="Chromosome 2">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ORIENT1,Number=1,Type=String,Description="Orientation 1">
##INFO=<ID=ORIENT2,Number=1,Type=String,Description="Orientation 2">
##INFO=<ID=POS2,Number=1,Type=String,Description="Position 2">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variation type">
##breakdancer.command=bdfast -o21 inv_del_bam_config
##breakdancer.version=1.4.1-unstable-10-fdfe9f2-dirty (commit fdfe9f2-dirty)
##breakdancer2vcf.meta=compilation:20200511120615 githash:c830e0b htsjdk:2.21.3 date:20200511120941 cmd:-R src/test/resources/human_b37.dict
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
(...)
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_IJ-NA19238-NA19238-extlibs	H_IJ-NA19240-NA19240-extlibs
21	29185056	.	N	<INS>	99	.	DP=2;END=29185377;SVLEN=226;SVTYPE=INS	GT:DP:GQ	0/1:1:99	0/1:1:99
21	29185462	.	N	<DEL>	99	.	DP=21;END=29186122;SVLEN=-545;SVTYPE=DEL	GT:AF:DP:GQ	0/1:176.58:21:99	0/1:167.89:.:99
21	34807694	.	N	<INS>	99	.	DP=3;END=34808852;SVLEN=304;SVTYPE=INS	GT:DP:GQ	0/1:1:99	0/1:2:99
21	34808937	.	N	<INV>	99	.	DP=2;END=34809799;SVLEN=-737;SVTYPE=INV	GT:AF:DP:GQ	0/1:847.39:.:99	0/1:878.83:2:99
```



END_DOC

 */
@Program(name="gridsspostprocessvcf",
description="Postprocess gridss (change SVTYPE, SVLEN.. etc...)",
keywords= {"cnv","sv","gidss","vcf"},
generate_doc=false,
creationDate="20200515",
modificationDate="20200518"
)
public class GridssPostProcessVcf extends OnePassVcfLauncher {
	private static final String EVENT_KEY="EVENT";
	private static final Logger LOG = Logger.build(GridssPostProcessVcf.class).make();
	@Parameter(names={"-D","--debug-file"},description="Debug File.",hidden=true)
	private Path debugFile=null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private enum BndType { DNA_OPEN,OPEN_DNA, DNA_CLOSE,CLOSE_DNA};
	
	private String toString(final VariantContext ctx) {
		final StringBuilder sb = new StringBuilder(ctx.getContig()).
				append(":").
				append(ctx.getStart()).
				append("-").
				append(ctx.getEnd()).
				append(" ").
				append(ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))).
				append(" ").
				append(ctx.getAttributeAsString(EVENT_KEY, ""));
		return sb.toString();
	}
	
	private BndType getBndType(final Breakend brk) {
	if(brk.getDelimiter() =='[' ) {
		return brk.isRight() ? BndType.OPEN_DNA : BndType.DNA_OPEN ;
		}
	else // ']'
		{
		return brk.isRight() ? BndType.CLOSE_DNA : BndType.DNA_CLOSE ;
		}
	}
	
	/* cannot use AcidNucleics.isATGC because there can be a dot */
	private boolean isATGC(final Allele a) {
	// do not test isSymbolic because '.' at the beginning of DNA string makes it symbolic ??
	final String s = a.getDisplayString();
	for(int i=0;i< s.length();i++) {
		switch(s.charAt(i)) {
		case '.' : case 'A': case 'T' : case 'C': case 'G': case 'N': break;
		default: return false;
		}
	}
	return true;
	}
	
	/* symbolic allele like '.AATCGTACGAT' have length =0, so I cannot use Allele.length() */
	private int length(final Allele a) {
	int n=0;
	final String s = a.getDisplayString();
	for(int i=0;i< s.length();i++) {
		switch(s.charAt(i)) {
		case '.' : break;
		case 'A': case 'T' : case 'C': case 'G': case 'N': n++; break;
		default: return 0;
		}
	}
	return n;
	}
	
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {	
		SortingCollection<VariantContext> sorter1 = null;
		SortingCollection<VariantContext> sorter2 = null;
		PrintWriter debug = null;
		try {
			debug = this.debugFile == null?
				new PrintWriter(new NullOuputStream()):
				new PrintWriter(Files.newBufferedWriter(this.debugFile))
				;
		
		
			final VCFHeader header= iterin.getHeader();
			LOG.info("reading input.");
			final Comparator<VariantContext> comparator = (A,B)->  A.getAttributeAsString(EVENT_KEY, "").compareTo(B.getAttributeAsString(EVENT_KEY, ""));
			sorter1 = SortingCollection.newInstance(
					VariantContext.class,
                    new VCFRecordCodec(header,true),
                    comparator,
                    this.writingSortingCollection.getMaxRecordsInRam(),
                    this.writingSortingCollection.getTmpPaths()
                    );
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				if(!ctx.hasAttribute(EVENT_KEY)) {
					LOG.warn("skipping variant without INFO/"+EVENT_KEY+" at "+toString(ctx));
					continue;
					}
				sorter1.add(ctx);
				}
			sorter1.doneAdding();
			sorter1.setDestructiveIteration(true);
			LOG.info("done adding");
			
			CloseableIterator<VariantContext> iter2=sorter1.iterator();
			@SuppressWarnings("resource")
			final EqualRangeIterator<VariantContext> equal_range = new EqualRangeIterator<>(iter2, comparator);
			
			sorter2 = SortingCollection.newInstance(
					VariantContext.class,
                    new VCFRecordCodec(header,true),
                    header.getVCFRecordComparator(),
                    this.writingSortingCollection.getMaxRecordsInRam(),
                    this.writingSortingCollection.getTmpPaths()
                    );
			
			while(equal_range.hasNext()) {
				final List<VariantContext> variants = equal_range.next();
				if(variants.size()==1) {
					final VariantContext ctx = variants.get(0);
					final Allele alt = ctx.getAlleles().get(1);
					// INSERTION LIKE. chr22:1234 A/.AAAACAAGGAG
					if(	isATGC(ctx.getReference()) &&
						isATGC(alt) &&
						length(ctx.getReference())< length(alt) ) {
						final VariantContextBuilder vcb1 = new VariantContextBuilder(ctx);		
						vcb1.attribute(VCFConstants.SVTYPE, "INS");
						vcb1.attribute("SVLEN",length(alt) - length(ctx.getReference()) );
						sorter2.add(vcb1.make());
						}
					//STRANGE ? no ? change
					else if(	
						isATGC(ctx.getReference()) &&
						isATGC(alt) &&
						length(ctx.getReference()) ==1 &&
						length(alt)==1 ) {
						sorter2.add(ctx);
						}
					else
						{
						sorter2.add(ctx);
						debug.println("ALONE "+toString(ctx));
						}
					}
				else if(variants.size()!=2) {
					for(final VariantContext ctx:variants) {
						debug.println("SIZE>2 "+toString(ctx));
						sorter2.add(ctx);
						}
					}
				else  {
					final VariantContext ctx1 = variants.get(0);
					final VariantContext ctx2 = variants.get(1);
					final Breakend brk1 = Breakend.parse(ctx1.getAlleles().get(1)).orElse(null);
					final Breakend brk2 = Breakend.parse(ctx2.getAlleles().get(1)).orElse(null);
					if(brk1==null || brk2==null) {
						debug.println("expected two bnd "+ toString(ctx1)+" "+toString(ctx2));
						sorter2.add(ctx1);
						sorter2.add(ctx2);
						return -1;
						//continue;
						}
					/* should be the same breakend */
					if( !ctx1.getContig().equals(brk2.getContig()) ||
						!ctx2.getContig().equals(brk1.getContig()) ||
						ctx1.getStart()!=brk2.getStart()  ||
						ctx2.getStart()!=brk1.getStart()) {
						debug.println("expected concordant bnd "+ toString(ctx1)+" "+toString(ctx2));
						sorter2.add(ctx1);
						sorter2.add(ctx2);
						return -1;
						//continue;
						}
					final VariantContextBuilder vcb1 = new VariantContextBuilder(ctx1);						
					final VariantContextBuilder vcb2 = new VariantContextBuilder(ctx2);

					if(!ctx1.contigsMatch(ctx2)) {
						vcb1.attribute(VCFConstants.SVTYPE, "CTX");
						vcb2.attribute(VCFConstants.SVTYPE, "CTX");
						sorter2.add(vcb1.make());
						sorter2.add(vcb2.make());
						continue;
						}
					
					if(ctx1.getStart()>ctx2.getStart()) {
						debug.println("CTX1>CTX2?" +toString(ctx1)+" "+toString(ctx2));
						sorter2.add(vcb1.make());
						sorter2.add(vcb2.make());
						continue;
						}
						
					final BndType bndType1  = getBndType(brk1);
					final BndType bndType2  = getBndType(brk2);
						
					if(bndType1.equals(BndType.DNA_OPEN) && bndType2.equals(BndType.CLOSE_DNA)) {
						final Allele ctx1_alt = ctx1.getAlleles().get(1);
						final List<Allele> alleles = Arrays.asList(ctx1.getReference(),Allele.create("<DEL>",false));
						vcb1.attribute(VCFConstants.SVTYPE, "DEL");
						vcb1.alleles(alleles);
						vcb1.attribute("SVLEN",ctx1.getEnd() - ctx2.getStart()  +1 );
						vcb1.genotypes(ctx1.getGenotypes().stream().map(GT->new GenotypeBuilder(GT).
								alleles(GT.getAlleles().
									stream().
									map(A->A.equals(ctx1_alt)?alleles.get(1):A).
									collect(Collectors.toList())).
									make()
							).collect(Collectors.toList()));
						sorter2.add(vcb1.make());
						}
					else if(bndType1.equals(BndType.CLOSE_DNA) && bndType2.equals(BndType.DNA_OPEN)) {
						final Allele ctx1_alt = ctx1.getAlleles().get(1);
						final List<Allele> alleles = Arrays.asList(ctx1.getReference(),Allele.create("<DUP>",false));
						vcb1.attribute(VCFConstants.SVTYPE, "DUP");
						vcb1.alleles(alleles);
						vcb1.attribute("SVLEN",ctx2.getEnd() - ctx1.getStart()  +1 );
						vcb1.genotypes(ctx1.getGenotypes().stream().map(GT->new GenotypeBuilder(GT).
								alleles(GT.getAlleles().
									stream().
									map(A->A.equals(ctx1_alt)?alleles.get(1):A).
									collect(Collectors.toList())).
									make()
							).collect(Collectors.toList()));
						sorter2.add(vcb1.make());
						}
					else if(bndType1.equals(BndType.OPEN_DNA) && bndType2.equals(BndType.OPEN_DNA)) {
						final Allele ctx1_alt = ctx1.getAlleles().get(1);
						final List<Allele> alleles = Arrays.asList(ctx1.getReference(),Allele.create("<INV>",false));
						vcb1.attribute(VCFConstants.SVTYPE, "INV");
						vcb1.alleles(alleles);
						vcb1.attribute("SVLEN",ctx2.getEnd() - ctx1.getStart()  +1 );
						vcb1.genotypes(ctx1.getGenotypes().stream().map(GT->new GenotypeBuilder(GT).
								alleles(GT.getAlleles().
									stream().
									map(A->A.equals(ctx1_alt)?alleles.get(1):A).
									collect(Collectors.toList())).
									make()
							).collect(Collectors.toList()));
						sorter2.add(vcb1.make());
						}
					else if(bndType1.equals(BndType.DNA_CLOSE) && bndType2.equals(BndType.DNA_CLOSE)) {
						final Allele ctx1_alt = ctx1.getAlleles().get(1);
						final List<Allele> alleles = Arrays.asList(ctx1.getReference(),Allele.create("<SV>",false));
						vcb1.attribute(VCFConstants.SVTYPE, "SV");
						vcb1.alleles(alleles);
						vcb1.attribute("SVLEN",ctx2.getEnd() - ctx1.getStart()  +1 );
						vcb1.genotypes(ctx1.getGenotypes().stream().map(GT->new GenotypeBuilder(GT).
								alleles(GT.getAlleles().
									stream().
									map(A->A.equals(ctx1_alt)?alleles.get(1):A).
									collect(Collectors.toList())).
									make()
							).collect(Collectors.toList()));
						sorter2.add(vcb1.make());
						}
					else
						{
						debug.println("How to handle "+toString(ctx1)+" "+toString(ctx2));
						sorter2.add(ctx1);
						sorter2.add(ctx2);
						}
						
					}
				}
			equal_range.close();
			iter2.close();
			sorter1.cleanup();
			sorter1= null;
			
			sorter2.doneAdding();
			sorter2.setDestructiveIteration(true);

			iter2=sorter2.iterator();

			
			final VCFHeader header2 = new VCFHeader(header);				
			JVarkitVersion.getInstance().addMetaData(this, header2);
		
		
			final VCFInfoHeaderLine originalAlt=new VCFInfoHeaderLine("BND",1,VCFHeaderLineType.String,"Original ALT allele.");			
			header.addMetaDataLine(originalAlt);
			header.addMetaDataLine(new VCFInfoHeaderLine(originalAlt.getID(),1, VCFHeaderLineType.Integer, "SV length"));
			header.addMetaDataLine(new VCFInfoHeaderLine("SVLEN",1, VCFHeaderLineType.Integer, "SV length"));
			header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY, true));
			
			out.writeHeader(header2);
			while(iter2.hasNext()) {
				out.add(iter2.next());
				}
		    iter2.close();
		    sorter2.cleanup();
		    sorter2=null;
		    debug.flush();
		    debug.close();
		    debug=null;
		return 0;
	} catch(final Throwable err) {
		err.printStackTrace();
		LOG.error(err);
		return -1;
	} finally {
		if(sorter1!=null) sorter1.cleanup();
		if(sorter2!=null) sorter2.cleanup();
		if(debug!=null) try{debug.close();} catch(final Throwable err2) {}
	}
}
	
	public static void main(final String[] args) {
		new GridssPostProcessVcf().instanceMainWithExit(args);
	}

}
