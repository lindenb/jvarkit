package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.semontology.Term;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/*
BEGIN_DOC

## Example

```
$ java -jar dist/biostar251649.jar -n 10 -R tests/ref.fa tests/mutations.vcf
##INFO=<ID=SEQ3_10,Number=1,Type=String,Description="Sequence on the 3' of mutation">
##INFO=<ID=SEQ5_10,Number=1,Type=String,Description="Sequence on the 5' of mutation">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SEQ3_10=GATGGTAAGC;SEQ5_10=TCTACTCAGC;SGB=-61.9012;VDB=3.53678e-05	GT:PL	0/0:0,255,134	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SEQ3_10=GTTGTTGCTG;SEQ5_10=TTGAAGCTGC;SGB=-369.163;VDB=0.313337	GT:PL	0/0:0,255,133	0/1:40,0,31	0/0:0,255,134	0/0:0,255,82
```

END_DOC
*/
@Program(name="biostar251649",
	description=" Annotating the flanking bases of SNPs in a VCF file",
	biostars=251649,
	keywords={"vcf","annotation","sequence","reference"},
	terms = Term.ID_0000015
	)
public class Biostar251649 extends Launcher
	{
	private static final Logger LOG= Logger.build(Biostar251649.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names="-5",description="Left tag")
	private String leftTag="SEQ5_";
	@Parameter(names="-3",description="Right tag")
	private String rightTag="SEQ3_";
	@Parameter(names="-n",description="number of bases")
	private int extend=1;
	@Parameter(names={"-r","-R","--reference"},
			description=INDEXED_FASTA_REFERENCE_DESCRIPTION,
			required=true,
			converter=Launcher.IndexedFastaSequenceFileConverter.class
			)
	private IndexedFastaSequenceFile faidx = null;
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in,
			VariantContextWriter w)
		{
		try {
			final VCFHeader header = new VCFHeader(in.getHeader());
			VCFInfoHeaderLine info5 = new VCFInfoHeaderLine(leftTag+extend,
					1, VCFHeaderLineType.String,"Sequence on the 5' of mutation");
			VCFInfoHeaderLine info3 = new VCFInfoHeaderLine(rightTag+extend,
					1, VCFHeaderLineType.String,"Sequence on the 3' of mutation");
		
			if(header.getInfoHeaderLine(info5.getID())!=null)
				{
				LOG.error("tag "+info5.getID()+" already present in VCF header");
				return -1;
				}
			if(header.getInfoHeaderLine(info3.getID())!=null)
				{
				LOG.error("tag "+info3.getID()+" already present in VCF header");
				return -1;
				}
			
			header.addMetaDataLine(info5);
			header.addMetaDataLine(info3);
			GenomicSequence chrom= null;
			w.writeHeader(header);
			while(in.hasNext()) {
				final VariantContext ctx = in.next();
				if(chrom==null || !chrom.getChrom().equals(ctx.getContig())) {
					chrom = new GenomicSequence(this.faidx,ctx.getContig());
					}
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				if(ctx.getStart()>0)
					{
					final StringBuilder sb = new StringBuilder(this.extend);
					for(int i=0;i< this.extend;++i)
						{
						final int pos0 = (ctx.getStart()-2)-i;
						if(pos0<= 0) continue;
						sb.insert(0, chrom.charAt(pos0));
						}
					if(sb.length()>0) vcb.attribute(info5.getID(),sb.toString());
					}
				
					{
					final StringBuilder sb = new StringBuilder(this.extend);
					for(int i=0;i< this.extend;++i)
						{
						int pos0 = ctx.getEnd()+i;
						if(pos0>= chrom.length()) break;
						sb.append(chrom.charAt(pos0));
						}
					if(sb.length()>0) vcb.attribute(info3.getID(),sb.toString());
					}
				w.add(vcb.make());
				}
			
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(faidx);
			}	
		}
	@Override
	public int doWork(List<String> args)
		{
		return doVcfToVcf(args, outputFile);
		}
	
	public static void main(String[] args)
		{
		new Biostar251649().instanceMainWithExit(args);
		}
	}
