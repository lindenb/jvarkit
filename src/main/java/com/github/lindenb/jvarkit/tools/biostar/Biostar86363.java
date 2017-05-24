package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.semontology.Term;


/**
BEGIN_DOC

## Example
```bash
$ cat reset.txt
20	14370	NA00001
20	1234567	NA00003
20	1110696	NA00002

$ curl "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/example-4.1.vcf" |\
  java -jar dist/biostar86363.jar -G reset.txt 

##fileformat=VCFv4.1
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GR,Number=1,Type=Integer,Description="(1) = Genotype was reset by Biostar86363:Set genotype of specific sample/genotype comb to unknown in multisample vcf fi
le. See http://www.biostars.org/p/86363/">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##fileDate=20090805
##phasing=partial
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##source=myImputationProgramV3.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	rs6054257	G	A	29	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:DP:GQ:GR:HQ	.|.:1:48:1:51,51	1|0:8:48:0:51,51	1/1:5:43:0
20	17330	.	T	A	3	q10	AF=0.017;DP=11;NS=3	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
20	1110696	rs6040355	A	G,T	67	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:DP:GQ:GR:HQ	1|2:6:21:0:23,27	.|.:0:2:1:18,2	2/2:4:35:0
20	1230237	.	T	.	47	PASS	AA=T;DP=13;NS=3	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
20	1234567	microsat1	GTC	G,GTCT	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ:GR	0/1:4:35:0	0/2:2:17:0	./.:3:40:1
```

END_DOC

*/
@Program(name="biostar86363",
	biostars=86363,
	keywords={"sample","genotype","vcf"},
			terms=Term.ID_0000015,
	description="Set genotype of specific sample/genotype comb to unknown in multisample vcf file. See http://www.biostars.org/p/86363/")
public class Biostar86363 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar86363.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	
	
	private Map<ContigPosRef,Set<String>> pos2sample=new HashMap<ContigPosRef, Set<String>>();

	private Biostar86363()
		{
		}

	@Parameter(names="-G",description="genotypes to reset. Format :CHROM(tab)POS(tab)ref(tab)SAMPLE. REQUIRED.",required=true)
	private File genotypeFile=null;
	
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) 
		{
		final List<Allele> empty_g=new ArrayList<Allele>(2);
		empty_g.add(Allele.NO_CALL);
		empty_g.add(Allele.NO_CALL);
		VCFHeader h=in.getHeader();
		final List<String> vcf_samples=h.getSampleNamesInOrder();
		h.addMetaDataLine(new VCFFormatHeaderLine("GR", 1, VCFHeaderLineType.Integer, "(1) = Genotype was reset by "+getProgramName()));
		out.writeHeader(h);
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			ContigPosRef cap=new ContigPosRef(ctx);
			Set<String> samplesToReset=this.pos2sample.get(cap);
			if(samplesToReset!=null)
				{
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				List<Genotype> genotypes=new ArrayList<Genotype>();
				for(String sample:vcf_samples)
					{
					Genotype g=ctx.getGenotype(sample);
					if(g==null) continue;
					GenotypeBuilder gb=new GenotypeBuilder(g);
					
					if(samplesToReset.contains(sample))
						{
						gb.alleles(empty_g);
						gb.attribute("GR",1);
						}
					else
						{
						gb.attribute("GR",0);
						}
					g=gb.make();
					genotypes.add(g);
					
					}
				vcb.genotypes(genotypes);
				ctx=VCFUtils.recalculateAttributes(vcb.make());
				}
			out.add(ctx);
			}
		return 0;
		}

	
	@Override
	public int doWork(List<String> args) {
		if(genotypeFile==null)
			{
			LOG.error("undefined genotype file");
			return -1;
			}
		BufferedReader in=null;
		try
			{
			Pattern tab=Pattern.compile("[\t]");
			in=IOUtils.openFileForBufferedReading(this.genotypeFile);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=tab.split(line);
				if(tokens.length<4)
					{
					LOG.error("Bad line in "+line);
					in.close();
					in=null;
					return -1;
					}
				ContigPosRef cap=new ContigPosRef(
					tokens[0],
					Integer.parseInt(tokens[1]),
					Allele.create(tokens[2],true)
					);
				Set<String> samples=this.pos2sample.get(cap);
				if(samples==null)
					{
					samples=new HashSet<String>();
					this.pos2sample.put(cap,samples);
					}
				samples.add(tokens[3]);
				}
			in.close();
			return doVcfToVcf(args, outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(in);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar86363().instanceMainWithExit(args);
		}

	}
