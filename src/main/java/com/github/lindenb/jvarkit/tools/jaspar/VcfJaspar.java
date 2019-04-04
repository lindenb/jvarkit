/**
 * 
 */
package com.github.lindenb.jvarkit.tools.jaspar;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import htsjdk.variant.vcf.VCFIterator;
/**

## Example

```bash

$ java -jar dist/vcfjaspar.jar  \
   -J jaspar.pfm \
   -R ref.fasta input.vcf.gz |\
grep -E '(JASPA|CHROM)' | cut -f 1-8


##INFO=<ID=JASPAR,Number=.,Type=String,Description="Jaspar pattern overlapping: Format: (Name|length|Score/1000|pos|strand)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
21	16133900	.	G	A	999	.	AC1=1;AF1=0.0009901;DP=111730;DP4=55553,55260,56,71;FQ=999;JASPAR=(MA0007.1/Arnt::Ahr|6|957|16133899|+);MQ=50;PV4=0.18,1.5e-08,5.5e-107,1;RPB=6.741643e+00;VDB=4.316638e-01
21	16147540	.	G	A	999	.	AC1=2;AF1=0.00198;DP=134303;DP4=33057,100000,135,399;FQ=999;JASPAR=(MA0189.1/NFE2L1::MafG|6|986|16147540|-);MQ=50;PV4=0.8,0.01,0,1;RPB=1.275421e+01;VDB=5.201277e-01
21	16149050	.	G	A	999	.	AC1=1;AF1=0.0009901;DP=379132;DP4=184618,191640,225,217;FQ=999;JASPAR=(MA0003.1/Arnt::Ahr|6|957|16149045|-);MQ=49;PV4=0.45,7.9e-08,9.1e-178,1;RPB=8.890954e+00;VDB=9.818940e-11
21	16155525	.	G	A	145	.	AC1=1;AF1=0.0009901;DP=145035;DP4=81479,55899,197,36;FQ=145;JASPAR=(MA0589.1/KLF5|10|1000|16155525|+);MQ=49;PV4=1.2e-16,9.2e-71,4.3e-229,1;RPB=2.359170e+00;VDB=5.181676e-01
21	16156637	.	G	A	999	.	AC1=1;AF1=0.0009901;DP=319931;DP4=222355,92779,313,22;FQ=999;JASPAR=(MA0003.1/Arnt|6|1000|16156632|+);MQ=49;PV4=2.4e-25,2.9e-28,1.1e-171,1;RPB=1.627989e+01;VDB=7.589921e-21
21	16156868	.	G	A	131	.	AC1=2;AF1=0.002035;DP=12301;DP4=895,11327,1,14;FQ=131;JASPAR=(MA0461.1/Atoh1|8|979|16156866|+);MQ=49;PV4=1,0.26,0.00026,1;RPB=5.751797e+00;VDB=2.956324e-03
21	16158562	.	A	G	999	.	AC1=619;AF1=0.6129;DP=555112;DP4=158228,47130,256836,81463;FQ=999;JASPAR=(MA1089.1/NFE2L1::MafG|6|986|16158559|-);MQ=47;PV4=1.8e-21,0,0,1;RPB=-4.773974e+02;VDB=6.345308e-27

```


 */
@Program(name="vcfjaspar",
description="Finds JASPAR profiles in VCF",
keywords={"vcf","matrix","jaspar"})
public class VcfJaspar extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfJaspar.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names="-J",description=" jaspar PFM uri. required. example: http://jaspar2016.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt",required=true)
	private String jasparUri=null;
	@Parameter(names="-f",description="(0<ratio<1) fraction of best score")
	private double fraction_of_max=0.95;
	@Parameter(names={"-R","-r","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File fasta = null;
	@Parameter(names={"-T","--tag"},description="VCF tag")
	private String ATT = "JASPAR";

	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private List<Matrix> jasparDb=new ArrayList<Matrix>();
	

	public VcfJaspar() {
		}
	
	@Override
	protected int doVcfToVcf(final String inputName,final VCFIterator in,final VariantContextWriter out) {
		GenomicSequence genomicSequence=null;
		final VCFHeader header=new VCFHeader(in.getHeader());
		addMetaData(header);

		header.addMetaDataLine(new VCFInfoHeaderLine(this.ATT,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String, 
				"Jaspar pattern overlapping: Format: (Name|length|Score/1000|pos|strand)"
				));
		JVarkitVersion.getInstance().addMetaData(getProgramName(), header);
		out.writeHeader(header);
		while(in.hasNext())
			{
			VariantContext var=in.next();

			if(genomicSequence==null || !genomicSequence.getChrom().equals(var.getContig()))
				{
				LOG.info("Loading sequence "+var.getContig());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,var.getContig());
				}
			
			final Set<String> hits=new HashSet<String>();
		
			for(final Matrix matrix:this.jasparDb)
				{
					int start0=Math.max(0, var.getStart() - matrix.length());
					for(int y=start0;y<var.getStart() && y+matrix.length() <= genomicSequence.length();++y)
						{
						final CharSequence forward=new SubSequence(genomicSequence,y,y+matrix.length());
						final CharSequence revcomp=new RevCompCharSequence(forward);
						
						//run each strand
						for(int strand=0;strand<2;++strand)
							{
							double score= matrix.score(strand==0?forward:revcomp);
							if(score<=0) continue;
							
							if(score>= matrix.max()*this.fraction_of_max)
								{
								StringBuilder b=new StringBuilder("(");
								b.append(matrix.getName().replaceAll("[ \t;=]+", "/"));
								b.append("|");
								b.append(matrix.length());
								b.append("|");
								b.append((int)(1000.0*(score/matrix.max())));
								b.append("|");
								b.append(y+1);
								b.append("|");
								b.append(strand==0?'+':'-');
								b.append(")");
								hits.add(b.toString());
								break;
								}
							}
						
					}
				}
			if(hits.isEmpty())
				{
				out.add(var);
				continue;
				}
			final VariantContextBuilder vcb=new VariantContextBuilder(var);
			vcb.attribute(ATT, hits.toArray(new String[hits.size()]));
			out.add(vcb.make());
			}
		return RETURN_OK;
		}
	@Override
	public int doWork(final List<String> args) {
		
		if(this.jasparUri==null)
			{
			LOG.error("Undefined jaspar-uri");
			return -1;
			}		
		
		if(this.fasta==null)
			{
			LOG.error("Undefined fasta sequence");
			return -1;
			}
		try
			{
			LOG.info("Reading JASPAR: "+jasparUri);
			final LineIterator liter = IOUtils.openURIForLineIterator(this.jasparUri);
			final Iterator<Matrix> miter=Matrix.iterator(liter);
			while(miter.hasNext())
				{
				final Matrix matrix = miter.next();
				this.jasparDb.add(matrix.convertToPWM());
				}
			CloserUtil.close(liter);
			LOG.info("JASPAR size: "+this.jasparDb.size());

			if(jasparDb.isEmpty())
				{
				LOG.warn("JASPAR IS EMPTY");
				}
			
			LOG.info("opening "+fasta);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(fasta);
			return doVcfToVcf(oneFileOrNull(args),outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.indexedFastaSequenceFile=null;
			}
		}

	public static void main(String[] args)
		{
		new VcfJaspar().instanceMainWithExit(args);
		}
	}
