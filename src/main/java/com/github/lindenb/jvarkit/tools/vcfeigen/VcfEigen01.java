/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

* 2016 creation

*/
package com.github.lindenb.jvarkit.tools.vcfeigen;


import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;





import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC


### About Eigen

"Eigen is a spectral approach to the functional annotation of genetic variants in coding and noncoding regions. Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance. Eigen is an unsupervised approach, and, unlike most existing methods, is not based on any labelled training data. Eigen produces estimates of predictive accuracy for each functional annotation score, and subsequently uses these estimates of accuracy to derive the aggregate functional score for variants of interest as a weighted linear combination of individual annotations. We show that the resulting meta-score has good discriminatory ability using disease associated and putatively benign variants from published studies (for both Mendelian and complex diseases). The Eigen score is particularly useful in prioritizing likely causal variants in a region of interest when it is combined with population-level genetic data in the framework of a hierarchical model. Furthermore, an important advantage of the Eigen score is that it can be easily adapted to a specific tissue or cell type. More information about the Eigen score can be found in the accompanying manuscript: [A spectral approach integrating functional genomic annotations for coding and noncoding variants (Iuliana Ionita-Laza, Kenneth McCallum, Bin Xu, Joseph Buxbaum).A spectral approach integrating functional genomic annotations for coding and noncoding variants (Iuliana Ionita-Laza, Kenneth McCallum, Bin Xu, Joseph Buxbaum).](http://www.columbia.edu/~ii2135/Eigen_11_24.pdf) 

### Example

```

$ gunzip -c input.vcf.gz | cut -f 1-8 | java -jar dist/vcfeigen01.jar -D eigen_data_dir  | grep -i eigen

##INFO=<ID=EIGEN_CODING_Consequence,Number=A,Type=String,Description="Consequence">
##INFO=<ID=EIGEN_CODING_Eigen_PC_phred,Number=A,Type=Float,Description="Eigen-PC-phred">
##INFO=<ID=EIGEN_CODING_Eigen_PC_raw,Number=A,Type=Float,Description="Eigen-PC-raw">
##INFO=<ID=EIGEN_CODING_Eigen_phred,Number=A,Type=Float,Description="Eigen-phred">
##INFO=<ID=EIGEN_CODING_Eigen_raw,Number=A,Type=Float,Description="Eigen-raw">
##INFO=<ID=EIGEN_CODING_GERP_NR,Number=A,Type=Float,Description="GERP_NR">
##INFO=<ID=EIGEN_CODING_GERP_RS,Number=A,Type=Float,Description="GERP_RS">
##INFO=<ID=EIGEN_CODING_MA,Number=A,Type=Float,Description="MA">
##INFO=<ID=EIGEN_CODING_PhastPla,Number=A,Type=Float,Description="PhastPla">
##INFO=<ID=EIGEN_CODING_PhastPri,Number=A,Type=Float,Description="PhastPri">
##INFO=<ID=EIGEN_CODING_PhastVer,Number=A,Type=Float,Description="PhastVer">
##INFO=<ID=EIGEN_CODING_PhyloPla,Number=A,Type=Float,Description="PhyloPla">
##INFO=<ID=EIGEN_CODING_PhyloPri,Number=A,Type=Float,Description="PhyloPri">
##INFO=<ID=EIGEN_CODING_PhyloVer,Number=A,Type=Float,Description="PhyloVer">
##INFO=<ID=EIGEN_CODING_PolyPhenDIV,Number=A,Type=Float,Description="PolyPhenDIV">
##INFO=<ID=EIGEN_CODING_PolyPhenVar,Number=A,Type=Float,Description="PolyPhenVar">
##INFO=<ID=EIGEN_CODING_SIFT,Number=A,Type=Float,Description="SIFT">
##INFO=<ID=EIGEN_NC_DnasePval,Number=A,Type=Float,Description="DnasePval">
##INFO=<ID=EIGEN_NC_DnaseSig,Number=A,Type=Float,Description="DnaseSig">
##INFO=<ID=EIGEN_NC_Eigen_PC_phred,Number=A,Type=Float,Description="Eigen-PC-phred">
##INFO=<ID=EIGEN_NC_Eigen_PC_raw,Number=A,Type=Float,Description="Eigen-PC-raw">
##INFO=<ID=EIGEN_NC_Eigen_phred,Number=A,Type=Float,Description="Eigen-phred">
##INFO=<ID=EIGEN_NC_Eigen_raw,Number=A,Type=Float,Description="Eigen-raw">
##INFO=<ID=EIGEN_NC_FairePval,Number=A,Type=Float,Description="FairePval">
##INFO=<ID=EIGEN_NC_FaireSig,Number=A,Type=Float,Description="FaireSig">
##INFO=<ID=EIGEN_NC_GERP_NR,Number=A,Type=Float,Description="GERP_NR">
##INFO=<ID=EIGEN_NC_GERP_RS,Number=A,Type=Float,Description="GERP_RS">
##INFO=<ID=EIGEN_NC_H3K27ac,Number=A,Type=Float,Description="H3K27ac">
##INFO=<ID=EIGEN_NC_H3K4Me1,Number=A,Type=Float,Description="H3K4Me1">
##INFO=<ID=EIGEN_NC_H3K4Me3,Number=A,Type=Float,Description="H3K4Me3">
##INFO=<ID=EIGEN_NC_OCPval,Number=A,Type=Float,Description="OCPval">
##INFO=<ID=EIGEN_NC_PhastPla,Number=A,Type=Float,Description="PhastPla">
##INFO=<ID=EIGEN_NC_PhastPri,Number=A,Type=Float,Description="PhastPri">
##INFO=<ID=EIGEN_NC_PhastVer,Number=A,Type=Float,Description="PhastVer">
##INFO=<ID=EIGEN_NC_PhyloPla,Number=A,Type=Float,Description="PhyloPla">
##INFO=<ID=EIGEN_NC_PhyloPri,Number=A,Type=Float,Description="PhyloPri">
##INFO=<ID=EIGEN_NC_PhyloVer,Number=A,Type=Float,Description="PhyloVer">
##INFO=<ID=EIGEN_NC_PolIIPval,Number=A,Type=Float,Description="PolIIPval">
##INFO=<ID=EIGEN_NC_PolIISig,Number=A,Type=Float,Description="PolIISig">
##INFO=<ID=EIGEN_NC_TFBS_max,Number=A,Type=Float,Description="TFBS_max">
##INFO=<ID=EIGEN_NC_TFBS_num,Number=A,Type=Float,Description="TFBS_num">
##INFO=<ID=EIGEN_NC_TFBS_sum,Number=A,Type=Float,Description="TFBS_sum">
##INFO=<ID=EIGEN_NC_cmycPval,Number=A,Type=Float,Description="cmycPval">
##INFO=<ID=EIGEN_NC_cmycSig,Number=A,Type=Float,Description="cmycSig">
##INFO=<ID=EIGEN_NC_ctcfPval,Number=A,Type=Float,Description="ctcfPval">
##INFO=<ID=EIGEN_NC_ctcfSig,Number=A,Type=Float,Description="ctcfSig">
(...)
14	741	.	A	C	1	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0;EIGEN_NC_Eigen_PC_phred=7.69892;EIGEN_NC_Eigen_PC_raw=-0.05363416;EIGEN_NC_Eigen_phred=1.64761;EIGEN_NC_Eigen_raw=-0.27013662;EIGEN_NC_FairePval=0.0;EIGEN_NC_FaireSig=0.0;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=62.64;EIGEN_NC_H3K4Me1=9.0;EIGEN_NC_H3K4Me3=13.04;EIGEN_NC_OCPval=0.0;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.0;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0;SNVHPOL=4;SNVSB=0.0
14	1142	.	T	A	19	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.9238;EIGEN_NC_Eigen_PC_raw=0.7665708;EIGEN_NC_Eigen_phred=3.47555;EIGEN_NC_Eigen_raw=-0.10196207;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=56.68;EIGEN_NC_H3K4Me1=7.48;EIGEN_NC_H3K4Me3=13.0;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=0.0
14	1195	.	T	A	2	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.9238;EIGEN_NC_Eigen_PC_raw=0.7665708;EIGEN_NC_Eigen_phred=3.47555;EIGEN_NC_Eigen_raw=-0.10196207;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=56.68;EIGEN_NC_H3K4Me1=7.48;EIGEN_NC_H3K4Me3=13.0;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=0.0
14	1437	.	T	A	1	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.9751;EIGEN_NC_Eigen_PC_raw=0.7803033;EIGEN_NC_Eigen_phred=3.53009;EIGEN_NC_Eigen_raw=-0.09754433;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=61.56;EIGEN_NC_H3K4Me1=7.64;EIGEN_NC_H3K4Me3=14.64;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=0.0
14	1659	.	G	C	8	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=14.1281;EIGEN_NC_Eigen_PC_raw=0.82222354;EIGEN_NC_Eigen_phred=3.69067;EIGEN_NC_Eigen_raw=-0.08474564;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=78.24;EIGEN_NC_H3K4Me1=7.56;EIGEN_NC_H3K4Me3=18.84;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=3;SNVSB=0.0
14	1910	.	T	C	58	.	EIGEN_NC_DnasePval=0.0;EIGEN_NC_DnaseSig=0.0059;EIGEN_NC_Eigen_PC_phred=13.8695;EIGEN_NC_Eigen_PC_raw=0.75230044;EIGEN_NC_Eigen_phred=3.41565;EIGEN_NC_Eigen_raw=-0.106857374;EIGEN_NC_FairePval=5.04;EIGEN_NC_FaireSig=0.033;EIGEN_NC_GERP_NR=0.0;EIGEN_NC_GERP_RS=0.0;EIGEN_NC_H3K27ac=52.16;EIGEN_NC_H3K4Me1=5.84;EIGEN_NC_H3K4Me3=12.16;EIGEN_NC_OCPval=3.94;EIGEN_NC_PolIIPval=0.0;EIGEN_NC_PolIISig=0.004;EIGEN_NC_TFBS_max=0.0;EIGEN_NC_TFBS_num=0.0;EIGEN_NC_TFBS_sum=0.0;EIGEN_NC_cmycPval=0.0;EIGEN_NC_cmycSig=0.0254;EIGEN_NC_ctcfPval=0.0;EIGEN_NC_ctcfSig=0.0107;SNVHPOL=4;SNVSB=-10.8


```

END_DOC
 */
@Program(name="vcfeigen",
description="Annotator for the data of https://xioniti01.u.hpc.mssm.edu/v1.1/ : Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance.",
keywords={"vcf","annotation","variant"})
public class VcfEigen01
	extends Launcher
	{	
	private static final Logger LOG = Logger.build(VcfEigen01.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	
	@XmlType(name="vcfeigen")
	@XmlRootElement(name="vcfeigen")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			private class CtxWriter extends DelegateVariantContextWriter
				{
				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				@Override
				public void writeHeader(final VCFHeader header) {
					final  VCFHeader h2 = new VCFHeader(header);
					//addMetaData(h2);
					for(final VCFInfoHeaderLine vihl: CtxWriterFactory.this.annotator.getInfoHeaderLines()) {
						if(h2.getInfoHeaderLine(vihl.getID())!=null) {
							throw new JvarkitException.DuplicateVcfHeaderInfo(header,vihl.getID());
							}
						h2.addMetaDataLine(vihl);
						}
					super.writeHeader(h2);
					}
				@Override
				public void add(final VariantContext ctx) {
					final Map<String,Object> m  =  CtxWriterFactory.this.annotator.getAnnotations(ctx);
					if(m==null || m.isEmpty())
						{
						super.add(ctx);
						}
					else
						{
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						for(final String key: m.keySet()) {
							vcb.attribute(key, m.get(key));
							}
						super.add(vcb.make());
						}
					}
				}
			
			@Parameter(names={"-D","--directory"},
						description="Eigen directory containing the tabix files *.tab.gz",
						required=true
						)
			@XmlElement(name="directory")
			private File eigenDir = null;
			
			@Parameter(names={"-p","--prefix","--tabixFilePrefix"},
					description="prefix of the files in the eigen directory"
					)
			@XmlElement(name="prefix")
			private String tabixFilePrefix = "Eigen_hg19_";

			
			@XmlTransient
			private EigenInfoAnnotator annotator = null;
			
			
			@Override
			public int initialize() {
				try
					{
					LOG.info("loading eigen directory "+this.eigenDir);
					this.annotator = new EigenInfoAnnotator(this.eigenDir);
					this.annotator.setTabixFilePrefix(this.tabixFilePrefix);
					}
				catch(final Exception err) {
					LOG.error(err);
					return -1;
					}
				return 0;
				}
			
			@Override
			public VariantContextWriter open(final VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			@Override
			public void close() throws IOException {
				this.annotator = null;
				}
			}
	
	public VcfEigen01() {
		}

	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter delegate
			) {	
		final VariantContextWriter out = this.component.open(delegate);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(iter.getHeader()).logger(LOG);
		out.writeHeader(iter.getHeader());
		while(iter.hasNext())
			{
			out.add(progress.watch(iter.next()));
			}
		out.close();
		progress.finish();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.component.initialize()!=0) return -1;
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}

	
	public static void main(String[] args) throws Exception
		{
		new VcfEigen01().instanceMainWithExit(args);
		}

	}
