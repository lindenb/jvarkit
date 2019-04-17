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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfrebase;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.PostponedVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;

/**
 BEGIN_DOC
 
 ## Example
 

 ```
 $  curl -s  "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz" |\
   gunzip -c  |\
   java -jar dist/vcfrebase.jar -R  human_g1k_v37.fasta -w 8   |\
   grep -E '(#|ENZ=)' 

##fileformat=VCFv4.1
##INFO=<ID=ENZ,Number=.,Type=String,Description="Enzyme overlapping: Format: (Name,Site,Sequence,pos-1,strand)">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	256138	rs372711491	T	G	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAAT|256139|+),(SwaI|ATTT^AAAT|ATTTAAAt|256131|+);OTHERKG;RS=372711491;RSPOS=256138;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=138
1	744076	rs142181107	A	C	.	.	CAF=[0.9775,0.0225];COMMON=1;ENZ=(SwaI|ATTT^AAAT|ATTTAAaT|744070|+);KGPROD;KGPhase1;RS=142181107;RSPOS=744076;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100014000100;WGT=1;dbSNPBuildID=134
1	762592	rs71507462	C	G	.	.	CAF=[0.2755,0.7245];COMMON=1;ENZ=(MauBI|CG^CGCGCG|cGCGCGCG|762592|+);GNO;KGPROD;KGPhase1;OTHERKG;R5;RS=71507462;RSPOS=762592;SAO=0;SLO;SSR=0;VC=SNV;VP=0x050100020001100116000100;WGT=1;dbSNPBuildID=130
1	780347	rs202219272	TTTAA	T	.	.	CAF=[0.9867,0.01331];COMMON=1;ENZ=(PacI|TTAAT^TAA|ttaaTTAA|780348|+);INT;KGPROD;KGPhase1;KGPilot123;OTHERKG;RS=202219272;RSPOS=780348;SAO=0;SSR=0;VC=DIV;VP=0x05000008000110001e000200;WGT=1;dbSNPBuildID=137
1	780361	rs111333032	A	T	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAaT|780355|+);GNO;INT;OTHERKG;RS=111333032;RSPOS=780361;SAO=0;SLO;SSR=0;VC=SNV;VP=0x050100080001000102000100;WGT=1;dbSNPBuildID=132
1	786891	rs183914415	C	T	.	.	CAF=[0.9752,0.02479];COMMON=1;ENZ=(AscI|GG^CGCGCC|GGCGCGCC|786892|+);INT;KGPROD;KGPhase1;RS=183914415;RSPOS=786891;SAO=0;SSR=0;VC=SNV;VP=0x050000080001100014000100;WGT=1;dbSNPBuildID=135
1	820332	rs201213736	A	G,T	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAaT|820326|+);OTHERKG;RS=201213736;RSPOS=820332;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=137
1	820333	rs201439577	T	C	.	.	ENZ=(SwaI|ATTT^AAAT|ATTTAAAt|820326|+);OTHERKG;RS=201439577;RSPOS=820333;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=137
1	822560	rs200342299	C	T	.	.	ENZ=(SfiI|GGCCNNNN^NGGCC|GGCCAGcTTGGCC|822554|+);OTHERKG;RS=200342299;RSPOS=822560;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000002000100;WGT=1;dbSNPBuildID=137
1	822565	rs111427246	C	T	.	.	ENZ=(SfiI|GGCCNNNN^NGGCC|GGCCAGCTTGGcC|822554|+);GNO;OTHERKG;RS=111427246;RSPOS=822565;SAO=0;SLO;SSR=0;VC=SNV;VP=0x050100000001000102000100;WGT=1;dbSNPBuildID=132
1	837753	rs187865648	G	A	.	.	CAF=[0.9995,0.0004591];COMMON=0;ENZ=(SfiI|GGCCNNNN^NGGCC|gGCCATTCTGGCC|837753|+);KGPROD;KGPhase1;RS=187865648;RSPOS=837753;SAO=0;SSR=0;VC=SNV;VP=0x050000000001000014000100;WGT=1;dbSNPBuildID=135
1	839873	rs192553893	C	T	.	.	CAF=[0.6768,0.3232];COMMON=1;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|839872|+);KGPROD;KGPhase1;OTHERKG;RS=192553893;RSPOS=839873;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100016000100;WGT=1;dbSNPBuildID=135
1	839911	rs76652930	C	T	.	.	ASP;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|839910|+);GNO;OTHERKG;RS=76652930;RSPOS=839911;SAO=0;SSR=0;VC=SNV;VP=0x050000000005000102000100;WGT=1;dbSNPBuildID=131
1	839912	rs369394889	GGCC	G	.	.	ENZ=(NotI|GC^GGCCGC|GCggccGC|839910|+);OTHERKG;RS=369394889;RSPOS=839913;SAO=0;SSR=0;VC=DIV;VP=0x050000000001000002000200;WGT=1;dbSNPBuildID=138
1	839933	rs146045242	C	T	.	.	CAF=[0.8852,0.1148];COMMON=1;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|839932|+);KGPROD;KGPhase1;OTHERKG;RS=146045242;RSPOS=839933;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100016000100;WGT=1;dbSNPBuildID=134
1	840009	rs140080750	C	T	.	.	CAF=[0.9692,0.03076];COMMON=1;ENZ=(NotI|GC^GGCCGC|GcGGCCGC|840008|+);KGPROD;KGPhase1;OTHERKG;RS=140080750;RSPOS=840009;SAO=0;SSR=0;VC=SNV;VP=0x050000000001100016000100;WGT=1;dbSNPBuildID=134
```
 
 END_DOC
 
 */
@SuppressWarnings("unused")
@Program(name="vcfrefbase",
	description="Restriction sites overlaping variations in a vcf",
	keywords={"vcf","rebase","restriction","enzyme"})
public class VcfRebase extends Launcher {
	private static final Logger LOG = Logger.build(VcfRebase.class).make();
	@Parameter(names={"-o","--out"},required=false,description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@ParametersDelegate
	private PostponedVariantContextWriter.WritingVcfConfig writingVcfArgs = new PostponedVariantContextWriter.WritingVcfConfig();
	@Parameter(names={"-A","--attribute"},description="VCF INFO attribute")
	private String ATT="ENZ";
	@Parameter(names={"-R","-reference","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File referenceFile;
	@Parameter(names={"-E","-enzyme","--enzyme"},description="restrict to that enzyme name. Default: use all enzymes")
	private Set<String> selEnzymesStr = new HashSet<>();
	@Parameter(names={"-w","-weight","--weight"},description="min enzyme weight 6 = 6 cutter like GAATTC, 2 = 2 cutter like ATNNNNNNAT  ")
	private float weight= 5f;		

	
	private Rebase rebase = Rebase.createDefaultRebase();
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private GenomicSequence genomicSequence=null;			
	private ContigNameConverter contigNameConverter;
	
	public VcfRebase() {
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			) {	
		final VCFHeader header= iter.getHeader();
		
		

		final VCFHeader header2=new VCFHeader(header);
		header2.addMetaDataLine(
				new VCFInfoHeaderLine(this.ATT,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Enzyme overlapping: Format: (Name,Site,Sequence,pos-1,strand)")
					);
		this.genomicSequence = null;
		
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();

		out.writeHeader(iter.getHeader());
		while(iter.hasNext())
			{
			final VariantContext var = progress.apply(iter.next());
			final String refContig = this.contigNameConverter.apply(var.getContig());
			if(StringUtils.isBlank(refContig))
				{
				out.add(var);
				continue;
				}
			
			if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(refContig))
				{
				this.genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,refContig);
				}
			
			final Set<String> hits=new HashSet<String>();
			for(final Rebase.Enzyme enz:this.rebase)
				{
				int start0=Math.max(0, var.getStart() - enz.size());
				for(int y=start0;y<=var.getStart();++y)
					{
					//run each strand
					for(int strand=0;strand<2;++strand)
						{
						int x=0;
						//loop over bases of the enzyme
						for(x=0;x< enz.size() && y+x <  this.genomicSequence.length() ;++x )
							{
							final char c=(strand==0?
									enz.at(x):
									AcidNucleics.complement(enz.at((enz.size()-1)-x))
									);
							if(!Rebase.compatible(this.genomicSequence.charAt(y+x),c)) break;
							}
						// match found
						if(x==enz.size())
							{
							final StringBuilder b=new StringBuilder("(");
							b.append(enz.getName());
							b.append("|");
							b.append(enz.getDecl());
							b.append("|");
							for(x=0;x < enz.size();++x)
								{
								char c= this.genomicSequence.charAt(y+x);
								if(y+x>=var.getStart()-1 && y+x<=var.getEnd()-1)
									{
									c=Character.toLowerCase(c);
									}
								b.append(c);
								}
							b.append("|");
							b.append(y+1);
							b.append("|");
							b.append(strand==0?"+":"-");
							b.append(")");
							hits.add(b.toString());
							break;
							}
						if(enz.isPalindromic()) break;
						}
					}
				}
			if(hits.isEmpty())
				{
				out.add(var);
				}
			else
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(var);
				vcb.attribute( this.ATT, hits.toArray(new String[hits.size()]) );
				out.add(vcb.make());
				}
			}
		out.close();
		progress.close();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(!this.selEnzymesStr.isEmpty())
				{
				final Rebase rebase2=new Rebase();
				for(final String e:this.selEnzymesStr)
					{
					if(StringUtils.isBlank(e)) continue;
					final Rebase.Enzyme enz= this.rebase.getEnzymeByName(e);
					if(enz==null)
						{
						final StringBuilder msg= new StringBuilder();
						msg.append("Cannot find enzyme \""+e +"\" in RE list.\n");
						msg.append("Current list is:\n");
						for(final Rebase.Enzyme E: rebase)
							{
							msg.append("\t"+E+"\n");
							}
						LOG.error(msg.toString());
						return -1;
						}
					rebase2.getEnzymes().add(enz);
					}
				this.rebase=rebase2;
				}
		
			int i=0;
			while(i< this.rebase.size())
				{
				if(this.rebase.get(i).getWeight()< weight)
					{
					this.rebase.getEnzymes().remove(i);
					}
				else
					{
					++i;
					}
				}
		
			if(this.rebase.size()==0)
				{
				LOG.warn("REBASE IS EMPTY");
				}
			
			IOUtil.assertFileIsReadable(this.referenceFile);
			
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.referenceFile);
			this.contigNameConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(this.referenceFile));
			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		}
	
	public static void main(final String[] args)
		{
		new VcfRebase().instanceMainWithExit(args);
		}
	}
