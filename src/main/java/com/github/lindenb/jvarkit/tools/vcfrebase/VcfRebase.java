/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

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
public class VcfRebase {
	private static final Logger LOG = Logger.build().prefix("vcfRebase").make();
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private Rebase rebase=Rebase.createDefaultRebase();
	@Parameter(names={"-A","--attribute"},description="VCF INFO attribute")
	private String ATT="ENZ";
	public VcfRebase() {
		}
	
	public  VcfRebase indexedFastaSequenceFile(final IndexedFastaSequenceFile ref) {
		this.indexedFastaSequenceFile=ref;
		return this;
	}
	public  VcfRebase tag(final String tag) {
		this.ATT=tag;
		return this;
	}
	public  VcfRebase rebase(final  Rebase reb) {
		this.rebase=reb;
		return this;
	}

	public VariantContextWriter open(final VariantContextWriter delegate)
		{
		if(this.indexedFastaSequenceFile==null)
			throw new JvarkitException.UserError("Indexed fasta file was not provided.");
		if(this.rebase==null)
			throw new JvarkitException.UserError("Rebase missing.");
		return new Worker(delegate,
				this.indexedFastaSequenceFile,
				this.rebase,
				this.ATT
				);
		}
	
	private static class Worker extends DelegateVariantContextWriter
		{
		private	final String ATT;
		private GenomicSequence genomicSequence=null;
		private final IndexedFastaSequenceFile indexedFastaSequenceFile;
		private final Rebase rebase;
		Worker( final VariantContextWriter delegate,
				final IndexedFastaSequenceFile indexedFastaSequenceFile,
				final Rebase rebase,
				final String att
				) {
			super(delegate);
			this.indexedFastaSequenceFile=indexedFastaSequenceFile;
			this.rebase=rebase;
			this.ATT=att;
			}
		@Override
		public void writeHeader(final VCFHeader header) {
			final VCFHeader header2=new VCFHeader(header);
			header2.addMetaDataLine(new VCFInfoHeaderLine(ATT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Enzyme overlapping: Format: (Name,Site,Sequence,pos-1,strand)"));
			super.writeHeader(header2);
			}
		@Override
		public void add(VariantContext var)
			{
			if(genomicSequence==null || !genomicSequence.getChrom().equals(var.getContig()))
				{
				LOG.info("Current contig "+var.getContig());
				this.genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,var.getContig());
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
						for(x=0;x< enz.size() && y+x < this.genomicSequence.length() ;++x )
							{
							char c=(strand==0?
									enz.at(x):
									AcidNucleics.complement(enz.at((enz.size()-1)-x))
									);
							if(!Rebase.compatible(genomicSequence.charAt(y+x),c)) break;
							}
						// match found
						if(x==enz.size())
							{
							StringBuilder b=new StringBuilder("(");
							b.append(enz.getName());
							b.append("|");
							b.append(enz.getDecl());
							b.append("|");
							for(x=0;x < enz.size();++x)
								{
								char c=genomicSequence.charAt(y+x);
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
				super.add(var);
				}
			else
				{
				VariantContextBuilder vcb=new VariantContextBuilder(var);
				vcb.attribute(ATT, hits.toArray(new String[hits.size()]));
				super.add(vcb.make());
				}
			}
		
		}
		
	public static class Launcher extends com.github.lindenb.jvarkit.util.jcommander.Launcher
		{
		@ParametersDelegate
		private VcfRebase instance=new VcfRebase();
		@Parameter(names={"-R","-reference"},description="Reference fasta file")
		File referenceFile;
		@Parameter(names={"-E","-enzyme"},description="restrict to that enzyme name")
		private Set<String> selEnzymesStr = new HashSet<>();
		@Parameter(names={"-w","-weight"},description="min enzyme weight")
		private float weight= 5f;		
		@Parameter(names={"-o","--out"},required=false,description="Output vcf , ot stdin")
		private VariantContextWriter out=new VcfWriterOnDemand();

		@Override
		public int doWork(List<String> args) {
			Rebase rebase=Rebase.createDefaultRebase();

			if(!this.selEnzymesStr.isEmpty())
				{
				Rebase rebase2=new Rebase();
				for(String e:selEnzymesStr)
					{
					if(e.isEmpty()) continue;
					Rebase.Enzyme enz=rebase.getEnzymeByName(e);
					if(enz==null)
						{
						StringBuilder msg= new StringBuilder();
						msg.append("Cannot find enzyme \""+e +"\" in RE list.\n");
						msg.append("Current list is:\n");
						for(Rebase.Enzyme E: rebase)
							{
							msg.append("\t"+E+"\n");
							}
						throw new JvarkitException.UserError(msg.toString());
						}
					rebase2.getEnzymes().add(enz);
					}
				rebase=rebase2;
				}
		
			int i=0;
			while(i< rebase.size())
				{
				if(rebase.get(i).getWeight()< weight)
					{
					rebase.getEnzymes().remove(i);
					}
				else
					{
					++i;
					}
				}
		
			if(rebase.size()==0)
				{
				LOG.warn("REBASE IS EMPTY");
				}
	
			if(referenceFile==null)
				{
				throw new JvarkitException.ReferenceMissing("reference.undefined");
				}
			IOUtil.assertFileIsReadable(referenceFile);
			IndexedFastaSequenceFile indexedFastaSequenceFile=null;
			VcfIterator in=null;
			VariantContextWriter out=null;
			try
				{
				this.instance.
					rebase(rebase).
					indexedFastaSequenceFile(indexedFastaSequenceFile);
				indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.referenceFile);
				in = VCFUtils.createVcfIterator(super.oneFileOrNull(args));
				out=instance.open(this.out);
				VCFUtils.copyHeaderAndVariantsTo(in,out);
				in.close();
				out.close();
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException(err);
				}
			finally
				{
				CloserUtil.close(indexedFastaSequenceFile);
				CloserUtil.close(out);
				CloserUtil.close(in);
				}
			return super.doWork(args);
			}
		}
	
	public static void main(String[] args)
		{
		new VcfRebase.Launcher().instanceMainWithExit(args);
		}
	}
