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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfrebase;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
 BEGIN_DOC
 
 ## Example
 

 ```
$ java -jar dist/vcfrebase.jar -w 6 -R ~/data/human_g1k_v37.fasta src/test/resources/test_vcf01.vcf | bcftools annotate -x '^INFO/ENZ' | bcftools view --drop-genotypes | grep ENZ

##INFO=<ID=ENZ,Number=.,Type=String,Description="Enzyme overlapping: Format: (Name|Site|Sequence|pos|strand)">
##bcftools_annotateCommand=annotate -x ^INFO/ENZ; Date=Wed Nov 13 10:38:39 2019
1	852063	.	G	A	387	PASS	ENZ=PflMI|CCANNNN^NTGG|CCAGGCCCTGG|852064|+
1	866893	.	T	C	431	PASS	ENZ=SacI|GAGCT^C|GAGCtC|866889|+
1	875770	.	A	G	338	PASS	ENZ=ClaI|AT^CGAT|ATCGaT|875766|+
1	909238	.	G	C	229	PASS	ENZ=PmaCI|CAC^GTG|CACgTG|909235|+
1	913889	.	G	A	372	PASS	ENZ=BsaXI|(9/12)ACNNNNNCTCC(10/7)|GGAGGCCCCgT|913880|-
1	918384	.	G	T	489	PASS	ENZ=DraIII|CACNNN^GTG|CACgCCGTG|918381|+
1	933790	.	G	A	436	PASS	ENZ=BsaXI|(9/12)ACNNNNNCTCC(10/7)|GGAGGAGGGgT|933781|-
1	940005	.	A	G	188	PASS	ENZ=GsuI|CTGGAG(16/14)|CTGGAG|940006|+,BaeI|(10/15)ACNNNNGTAYC(12/7)|GGTaCTGGAGT|940002|-
1	940096	.	C	T	487	PASS	ENZ=BcgI|(10/12)CGANNNNNNTGC(12/10)|cGAGGTGGGTGC|940096|+
1	950113	.	GAAGT	G	1427	PASS	ENZ=Eco57I|CTGAAG(16/14)|CTgaag|950111|+
1	950243	.	A	C	182	PASS	ENZ=BclI|T^GATCA|TGaTCA|950241|+
1	951283	.	C	T	395	PASS	ENZ=NarI|GG^CGCC|GGcGCC|951281|+
1	951564	.	A	G	105	PASS	ENZ=BstXI|CCANNNNN^NTGG|CCaAGTAGTTGG|951562|+
1	952003	.	G	A	177	PASS	ENZ=Bpu10I|CCTNAGC(-5/-2)|CCTCAGC|952004|+,BbvCI|CCTCAGC(-5/-2)|CCTCAGC|952004|+
1	952428	.	G	A	456	PASS	ENZ=EciI|GGCGGA(11/9)|TCCgCC|952425|-
1	953952	.	G	A	490	PASS	ENZ=BsrDI|GCAATG(2/0)|CATTgC|953948|-
1	959155	.	G	A	370	PASS	ENZ=BarI|(7/12)GAAGNNNNNNTAC(12/7)|gAAGCCGCTCTAC|959155|+
1	959231	.	G	A	350	PASS	ENZ=BsaXI|(9/12)ACNNNNNCTCC(10/7)|GGAGGGTCCgT|959222|-
1	960409	.	G	C	357	PASS	ENZ=BseYI|CCCAGC(-5/-1)|CCCAgC|960405|+
1	962210	.	A	G	300	PASS	ENZ=NcoI|C^CATGG|CCaTGG|962208|+
1	964389	.	C	T	32	LowGQXHetSNP;LowGQXHomSNP	ENZ=BseYI|CCCAGC(-5/-1)|cCCAGC|964389|+
1	967658	.	C	T	515	PASS	ENZ=StuI|AGG^CCT|AGGCcT|967654|+
1	970215	.	G	C	379	PASS	ENZ=DrdI|GACNNNN^NNGTC|GACCCCTCGGTC|970216|+
1	972180	.	G	A	403	PASS	ENZ=AgeI|A^CCGGT|ACCgGT|972177|+
1	1004957	.	G	A	316	PASS	ENZ=BsgI|GTGCAG(16/14)|gTGCAG|1004957|+
1	1004980	.	G	A	292	PASS	ENZ=BsePI|G^CGCGC|gCGCGC|1004980|+
1	1011087	.	CG	C	1052	PASS	ENZ=Eam1105I|GACNNN^NNGTC|GACTCTCAGTc|1011077|+
1	1017170	.	C	G	507	PASS	ENZ=AloI|(7/12)GAACNNNNNNTCC(12/7)|GAACAGAGcATCC|1017162|+
```
 
 END_DOC
 
 */
@Program(name="vcfrebase",
	description="Restriction sites overlaping variations in a vcf",
	keywords={"vcf","rebase","restriction","enzyme"},
	creationDate="20131115",
	modificationDate="20200624"
	)
public class VcfRebase extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfRebase.class).make();
	@Parameter(names={"-A","--attribute"},description="VCF INFO attribute")
	private String ATT="ENZ";
	@Parameter(names={"-R","-reference","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path referenceFile;
	@Parameter(names={"-E","-enzyme","--enzyme"},description="restrict to that enzyme name. Default: use all enzymes")
	private Set<String> selEnzymesStr = new HashSet<>();
	@Parameter(names={"-w","-weight","--weight"},description="min enzyme weight 6 = 6 cutter like GAATTC, 2 = 2 cutter like ATNNNNNNAT  ")
	private float weight= 5f;		
	
	
	public VcfRebase() {
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iter, VariantContextWriter out) {
		 ReferenceSequenceFile indexedFastaSequenceFile=null;
		 GenomicSequence genomicSequence=null;			

		try
			{
			Rebase rebase = Rebase.createDefaultRebase();

			
			if(!this.selEnzymesStr.isEmpty())
				{
				final Rebase rebase2=new Rebase();
				for(final String e:this.selEnzymesStr)
					{
					if(StringUtils.isBlank(e)) continue;
					final Rebase.Enzyme enz= rebase.getEnzymeByName(e);
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
			
			IOUtil.assertFileIsReadable(this.referenceFile);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceFile);
			indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referenceFile);
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
			
			
			
			
			final VCFHeader header= iter.getHeader();
			
			

			final VCFHeader header2=new VCFHeader(header);
			header2.addMetaDataLine(
					new VCFInfoHeaderLine(this.ATT,
							VCFHeaderLineCount.UNBOUNDED,
							VCFHeaderLineType.String,
							"Enzyme overlapping: Format: (Name|Site|Sequence|pos|strand)")
						);
			

			out.writeHeader(header2);
			while(iter.hasNext())
				{
				final VariantContext var = iter.next();
				final String refContig = contigNameConverter.apply(var.getContig());
				if(StringUtils.isBlank(refContig))
					{
					out.add(var);
					continue;
					}
				
				if(genomicSequence==null || !genomicSequence.getChrom().equals(refContig))
					{
					genomicSequence=new GenomicSequence(indexedFastaSequenceFile,refContig);
					}
				
				final Set<String> hits=new HashSet<String>();
				for(final Rebase.Enzyme enz: rebase)
					{
					int start0=Math.max(0, var.getStart() - enz.size());
					for(int y=start0;y<=var.getStart();++y)
						{
						//run each strand
						for(int strand=0;strand<2;++strand)
							{
							int x=0;
							//loop over bases of the enzyme
							for(x=0;x< enz.size() && y+x <  genomicSequence.length() ;++x )
								{
								final char c=(strand==0?
										enz.at(x):
										AcidNucleics.complement(enz.at((enz.size()-1)-x))
										);
								if(!Rebase.compatible(genomicSequence.charAt(y+x),c)) break;
								}
							// match found
							if(x==enz.size())
								{
								final StringBuilder b=new StringBuilder();
								b.append(enz.getName());
								b.append("|");
								b.append(enz.getDecl());
								b.append("|");
								for(x=0;x < enz.size();++x)
									{
									char c= genomicSequence.charAt(y+x);
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
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			}
		}
	
	public static void main(final String[] args)
		{
		new VcfRebase().instanceMainWithExit(args);
		}
	}
