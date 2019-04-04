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
*/
package com.github.lindenb.jvarkit.tools.biostar;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**

BEGIN_DOC


### Example

```

$ java -jar dist-2.0.1/biostar175929.jar -x 2 -R ~/src/gatk-ui/testdata/ref.fa -b ~/src/gatk-ui/testdata/S1.vcf.gz  | more

>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(A)|rotavirus:961-961(T)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[A]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[T]tg
>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(A)|rotavirus:961-961(A)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[A]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[A]tg
>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(T)|rotavirus:961-961(T)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[T]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[T]tg
>rotavirus:127|rotavirus:130-130(T)|rotavirus:232-232(T)|rotavirus:267-267(C)|rotavirus:286-286(G)|rotavirus:536-536(A)|rotavirus:693-693(T)|rotavirus:833-833(G)|rotavirus:916-916
(T)|rotavirus:961-961(A)
at[T]caatatgattacaatgaagtatttaccagagttaaaagtaaatttgattatgtga
tggatgactctggtgttaaaaacaatcttttgggtaaagctataac[T]attgatcaggc
gttaaatggaaagtttagctcag[C]tattagaaatagaaattg[G]atgactgattcta
aaacggttgctaaattagatgaagacgtgaataaacttagaatgactttatcttctaaag
ggatcgaccaaaagatgagagtacttaatgcttgttttagtgtaaaaagaataccaggaa
aatcatcatcaataattaaatgcactagacttatgaaggataaaatagaacgtggagaag
ttgaggttgatgattcatatgttgatgagaaaatggaaattgatactattgattgg[A]a
atctcgttatgatcagttagaaaaaagatttgaatcactaaaacagagggttaatgagaa
atacaatacttgggtacaaaaagcgaagaaagtaaatgaaaatatgtactctcttcagaa
tgttatctcacaacagcaaaaccaaatagcagatc[T]tcaacaatattgtagtaaattg
gaagctgatttgcaaggtaaatttagttcattagtgtcatcagttgagtggtatctaagg
tctatggaattaccagatgatgtaaagaatgacattgaacagcagttaaattcaatt[G]
atttaattaatcccattaatgctatagatgatatcgaatcgctgattagaaatttaattc
aagattatgacagaacattttt[T]atgttaaaaggactgttgaagcaatgcaactatga
atatgcata[A]tg

```



END_DOC
*/


@Program(name="biostar175929",
	description="Construct a combination set of fasta sequences from a vcf",
	biostars=175929,
	keywords={"fasta","vcf"}
	)
public class Biostar175929 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar175929.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-x","--extend"},description="extend FASTA sequence by 'n' bases")
	private int extendBases = 100 ;

	@Parameter(names={"-b","--bracket"},description="Surround variant with '[' and ']'")
	private boolean bracket = false;
	
	@Parameter(names={"-R","--reference"},description="indexed Fasta reference",required=true)
	private File faidx = null;

	
	private PrintWriter pw;
	
	private void recursive(
			final GenomicSequence chromosome,
			final List<VariantContext> variants,
			int index,
			StringBuilder title,
			StringBuilder sequence
			
			)
		{
		if(index==variants.size())
			{
			int lastPos0= (variants.get(index-1).getEnd()-1);
			for(int i=0;i< this.extendBases && lastPos0+i< chromosome.length();++i)
				{
				sequence.append(Character.toLowerCase(chromosome.charAt(lastPos0+i)));
				}
			pw.print(">");
			pw.print(title);
			for(int i=0;i< sequence.length();++i)
				{
				if(i%60==0) pw.println();
				pw.print(sequence.charAt(i));
				}
			pw.println();
			return;
			}
		if(index==0)
			{
			int firstPos0= (variants.get(0).getStart()-1);
			int chromStart=Math.max(0, firstPos0-extendBases);
			title.append(variants.get(0).getContig()+":"+chromStart);
			for(int i=Math.max(0, firstPos0-extendBases);i< firstPos0 ;++i)
				{
				sequence.append(Character.toLowerCase(chromosome.charAt(i)));
				}
			}
		else
			{
			int endPos0= (variants.get(index-1).getEnd());
			int begPos0= (variants.get(index).getStart()-1);
			while(endPos0< begPos0)
				{
				sequence.append(Character.toLowerCase(chromosome.charAt(endPos0)));
				endPos0++;
				}
			}
		final int  title_length= title.length();
		final int  sequence_length= sequence.length();
		final VariantContext ctx = variants.get(index);
		for(final Allele allele: ctx.getAlleles())
			{
			if(allele.isNoCall()) continue;
			if(allele.isSymbolic())  continue;
			title.setLength(title_length);
			sequence.setLength(sequence_length);
			title.append("|"+ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+"("+allele.getBaseString()+")");
			if(this.bracket) sequence.append("[");
			sequence.append(allele.getBaseString().toUpperCase());
			if(this.bracket) sequence.append("]");
			recursive(chromosome, variants, index+1, title, sequence);
			}
		
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.faidx==null) 
			{
			LOG.error("fasta reference was not defined.");
			return -1;
			}
		IndexedFastaSequenceFile reference = null;
		VCFIterator iter=null;
		try {
			reference = new IndexedFastaSequenceFile(this.faidx);
			iter = super.openVCFIterator(oneFileOrNull(args));
			this.pw = openFileOrStdoutAsPrintWriter(this.outputFile);
			final List<VariantContext> variants = new ArrayList<>();
			for(;;)
				{
				VariantContext ctx = null;
				if(iter.hasNext()) {
					ctx = iter.next();
				}
				
				if( ctx == null || (!variants.isEmpty() && !ctx.getContig().equals(variants.get(0).getContig()))) {
					if(!variants.isEmpty()) 
						{
						LOG.info("chrom:" +variants.get(0).getContig()+ " N="+variants.size());
						final GenomicSequence genomic = new GenomicSequence(reference,variants.get(0).getContig());
						final StringBuilder title= new StringBuilder();
						final StringBuilder sequence= new StringBuilder();
						recursive(genomic,variants,0,title,sequence);
						variants.clear();
						}
					if( ctx == null) break;
					}
				variants.add(ctx);
				}
			iter.close();iter=null;
			this.pw.flush();
			this.pw.close();
			return RETURN_OK;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		}
			finally {
				CloserUtil.close(reference);
				CloserUtil.close(iter);
				CloserUtil.close(pw);
			}
		}
	
	public static void main(final String[] args)throws Exception
		{
		new Biostar175929().instanceMainWithExit(args);
		}
	}
