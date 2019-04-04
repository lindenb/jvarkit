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
package com.github.lindenb.jvarkit.tools.onekgenomes;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
BEGIN_DOC

## Motivation

For @isamtalves : annotate a VCF for ancestral allele using data from 1000 genomes: 

  * http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README

## Manifest

manifest is a tab delimited file.

empty lines are ignored

lines starting with '#' are ignored

three columns:

  * 1st column: REF chromosome name, multiple/aliases can be separated by a pipe '|'
  * 2nd column: ancestral contig name
  * 3td column: File PATH to the equivalent ancestral contig.

### note to self: building the manifest

```bash
 find ${PWD} -name "*.fa" | sort -V | while read F; do echo "$F" | sed 's%/commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_%%' | sed 's/.fa//' | tr -d "\n" && echo -ne "\t" &&  head -n 1 $F | cut -c 2- | tr -d '\n' && echo -en "\t" && echo $F  ; done | sed 's/^\([^\t]*\)/\1|chr\1/'
```

```
1|chr1    ANCESTOR_for_chromosome:GRCh37:1:1:249250621:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_1.fa
2|chr2    ANCESTOR_for_chromosome:GRCh37:2:1:243199373:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_2.fa
3|chr3    ANCESTOR_for_chromosome:GRCh37:3:1:198022430:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_3.fa
4|chr4    ANCESTOR_for_chromosome:GRCh37:4:1:191154276:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_4.fa
5|chr5    ANCESTOR_for_chromosome:GRCh37:5:1:180915260:1   /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/human_ancestor_5.fa
```

## Example:

```
$ java -jar dist/vcfancestralalleles.jar \
	-m /commun/data/pubdb/ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/human_ancestor_GRCh37_e59/manifest.mf \
	src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz |\
	bcftools annotate -x '^INFO/AA'

(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	905606	rs540662886	G	C,A	41743.9	PASS	AA=G
1	905608	rs770396126	G	T,A	1006.45	PASS	AA=G
1	905609	.	G	A	694.04	PASS	AA=G
1	905610	rs775689041	TG	AG,T,TGGGGGGCCCAG	4327.1	PASS	AA=TG
1	905611	rs749351425	G	T	2434.35	PASS	AA=G
1	905616	.	G	C	1801.88	PASS	AA=G
1	905617	rs376988925	C	T,G	4379.6	PASS	AA=C
1	905619	rs774441222	C	T	19350.2	PASS	AA=C
1	905621	rs368876607	G	A	14291.5	PASS	AA=G
1	905623	rs770778738	G	A	8255.09	PASS	AA=g
1	905626	.	G	A	257.48	PASS	AA=G
1	905627	rs755913930	GC	G,GCC	75797.7	PASS	AA=GC
1	905628	rs776376856	C	A	76618.4	AC0;RF	AA=C
1	905629	rs759517260	C	T	13583.4	PASS	AA=C
1	905632	.	C	T	2152.62	PASS	AA=C
1	905634	rs765105613	C	T	4534.75	PASS	AA=C
1	905635	rs752637312	C	A	4947.89	PASS	AA=C
1	905636	.	C	T	888.44	PASS	AA=C
1	905637	rs762369764	C	T	313.28	PASS	AA=C
(...)
```

END_DOC 
 */
@Program(name="vcfancestralalleles",
description="Annotate a VCF with it's ancestral allele. Data from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.README",
keywords={"vcf","sort"}
)
public class VcfAncestralAllele
extends Launcher
	{
	private static final  Logger LOG = Logger.build(VcfAncestralAllele.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-m","--manifest"},description="Manifest file containing the path to the fasta files. See doc. ALL fasta files must be indexed with `samtools faidx`", required=true)
	private File manifestFile = null;
	@Parameter(names={"-t","--tag"},description="Ancestral allele INFO attribute")
	private String aaAttribute = "AA";

	private final Map<String,Mapping> contigName2ancestralFasta = new HashMap<>();
	
	private static class Mapping
		{
		final String ancestralContig;
		final File ancestralFile;
		Mapping(final String ancestralContig,final File ancestralFile) {
			this.ancestralContig = ancestralContig;
			this.ancestralFile = ancestralFile;
			}
		}
	
	VcfAncestralAllele() {
		
	}
	
	private void loadManifest() throws IOException
		{
		BufferedReader r = null;
		try {
			final CharSplitter tab = CharSplitter.TAB;
			final CharSplitter pipe = CharSplitter.PIPE;
			r = IOUtils.openFileForBufferedReading(this.manifestFile);
			String line;
			while((line=r.readLine())!=null) {
				if(line.isEmpty() || StringUtil.isBlank(line)) continue;
				final String tokens[] = tab.split(line);
				if(tokens.length<3) throw new JvarkitException.TokenErrors("expected two columns", tokens);
				final String ancestralContig  = tokens[1];
				final File fastaPath  = new File(tokens[2]);
				IOUtil.assertFileIsReadable(fastaPath);
				/* no, one fasta can contains more than one sequence
				if(contigName2ancestralFasta.values().contains(fastaPath)) {
					throw new IOException("fasta already defined for  "+fastaPath);
					}
				*/
				final File faiPath  = new File(tokens[2]+".fai");
				IOUtil.assertFileIsReadable(faiPath);
				for(final String sn:pipe.split(tokens[0])) {
					if( StringUtil.isBlank(line)) {
						throw new IOException("empty contig name in "+line);
						}
					if( this.contigName2ancestralFasta.containsKey(sn)) {
						throw new IOException("duplicate contig name in "+line);
						}
					this.contigName2ancestralFasta.put(sn, new Mapping(ancestralContig,fastaPath));
					}
				}
			if(this.contigName2ancestralFasta.isEmpty()) {
				LOG.warn("Manifest is empty "+this.manifestFile);
				}
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iterin,
			final VariantContextWriter out) {
		final Set<String> unmapped_contigs= new HashSet<>();
		IndexedFastaSequenceFile indexedFastaSequenceFile = null;
		try {
			final VCFInfoHeaderLine AAheaderLine = new VCFInfoHeaderLine(
					aaAttribute,
					1,
					VCFHeaderLineType.String,
					"Ancestral allele. Manifest file was "+this.manifestFile
					);
			final VCFHeader header = iterin.getHeader();
			if(header.getInfoHeaderLine(aaAttribute)!=null) {
				throw new JvarkitException.DuplicateVcfHeaderInfo(aaAttribute);
			}
			final VCFHeader header2 = new VCFHeader(header);
			header2.addMetaDataLine(AAheaderLine);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			out.writeHeader(header2);
			
			String prev_contig = null;
			Mapping prev_mapping = null;
	
			while(iterin.hasNext())
				{
				final VariantContext ctx = progress.watch(iterin.next());
				if(prev_contig==null || !prev_contig.equals(ctx.getContig()))
					{
					prev_contig = null;
					prev_mapping = null;
					CloserUtil.close(indexedFastaSequenceFile);
					indexedFastaSequenceFile = null;
					if(unmapped_contigs.contains(ctx.getContig())) {
						out.add(ctx);
						continue;
						}
					prev_mapping = this.contigName2ancestralFasta.get(ctx.getContig());
					if(prev_mapping==null) {
						LOG.warn("No mapping for contig "+ctx.getContig());
						unmapped_contigs.add(ctx.getContig());
						out.add(ctx);
						continue;
						}
					prev_contig = ctx.getContig();
					indexedFastaSequenceFile = new IndexedFastaSequenceFile(
							prev_mapping.ancestralFile
							);
					}
				
				
				String ancestral = null;
				try {
					final ReferenceSequence refseq = indexedFastaSequenceFile.getSubsequenceAt(
							prev_mapping.ancestralContig,
							ctx.getStart(),
							ctx.getEnd()
							);
					ancestral = refseq.getBaseString();
				} catch(final Exception err2) {
					LOG.warn(ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd()+" "+err2.getMessage());
					ancestral=VCFUtils.escapeInfoField("ERROR_"+err2.getMessage());
					}
				
				
				if(StringUtil.isBlank(ancestral)) {
					out.add(ctx);
					continue;
					}
				final VariantContextBuilder vcb= new VariantContextBuilder(ctx);
				vcb.attribute(AAheaderLine.getID(),ancestral);
				out.add(vcb.make());
				}
			out.close();
			progress.finish();
			if(!unmapped_contigs.isEmpty()) {
				LOG.warn("UNMAPPED CONTIGS :"+unmapped_contigs);
				}
			return 0;
			}
		catch(final Exception error) {
			LOG.error(error);
			return -1;
			}
		finally
			{
			CloserUtil.close(indexedFastaSequenceFile);
			}
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			loadManifest();
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	
	public static void main(final String[] args)
	{
	new VcfAncestralAllele().instanceMainWithExit(args);
	}
	
}
