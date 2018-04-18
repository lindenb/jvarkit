/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

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

two columns:

  * 1st column: REF chromosome name, multiple/aliases can be separated by a pipe '|'
  * 2st column: ancestral contig name
  * 3nd column: File PATH to the equivalent ancestral contig.


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
			final Pattern tab =Pattern.compile("[\t]");
			final Pattern pipe =Pattern.compile("[\\|]");
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
			final VcfIterator iterin,
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
