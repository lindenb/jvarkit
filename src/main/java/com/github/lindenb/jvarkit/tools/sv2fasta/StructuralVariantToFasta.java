/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.sv2fasta;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**
BEGIN_DOC

# motivation

create a fasta input for `pggb` (pangenome graph builder)

# input

input is a set of indexed VCF files or a file with the suffix '.list' containing the path to the VCFs

# example

```
java -jar dist/jvarkit.jar sv2fasta --interval "chr1:204091-325320" --skip-no-sv --exclude-filtered -R ref.fasta vcfs.list > tmp.fa

samtools faidx tmp.fa

pggb -i tmp.fa \
	-o OUT \
    -n `grep ">" tmp.fa | wc -l` \
	-p 90 -s 100 \
	-t 5 1>&2

```


END_DOC
 */
@Program(
	name="sv2fasta",
	description="convert VCF of structural variant(s) to fasta for pggb",
	keywords={"vcf","cnv","fasta"},
	creationDate="20230403",
	modificationDate="20230403",
	jvarkit_amalgamion = true
	)
public class StructuralVariantToFasta extends Launcher {
	private static final Logger LOG=Logger.of(StructuralVariantToFasta.class);
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-r","--region","--interval"},description="interval. "+ IntervalParser.OPT_DESC,required=true)
	private String intervalStr = null;
	@Parameter(names={"--exclude-filtered"},description="Exclude FILTER-ed variants.")
	private boolean exclude_filtered = false;
	@Parameter(names={"--skip-no-sv"},description="Exclude VCF without structural variant")
	private boolean skip_no_sv = false;


	private static class GenomicSeq extends AbstractCharSequence{
		final SAMSequenceRecord ssr;
		final Locatable interval;
		final byte[] bases;
		GenomicSeq(final SAMSequenceRecord ssr,final Locatable interval,final byte[] bases) {
			this.ssr = ssr;
			this.interval = interval;
			this.bases= bases;
			}
		@Override
		public char charAt(int index) {
			return (char)Character.toUpperCase(bases[index - (this.interval.getStart()-1)]);
			}
		
		@Override
		public int length() {
			return ssr.getLengthOnReference();
			}
		}
	
	private CharSequence fetchFasta(
			final List<VariantContext> variants,
			final GenomicSeq genomicSeq,
			final int chromStart,
			final int chromEnd
			) throws IOException {
		
		
		final StringBuilder sb=new StringBuilder(CoordMath.getLength(chromStart, chromEnd));
		int pos = chromStart;
		while(pos < chromEnd) {
			VariantContext ctx = null;
			for(int i=0;i< variants.size();i++) {
				final VariantContext ctx2 = variants.get(i);
				if(ctx2.getStart() > pos) break;
				if(ctx2.getStart() == pos) {
					if(!CoordMath.encloses(chromStart, chromEnd, ctx2.getStart(), ctx2.getEnd())) continue;
					ctx = variants.remove(i);
					break;
					}
				}
			
			if(ctx==null) {
				sb.append(genomicSeq.charAt(pos-1));
				pos++;
				continue;
				}
			
			final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
			if(StringUtils.isBlank(svType)) {
				final Allele alt = ctx.getAlternateAlleles().stream().filter(A->AcidNucleics.isATGCN(A)).findFirst().orElse(null);
				if(alt!=null) {
					sb.append(alt.getDisplayString());
					pos = ctx.getEnd()+1;
					}
				else
					{
					sb.append(genomicSeq.charAt(pos-1));
					pos++;
					}
				}
			else if(svType.equals("BND"))
				{
				sb.append(genomicSeq.charAt(pos-1));
				pos++;
				}
			else if(svType.equals("DEL"))
				{
				pos = ctx.getEnd()+1;
				}
			else if(svType.equals("DUP"))
				{
				final CharSequence subseq = fetchFasta(
						variants,
						genomicSeq,
						ctx.getStart(),
						ctx.getEnd()
						);
				sb.append(subseq);// first
				sb.append(subseq);// twice
				pos = ctx.getEnd()+1;
				}
			else if(svType.equals("INV"))
				{
				final CharSequence subseq = fetchFasta(
						variants,
						genomicSeq,
						ctx.getStart(),
						ctx.getEnd()
						);
				pos = ctx.getEnd()+1;
				sb.append(AcidNucleics.reverseComplement(subseq));
				}
			else if(svType.equals("INS"))
				{
				final Allele alt = ctx.getAlternateAlleles().stream().filter(A->AcidNucleics.isATGCN(A)).findFirst().orElse(null);
				if(alt!=null) {
					sb.append(alt.getDisplayString());
					pos = ctx.getEnd()+1;
					}
				else if(ctx.hasAttribute("CONSENSUS")) /* delly */ {
					sb.append(ctx.getAttributeAsString("CONSENSUS",""));
					pos = ctx.getEnd()+1;
					}
				else if(ctx.hasAttribute("LEFT_SVINSSEQ") && ctx.hasAttribute("RIGHT_SVINSSEQ"))  /* manta */ {
					sb.append(ctx.getAttributeAsString("LEFT_SVINSSEQ",""));
					sb.append(StringUtils.repeat(50, 'N'));
					sb.append(ctx.getAttributeAsString("RIGHT_SVINSSEQ",""));
					pos = ctx.getEnd()+1;
					}
				else if(ctx.hasAttribute("DUPHOMSEQ") && ctx.hasAttribute("DUPHOMLEN")) /* manta */ {
					final String subseq = ctx.getAttributeAsString("DUPHOMSEQ","");
					final int homlen =  ctx.getAttributeAsInt("DUPHOMLEN",0);
					for(int i=0;i< homlen;i++) {
						sb.append(subseq);
						}
					pos = ctx.getEnd()+1;
					}
				else if(ctx.hasAttribute("HOMSEQ") && ctx.hasAttribute("HOMLEN")) /* manta */ {
					final String subseq = ctx.getAttributeAsString("HOMSEQ","");
					final int homlen =  ctx.getAttributeAsInt("HOMLEN",0);
					for(int i=0;i< homlen;i++) {
						sb.append(subseq);
						}
					pos = ctx.getEnd()+1;
					}
				else
					{
					sb.append(genomicSeq.charAt(pos-1));
					pos++;
					}
				}
			else
				{
				sb.append(genomicSeq.charAt(pos-1));
				pos++;
				}
			}
		
		return sb;
		}
	
	private void writeFasta(PrintWriter pw, String name, final CharSequence fasta) {
		pw.print(">");
		pw.print(name);
		for(int i=0;i< fasta.length();i++) {
			if(i%60==0) pw.println();
			pw.print(fasta.charAt(i));
		}
		pw.println();
	}
	
	private void applyVcf(
			final PrintWriter pw,
			final Path path,
			final GenomicSeq genomicSeq,
			final SAMSequenceDictionary dict
			) throws IOException {
		
		final int chromStart = genomicSeq.interval.getStart();
		final int chromEnd = genomicSeq.interval.getEnd();
		
		final Predicate<VariantContext> acceptVariant = V-> {
				if(!CoordMath.encloses(chromStart, chromEnd, V.getStart(), V.getEnd())) return false;
				if(exclude_filtered && V.isFiltered()) return false;
				return true;
				};
				
		try(VCFReader reader = VCFReaderFactory.makeDefault().open(path, true)) {
			final VCFHeader header = reader.getHeader();
			final SAMSequenceDictionary dict2 = SequenceDictionaryUtils.extractRequired(header);
			SequenceUtil.assertSequenceDictionariesEqual(dict, dict2);
			
			
			
			if(!header.hasGenotypingData()) {
				LOG.info(path.toString());
				
				final List<VariantContext> variants = new ArrayList<>();
				try(CloseableIterator<VariantContext> iter= reader.query(genomicSeq.ssr.getContig(),chromStart,chromEnd)) {
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						if(!acceptVariant.test(ctx)) continue;
						variants.add(new VariantContextBuilder(ctx).noGenotypes().make());
						}
					}
				
				if(skip_no_sv && variants.isEmpty()) return;
				
				final CharSequence seq =  fetchFasta(
						variants,
						genomicSeq, 
						chromStart,
						chromEnd
						);
				
				
				
				writeFasta(pw, IOUtils.getFilenameWithoutCommonSuffixes(path),seq);
				}
			else
				{
				for(final String sn: header.getGenotypeSamples()) {
					LOG.info(sn+" in "+path);
					final List<VariantContext> variants = new ArrayList<>();
					
					try(CloseableIterator<VariantContext> iter= reader.query(genomicSeq.ssr.getContig(),chromStart,chromEnd)) {
						while(iter.hasNext()) {
							final VariantContext ctx = iter.next();
							if(!acceptVariant.test(ctx)) continue;
							final Genotype gt = ctx.getGenotype(sn);
							if(gt.getAlleles().stream().allMatch(A->A.isNoCall() || A.isReference())) continue;
							variants.add(new VariantContextBuilder(ctx).noGenotypes().make());
							}
						}
					if(skip_no_sv && variants.isEmpty()) continue;
					
					final CharSequence seq =  fetchFasta(
							variants,
							genomicSeq,
							chromStart,
							chromEnd
							);
					writeFasta(pw, sn,seq);
					}
				}
			}
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final List<Path> paths= IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				LOG.error("No vcf was provided.");
				return -1;
				}
			try(ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(reference);
				final Locatable loc = new IntervalParser(dict).apply(this.intervalStr).orElse(null);
				if(loc==null) {
					LOG.error("Cannot parse interval "+loc);
					return -1;
					}
				final SAMSequenceRecord ssr = dict.getSequence(loc.getContig());
				final byte[] bases = reference.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getEnd()).getBases();
				final GenomicSeq genomicSeq = new GenomicSeq(ssr,loc,bases);
				
				
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					for(final Path path: paths) {
						applyVcf(pw, path,genomicSeq, dict);
						}
					pw.flush();
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new StructuralVariantToFasta().instanceMainWithExit(args);
	}
	
}
