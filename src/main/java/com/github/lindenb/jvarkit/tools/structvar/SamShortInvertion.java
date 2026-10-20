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
package com.github.lindenb.jvarkit.tools.structvar;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Motivation

finds the regions having some short inversions.


output is a VCF file

## Example:

```
$ java -jar ${JVARKIT_DIST}/samshortinvert.jar -R ref.fasta bams.list |\
 	bcftools view -i 'INFO/DPMAX>10' > out.vcf
```

END_DOC

**/

@Program(name="samshortinvert",
	description="Scan short inversions in reads.",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20261020",
	creationDate="20140228"
	)
public class SamShortInvertion extends Launcher
	{
	private static final Logger LOG = Logger.of(SamShortInvertion.class);	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFaidx = null;
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int mapq=1;
	@Parameter(names={"--min-inversion-length"},description="min inversion length")
	private int min_inv_length = 8;

	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();

	private static abstract class MappedSequence extends AbstractCharSequence {
		abstract CigarOperator getOperator();
		}
	
	private static class SubReadSequence extends MappedSequence {
		private final byte[] sequence;
		private final int start;
		private final int length;
		SubReadSequence( byte[] sequence,int start,int length) {
			this.sequence= sequence;
			this.start = start;
			this.length=length;
			}
		@Override
		public char charAt(int index) {
			return (char)sequence[start+index];
			}
		@Override
		public int length() {
			return this.length;
			}
		@Override
		CigarOperator getOperator() {
			return CigarOperator.I;
			}
		}

	private static class SubGenomeSequence extends MappedSequence {
		private final Locatable genomicPos;
		private final CharSequence seq;
		SubGenomeSequence(Locatable genomicPos,GenomicSequence genomicSeq) {
			this.genomicPos = genomicPos;
			this.seq = genomicSeq.subSequence(genomicPos); 
			}
		@Override
		public char charAt(int index) {
			return (char)seq.charAt(index);
			}
		@Override
		public int length() {
			return seq.length();
			}
		@Override
		CigarOperator getOperator() {
			return CigarOperator.D;
			}
		}

	
	private static class Inversion {
		String contig;
		int start;
		String sequence;
		int count=0;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input =oneFileOrNull(args);
			final SamReaderFactory srf = super.createSamReaderFactory()
					.validationStringency(ValidationStringency.LENIENT)
					.referenceSequence(this.referenceFaidx)
					;
			try(ReferenceSequenceFile referenceSequenceFile = new IndexedFastaSequenceFile(this.referenceFaidx)) {
				GenomicSequence genomicSequence= null;
				final SAMSequenceDictionary dict =SequenceDictionaryUtils.extractRequired(referenceSequenceFile);
				try(SamReader samReader =  srf.open(StringUtils.isBlank(input)?SamInputResource.of(stdin()):SamInputResource.of(input))) {
					final SAMFileHeader samHeader = samReader.getFileHeader();
					final String sampleName = samHeader.getReadGroups().
							stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().orElseThrow();
	
					
					final Set<VCFHeaderLine> meta=new HashSet<>();
					VCFStandardHeaderLines.addStandardFormatLines(meta,true,
							VCFConstants.GENOTYPE_KEY,
							VCFConstants.DEPTH_KEY
							);
					VCFStandardHeaderLines.addStandardInfoLines(meta,true,
							VCFConstants.DEPTH_KEY,
							VCFConstants.END_KEY
							
							);
					meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
					meta.add(new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Structural variant type"));
					meta.add(new VCFInfoHeaderLine("SEQ", 1, VCFHeaderLineType.String,"DNA Sequence"));
					
					
					final VCFHeader header=new VCFHeader(meta,Collections.singletonList(sampleName));
					JVarkitVersion.getInstance().addMetaData(this, header);
					header.setSequenceDictionary(dict);
					try(VariantContextWriter vcw = this.writingVariants.open(this.outputFile)) {
						vcw.writeHeader(header);
						final Map<Interval,Inversion> all_inversions = new TreeMap<>();
						String prevContig=null;
						try(CloseableIterator<SAMRecord> iter = samReader.iterator()) {
						for(;;)
							{
							final SAMRecord rec = iter.hasNext()?iter.next():null;
							if(rec!=null && !SAMRecordDefaultFilter.accept(rec,this.mapq)) continue;
							
							if(rec==null || (prevContig!=null && !rec.getContig().equals(prevContig))) {
								for(Inversion inv: all_inversions.values()) {
									final Allele alt = Allele.create("<INV>",false);
									final VariantContextBuilder vcb = new VariantContextBuilder(
											"", 
											inv.contig, 
											inv.start,
											inv.start + inv.sequence.length()-1,
											Arrays.asList(Allele.REF_N,alt)
											);
									vcb.attribute("SVLEN",inv.sequence.length());
									vcb.attribute("SEQ",inv.sequence);
									vcb.attribute(VCFConstants.END_KEY,inv.start + inv.sequence.length()-1);
									vcb.attribute("SVTYPE","INV");
									
									vcb.genotypes(new GenotypeBuilder(sampleName, Arrays.asList(Allele.REF_N,alt))
											.DP(inv.count)
											.make());
									
									vcw.add(vcb.make());
									}
								
								
								if(rec==null) {
									break;
									}
								prevContig = rec.getContig();
								all_inversions.clear();
								}
							
							
							final byte[] sequence = rec.getReadBases();
							if(SAMRecord.NULL_SEQUENCE.equals(sequence)) continue;
							final Cigar cigar = rec.getCigar();
							if(cigar.numCigarElements()< 4 ) continue;
							
							
							final List<MappedSequence> subSequences= new ArrayList<>();
							
							int refpos1 = rec.getAlignmentStart();
							int readpos0=0;
							for(CigarElement ce: cigar.getCigarElements()) {
								final CigarOperator op = ce.getOperator();
								switch(op) {
									case S: case H: break;
									case M: case X: case EQ: break;
									case N: case D: {
										if(genomicSequence==null || !genomicSequence.hasName(rec.getContig())) {
											genomicSequence = new GenomicSequence(referenceSequenceFile, rec.getContig());
											}
										final Locatable loc = new SimpleInterval(
												rec.getContig(),
												refpos1,
												refpos1+ce.getLength()-1
												);
										subSequences.add(new SubGenomeSequence(loc, genomicSequence));
										break;
										}
									case I: {
										subSequences.add(new SubReadSequence(sequence, readpos0, ce.getLength()));
										break;
										}
									case P:break;
									}

								if(ce.getOperator().consumesReadBases()) {
									readpos0+=ce.getLength();
									}
								if(ce.getOperator().consumesReferenceBases()) {
									refpos1+=ce.getLength();
									}
								} // end loop cigar elements
							
							for(int x=0;x+1< subSequences.size();++x) {
								final MappedSequence seqx = subSequences.get(x);
								if(seqx.length()< this.min_inv_length) continue;
								for(int y=x+1;y< subSequences.size();++y) {
									final MappedSequence seqy = subSequences.get(y);
									if(seqx.length()!=seqy.length()) continue;
									if(seqx.getOperator().equals(seqy.getOperator())) continue;
									if(!AcidNucleics.areReverseComplement(seqx,seqy)) continue;
									Locatable loc=null;
									if(seqx instanceof SubGenomeSequence) {
										loc = SubGenomeSequence.class.cast(seqx).genomicPos;
										}
									else if(seqy instanceof SubGenomeSequence) {
										loc = SubGenomeSequence.class.cast(seqy).genomicPos;
										}
									if(loc!=null) {
										final Interval r= new Interval(loc);
										Inversion inv = all_inversions.get(r);
										if(inv==null) {
											inv=new Inversion();
											inv.contig = loc.getContig();
											inv.start = loc.getStart();
											inv.sequence = seqx.toString();
											all_inversions.put(r,inv);
											}
										inv.count++;
										}
									}
								
								}
							
							}
							
							}
						}
					}
				}
			return 0;
			}
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		}

	
	public static void main(final String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
