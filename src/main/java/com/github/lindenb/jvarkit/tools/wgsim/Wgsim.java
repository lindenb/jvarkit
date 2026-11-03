package com.github.lindenb.jvarkit.tools.wgsim;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriterFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.Launcher.WritingBamArgs;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMFileHeader;

/*


*/

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
@Program(name="wgsim",
description="generate fastqs from VCF",
keywords={"fastq","vcf","fasta"},
creationDate="20261103",
modificationDate="20261103"
)
public class Wgsim extends Launcher{
	private static final Logger LOG = Logger.of(Wgsim.class);
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	final Path faidx=null;
	@Parameter(names={"--output","-o"},description="fastq R1 output (or standout)")
	final Path output=null;
	@Parameter(names={"--read-length"},description="read length")
	final int read_length = 150;
	@Parameter(names={"--insert-size"},description="insert length")
	final int insert_length = 500;
	@Parameter(names={"--seed"},description="random seed generator -1: current time")
	final long seed= -1L;

	@ParametersDelegate
	protected WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	private static abstract class IntervalSeq extends DelegateCharSequence {
		protected Interval interval;
		IntervalSeq(final Interval r,CharSequence delegate)
			{
			super(delegate);
			interval=r;
			}
		boolean in_interval(int index0) {
			return index0+1>=this.interval.getStart() && index0+1<= interval.getEnd();
			}
		}
	
	private static class InversionSeq extends IntervalSeq {
		InversionSeq(final Interval r,CharSequence delegate)
			{
			super(r,delegate);
			}
		public char charAt(int index0) {
			if(!in_interval(index0)) {
				return getDelegate().charAt(index0);
				}
			else
				{
				char c= getDelegate().charAt((interval.getEnd()-1)-index0);
				return AcidNucleics.complement(c);
				}
			}
		}
	
	private static class DeletionSeq extends IntervalSeq {
		DeletionSeq(final Interval r,CharSequence delegate)
			{
			super(r,delegate);
			}
		public int length() {
			return getDelegate().length()-super.interval.length();
			}
		public char charAt(int index0) {
			if(index0+1 >= super.interval.getStart()) {
				index0+=super.interval.getLengthOnReference();
				}
			return getDelegate().charAt(index0);
			}
		}
	
	
	private static class DuplicationSeq extends IntervalSeq {
		DuplicationSeq(final Interval r,CharSequence delegate)
			{
			super(r,delegate);
			}
		public int length() {
			return getDelegate().length()+super.interval.length();
			}
		public char charAt(int index0) {
			if(index0+1 <=  this.interval.getEnd() ) {
				return getDelegate().charAt(index0);
				}
			return getDelegate().charAt(index0-super.interval.getLengthOnReference());
			}
		}
	
	private static class SnvSeq extends DelegateCharSequence {
		final VariantContext ctx;
		SnvSeq(final VariantContext ctx,CharSequence delegate)
			{
			super(delegate);
			this.ctx = ctx;
			}
		public String getAlt() {
			return ctx.getAlleles().get(1).getBaseString();
			}
		public int length() {
			return getDelegate().length() - 
					ctx.getReference().length() +
					getAlt().length()
					;
			}
		public char charAt(int index0) {
			if(index0+1 <  ctx.getStart()) {
				return getDelegate().charAt(index0);
				}
			index0-=(ctx.getStart()-1);
			final String alt = getAlt();
			if(index0<alt.length()) {
				return alt.charAt(index0); 
				}
			index0-= alt.length();
			return getDelegate().charAt(ctx.getStart()-1+index0);
			}
		}
	private void writeReads(
			SAMFileWriter samWriter,
			CharSequence genomic,
			Locatable interval
			)
		{
		
		}
		
	
	@Override
	public int doWork(List<String> args) {
		try {
			final Random random = new Random(seed<0?System.currentTimeMillis():seed);
			final String input=oneFileOrNull(args);
			final FastqPairedWriterFactory fqw = new FastqPairedWriterFactory();
			try(ReferenceSequenceFile reference=ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(reference);
				final IntervalTreeMap<List<VariantContext>> variants= new IntervalTreeMap<>();
				try( VCFIterator iter= new VCFIteratorBuilder().open(input)) {
					final VCFHeader header = iter.getHeader();
					while(iter.hasNext()) {
						final VariantContext vc1 = iter.next();
						final Interval r= new Interval(vc1);
						List<VariantContext> L = variants.get(r);
						if(L==null) {
							L=new ArrayList<>();
							variants.put(r, L);
							}
						L.add(vc1);
						}
					}
				try(OutputStream os = super.openPathOrStdoutAsStream(output)) {
					final SAMFileHeader samHeader=new SAMFileHeader(dict);
					samHeader.setSortOrder(SAMFileHeader.SortOrder.unknown);
					try(SAMFileWriter sfw = this.writingBamArgs.
							setReferencePath(this.faidx).
							openSamWriter(this.output, samHeader,true);) {
					
						for(SAMSequenceRecord ssr: dict.getSequences()) {
							if(!variants.containsOverlapping(ssr)) continue;
							final GenomicSequence genomic = new GenomicSequence(reference, ssr.getContig());
							
							
							}
						}
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
		}
	}
	
	public static void main(final String[] args) {
		new Wgsim().instanceMainWithExit(args);
	}

	
}
