/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.phased;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC


# Example
```
$ java -jar dist/bam2haplotypes.jar -V src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S5.bam

#CHROM	START	END	COUNT	N-VARIANTS	(POS\tALT)+
RF03	1221	1242	11	2	1221	C	1242	C
RF03	1688	1708	5	2	1688	G	1708	T
RF04	1900	1920	4	2	1900	C	1920	A
RF06	517	543	9	2	517	C	543	G
RF06	668	695	4	2	668	G	695	T
RF08	926	992	2	2	926	C	992	G
RF09	294	317	6	2	294	T	317	A
RF10	139	175	1	2	139	T	175	G
RF10	139	175	3	2	139	T	175	C
```

in **paired** mode

```
samtools collate -O -u src/test/resources/S5.bam TMP | java -jar dist/bam2haplotypes.jar --paired -V src/test/resources/rotavirus_rf.vcf.gz

#CHROM	START	END	COUNT	N-VARIANTS	(POS\tALT)+
RF02	251	578	1	2	251	A	578	G
RF03	1221	1688	1	2	1221	C	1688	G
RF03	1221	1242	7	2	1221	C	1242	C
RF03	1221	1688	1	2	1221	C	1688	G
RF03	1688	1708	1	2	1688	G	1708	T
RF03	1708	2150	1	2	1708	T	2150	T
RF03	1221	1708	1	3	1221	C	1688	G	1708	T
RF03	1221	1688	2	3	1221	C	1242	C	1688	G
RF03	1688	2150	1	3	1688	G	1708	T	2150	T
RF03	1221	1708	2	4	1221	C	1242	C	1688	G	1708	T
RF04	887	1241	1	2	887	A	1241	T
RF04	1900	1920	4	2	1900	C	1920	A
RF05	41	499	2	2	41	T	499	A
RF05	499	879	1	2	499	A	879	C
RF05	795	1297	2	2	795	A	1297	T
RF05	879	1297	2	2	879	C	1297	T
RF06	517	543	9	2	517	C	543	G
RF06	668	695	4	2	668	G	695	T
RF07	225	684	1	2	225	C	684	G
RF07	225	684	1	2	225	C	684	T
RF08	926	992	2	2	926	C	992	G
RF09	294	317	6	2	294	T	317	A
RF10	139	175	1	2	139	T	175	G
RF10	139	175	3	2	139	T	175	C
```

END_DOC
*/
@Program(name="bam2haplotypes",
description="Reconstruct SNP haplotypes from reads",
keywords={"vcf","phased","genotypes","bam"},
biostars=9493599,
creationDate="20211015",
modificationDate="20211020",
menu="BAM Manipulation"
)
public class BamToHaplotypes extends MultiBamLauncher {
	private static final Logger LOG=Logger.of(BamToHaplotypes.class);

	
	/* what shall we do with ALT that are not in the VCF */
	private enum AltHandler {
		skip,
		warn,
		error,
		N,
		all
		}

	
	private static class Change implements Comparable<Change> {
		int pos1;
		byte ref;
		byte alt;
		
		@Override
		public int compareTo(Change o) {
			return Integer.compare(pos1, o.pos1);
			}
		
		@Override
		public boolean equals(Object obj)
			{
			if(obj==this) return true;
			final Change other= (Change)obj;
			return other.pos1 == this.pos1 && other.ref==this.ref && other.alt==this.alt;
			}
		@Override
		public int hashCode()
			{
			int h = Integer.hashCode(pos1);
			h = h*31 + Integer.hashCode(ref);
			h = h*31 + Integer.hashCode(alt);
			return h;
			}
		}

	private static class Haplotype implements Comparable<Haplotype> {
		final int tid;
		final List<Change> changes = new ArrayList<>();
		Haplotype(final int tid) {
			this.tid = tid;
			}
		@Override
		public int compareTo(final Haplotype o)
			{
			int i= Integer.compare(this.tid,o.tid);
			if(i!=0) return i;
			final int n= this.changes.size();
			i = Integer.compare(n, o.changes.size());
			if(i!=0) return i;
			for(int x=0;x<n;x++) {
				i = Integer.compare(this.changes.get(x).pos1, o.changes.get(x).pos1);
				if(i!=0) return i;
				}
			for(int x=0;x<n;x++) {
				i = Integer.compare(this.changes.get(x).ref, o.changes.get(x).ref);
				if(i!=0) return i;
				}
			for(int x=0;x<n;x++) {
				i = Integer.compare(this.changes.get(x).alt, o.changes.get(x).alt);
				if(i!=0) return i;
				}
			return 0;
			}

		public int getStart() {
			return changes.get(0).pos1;
			}

		public int getEnd() {
			return changes.get(this.changes.size()-1).pos1;
			}
		public boolean equals(final Object obj)
			{
			if(obj==this) return true;
			final Haplotype other= (Haplotype)obj;
			return this.tid==other.tid && this.changes.equals(other.changes);
			}
		@Override
		public int hashCode()
			{
			int h = Integer.hashCode(this.tid);
			h = h*31 + this.changes.hashCode();
			return h;
			}

		}
	
	private static class HaplotypeCodec extends AbstractDataCodec<Haplotype>{
		@Override
		public Haplotype decode(DataInputStream dis) throws IOException
			{
			int ctg;
			try {
				ctg = dis.readInt();
				}
			catch(final EOFException err) {
				return null;
				}
			final Haplotype h = new Haplotype(ctg);
			final int n= dis.readInt();
			for(int i=0;i< n;i++) {
				final Change c = new Change();
				c.pos1 = dis.readInt();
				c.ref = dis.readByte();
				c.alt = dis.readByte();
				h.changes.add(c);
				}
			return h;
			}
		@Override
		public void encode(DataOutputStream dos, final Haplotype h)
				throws IOException
			{
			dos.writeInt(h.tid);
			dos.writeInt(h.changes.size());
			for(Change c:h.changes) {
				dos.writeInt(c.pos1);
				dos.writeByte(c.ref);
				dos.writeByte(c.alt);
				}
			
			}
		@Override
		public AbstractDataCodec<Haplotype> clone()
			{
			return new HaplotypeCodec();
			}
	
	
	}

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-V","--vcf"},description="Indexed VCf file. Only diallelic SNP will be considered.",required=true)
	protected Path vcfFile=null;
	@Parameter(names={"--buffer-size"},description=BufferedVCFReader.OPT_BUFFER_DESC,splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int buffSizeInBp = 1_000;
	@Parameter(names={"--paired"},description="Activate Paired-end mode. Variant can be supported by the read or/and is mate. Input must be sorted on query name using for example 'samtools collate'.")
	private boolean paired_mode = false;
	@Parameter(names={"--ignore-discordant-rg"},description="In paired mode, ignore discordant read-groups RG-ID.")
	private boolean ignore_discordant_rg = false;
	@Parameter(names={"--alt"},description="How shall we handle ALT allele that are not in the VCF. skip, warn (skip and warning), error (raise an error), N (replace with 'N')), all: use all alleles.")
	private AltHandler altHandler = AltHandler.all;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	
	
	
	private	VCFReader vcfReader = null;
	private BufferedVCFReader bufferedVCFReader = null;
	
	private VariantContext simplify(final VariantContext vc) {
		if(!vc.isSNP()) return null;
		if(!vc.isBiallelic()) return null;
		
		return new VariantContextBuilder(vc).noID().passFilters().
				log10PError(VariantContext.NO_LOG10_PERROR).
				noGenotypes().
				attributes(Collections.emptyMap()).
				make();
		}
	
	@Override
	protected int beforeSam() {
		this.vcfReader = VCFReaderFactory.makeDefault().open(this.vcfFile, true);
		this.bufferedVCFReader = new BufferedVCFReader(this.vcfReader, this.buffSizeInBp);
		this.bufferedVCFReader.setSimplifier(V->simplify(V));
		return 0;
		}
	@Override
	protected void afterSam() {
		CloserUtil.close(this.bufferedVCFReader);
		CloserUtil.close(this.vcfReader);
		}
	

	@Override
	protected int processInput(
			final SAMFileHeader headerIn,
			final CloseableIterator<SAMRecord> iter0) {
		SortingCollection<Haplotype> sorting  =null;
		try {
		final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(headerIn);
		final String sample = headerIn.getReadGroups().
				stream().map(RG->RG.getSample()).
				filter(R->!StringUtils.isBlank(R)).
				findFirst().
				orElse("SAMPLE");
		
		sorting = SortingCollection.newInstance(
				Haplotype.class,
				new HaplotypeCodec(),
				(A,B)->A.compareTo(B),
				this.writingSortingCollection.getMaxRecordsInRam(),
				this.writingSortingCollection.getTmpPaths()
				);
		
		if(this.paired_mode) {
			try(EqualIterator<SAMRecord> iter = new EqualIterator<>(iter0,(A,B)->A.getReadName().compareTo(B.getReadName()))) {
			while(iter.hasNext()) {
					final LinkedList<SAMRecord> buffer = new LinkedList<>(iter.next());
					SAMRecord R1= null;
					SAMRecord R2= null;
					while(!buffer.isEmpty()) {
						final SAMRecord rec = buffer.pop();
						if(rec.getReadUnmappedFlag() || rec.isSecondaryOrSupplementary()){
							continue;
							}
						else if(!rec.getReadPairedFlag()) {
							scanVariants(dict,Collections.singletonList(rec),sorting);
							}
						else if(R1==null && rec.getFirstOfPairFlag()) {
							R1 = rec;
							}
						else if(R2==null && rec.getSecondOfPairFlag()) {
							R2 = rec;
							}
						else
							{
							continue;
							}
						}
					if(R1!=null && R2!=null) {
						if(R1.contigsMatch(R2)) {
							scanVariants(dict,Arrays.asList(R1,R2),sorting);
							}
						else
							{
							scanVariants(dict,Collections.singletonList(R1),sorting);
							scanVariants(dict,Collections.singletonList(R2),sorting);
							}
						}
					else if(R1!=null && R2==null) {
						scanVariants(dict,Collections.singletonList(R1),sorting);
						}
					else if(R2!=null && R1==null) {
						scanVariants(dict,Collections.singletonList(R2),sorting);
						}
					}
				}
			}
		else
			{
			while(iter0.hasNext()) {
				final SAMRecord rec = iter0.next();
				if(rec.getReadUnmappedFlag()){
					continue;
					}
				scanVariants(dict,Collections.singletonList(rec),sorting);
				}
			}
		sorting.doneAdding();
		sorting.setDestructiveIteration(true);
		try(CloseableIterator<Haplotype> iter = sorting.iterator()) {
			PeekableIterator<Haplotype> peek = new PeekableIterator<Haplotype>(iter);
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				out.println("#CHROM\tSTART\tEND\tSAMPLE\tN-HAPLOTYPES\tN-VARIANTS\t(POS\\tALT)+");
				while(peek.hasNext()) {
					int n=1;
					final Haplotype hap = peek.next();
					while(peek.hasNext()) {
						final Haplotype hap2 = peek.peek();
						if(!hap.equals(hap2)) break;
						peek.next();//consumme
						n++;
						}
					out.print(dict.getSequence(hap.tid).getContig());
					out.print("\t");
					out.print(hap.getStart());
					out.print("\t");
					out.print(hap.getEnd());
					out.print("\t");
					out.print(sample);
					out.print("\t");
					out.print(n);
					out.print("\t");
					out.print(hap.changes.size());
					for(Change c:hap.changes) {
						out.print("\t");
						out.print(c.pos1);
						out.print("\t");
						out.print((char)c.alt);
						}
					out.println();

					}
				out.flush();
				}
			peek.close();
			}
		sorting.cleanup();
		return 0;
		}
	catch(Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

	
	private byte alleleToChar(final Allele a) {
		if(a==null || a.isSymbolic() || a.isNoCall() || a.length()!=1) throw new IllegalArgumentException("length!=1 for allele "+a);
		return a.getBases()[0];
	}
	
	private void scanVariants(
			final SAMSequenceDictionary dict,
			final List<SAMRecord> buffer,
			final SortingCollection<Haplotype> sorting
			) {
			if(buffer.isEmpty()) return;
			
			final Set<Change> change_to_test = new HashSet<>();

			for(final SAMRecord rec:buffer) {
				try(CloseableIterator<VariantContext> iter=this.bufferedVCFReader.query(rec)) {
				while(iter.hasNext()) {
					final VariantContext ctx = iter.next();
					if(!ctx.isBiallelic() || !ctx.isSNP()) continue;
					final Change c = new Change();
					c.pos1  = ctx.getStart();
					c.ref = alleleToChar(ctx.getReference());
					c.alt = alleleToChar(ctx.getAlleles().get(1));
					if(c.ref==c.alt) throw new IllegalArgumentException("REF==ALT in "+ctx);
					change_to_test.add(c);
					}
				}
			}
			if(change_to_test.size()<2) return;
			
			final SAMSequenceRecord ssr = dict.getSequence(buffer.get(0).getContig());
			if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(buffer.get(0).getContig(), dict);
			final Haplotype hap = new Haplotype(ssr.getSequenceIndex());

			
			
			for(final SAMRecord rec:buffer) {
				final byte[] bases = rec.getReadBases();
				if(bases==null || bases==SAMRecord.NULL_QUALS || bases.length==0) {
					continue;
					}

				for(AlignmentBlock ab:rec.getAlignmentBlocks()) {
					final int readPos1 = ab.getReadStart();
					final int refPos1 = ab.getReferenceStart();
					for(int x=0;x< ab.getLength();++x) {
						final int refPos1_x = refPos1+x;
						final Change change =  change_to_test.stream().filter(C->C.pos1==refPos1_x).findFirst().orElse(null);
						if(change==null) continue;
						//paired overlapping discordant reads, already inserted in haplotype
						if(hap.changes.stream().anyMatch(C->C.pos1==change.pos1)) continue;
						byte readBase = bases [ (readPos1-1) + x ];
						final boolean allele_is_known = (change.ref==readBase || change.alt==readBase);
						if(!allele_is_known) {
							final String msg = rec.getReadName()+" at "+change.pos1+" is "+(char)readBase+" while vcf is "+(char)change.ref+"/"+(char)change.alt;
							switch(this.altHandler) {
								case error: throw new RuntimeIOException(msg);
								case warn: LOG.warn(msg);continue;
								case skip: continue;
								case N:  readBase = (byte)'N'; break;
								case all: break;
								default: throw new IllegalStateException();
								}
							}

						final Change c = new Change();
						c.pos1 = change.pos1;
						c.ref = change.ref;
						c.alt = readBase;
						hap.changes.add(c);
						}
					}
				}
			
			if(hap.changes.size()<2) {
				return;
				}
			Collections.sort(hap.changes,(A,B)->A.compareTo(B));
			sorting.add(hap);
			}
	
	
	public static void main(final String[] args) {
		new BamToHaplotypes().instanceMainWithExit(args);

	}

}
