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

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
 * 
 BEGIN_DOC

## Description

input must be sorted on read name using `samtools sort -n` or ` samtools collate`

it will be **much** faster is the reads belong to one chromosome

## Example

```
$ samtools view -O BAM --reference "ref.fasta" in.cram "chr22:41201525-41490147" |\
	samtools collate -O -u - |\
	java -jar dist/biostar489074.jar --reference "ref.fasta"
```

 END_DOC
 *
 */

@Program(name="biostar489074",
keywords={"sam","bam","vcf","call"},
description="call variants for every paired overlaping read",
biostars= {489074},
creationDate="20200205",
modificationDate="20210412",
generate_doc=true
)
public class Biostar489074 extends MultiBamLauncher {			
private static final Logger LOG = Logger.build(Biostar489074.class).make();
@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private Path outputFile = null;
@Parameter(names={"--groupby","--partition"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
private SAMRecordPartition groupBy=SAMRecordPartition.sample;
@Parameter(names= {"--ploidy"},description="default ploidy")
private int ploidy=2;

@ParametersDelegate
private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
@ParametersDelegate
private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();


private ReferenceSequenceFile indexedFastaRef=null;

private static class Call {
	int sampleid;
	int tid;
	int pos;
	byte ref;
	byte alt;
	@Override
	public int hashCode()
		{
		final int prime = 31;
		int result = 1;
		result = prime * result + alt;
		result = prime * result + pos;
		result = prime * result + ref;
		result = prime * result + sampleid;
		result = prime * result + tid;
		return result;
		}
	@Override
	public boolean equals(Object obj)
		{
		if (this == obj) return true;
		if (obj == null) return false;
		if (!(obj instanceof Call)) return false;
		final Call other = (Call) obj;
		if (alt != other.alt) return false;
		if (pos != other.pos) return false;
		if (ref != other.ref) return false;
		if (sampleid != other.sampleid) return false;
		if (tid != other.tid) return false;
		return true;
		}
	int compare1(Call c) {
		int i= Integer.compare(this.tid, c.tid);
		if(i!=0) return i;
		i= Integer.compare(this.pos, c.pos);
		if(i!=0) return i;
		i= Byte.compare(this.ref, c.ref);
		return i;
		}
	int compare2(Call c) {
		int i= compare1(c);
		if(i!=0) return i;
		i= Integer.compare(this.sampleid, c.sampleid);
		if(i!=0) return i;
		i= Byte.compare(this.alt, c.alt);
		return i;
		}

	}

private static class CallCodec extends AbstractDataCodec<Call> {
	@Override
	public Call decode(DataInputStream dis) throws IOException
		{
		Call c= new Call();
		try {
			c.tid = dis.readInt();
			}
		catch(IOException err) {
			throw new EOFException();
		}
		c.pos = dis.readInt();
		c.sampleid = dis.readInt();
		c.ref = dis.readByte();
		c.alt = dis.readByte();
		return c;
		}
	@Override
	public void encode(DataOutputStream dos, Call c) throws IOException
		{
		dos.writeInt(c.tid);
		dos.writeInt(c.pos);
		dos.writeInt(c.sampleid);
		dos.writeByte(c.ref);
		dos.writeByte(c.alt);
		}
	@Override
	public CallCodec clone()
		{
		return new CallCodec();
		}
	}

@Override
protected Logger getLogger() {
	return LOG;
	}

@Override
protected int processInput(final SAMFileHeader header, final CloseableIterator<SAMRecord> iter0) {
	VariantContextWriter out=null;
	GenomicSequence genome = null;
	SortingCollection<Call> sorting  = null;
	try {
		this.indexedFastaRef = ReferenceSequenceFileFactory.getReferenceSequenceFile(getRequiredReferencePath());
		if(!(header.getSortOrder().equals(SAMFileHeader.SortOrder.unsorted) || header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname))) {
			getLogger().error("input should be sorted with 'samtools sort -n' or 'samtools collate' but got " + header.getSortOrder());
			return -1;
			}
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
		this.writingVariantsDelegate.dictionary(dict);
		
		sorting =  SortingCollection.newInstance(Call.class,
				new CallCodec(),
				(A,B)->A.compare2(B),
				this.writingSortingCollection.getMaxRecordsInRam(),
				this.writingSortingCollection.getTmpPaths()
				);
		sorting.setDestructiveIteration(true);
		
		final List<String> samples = header.getReadGroups().stream().
			map(R->this.groupBy.apply(R)).
			filter(S->!StringUtils.isBlank(S)).
			sorted().
			collect(Collectors.toSet()).
			stream().
			collect(Collectors.toList());
		
		final Map<String,Integer> rgid2idx = new HashMap<>();
		for(final SAMReadGroupRecord rg : header.getReadGroups()) {
			final String sn = this.groupBy.apply(rg);
			if(StringUtils.isBlank(sn)) continue;
			int idx = samples.indexOf(sn);
			if(idx==-1) continue;
			rgid2idx.put(rg.getId(),idx);
			}

			
		try(CloseableIterator<List<SAMRecord>> iter = new EqualIterator<>(
				iter0,
				(A,B)->A.getReadName().compareTo(B.getReadName()))
					) {
			while(iter.hasNext()) {
				final List<SAMRecord> buffer = iter.next();
				int read1_idx = -1;
				int read2_idx = -1;
				for(int i=0;i< buffer.size();i++) {
					final SAMRecord rec = buffer.get(i);
					if(!rec.getReadPairedFlag()) continue;
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getMateUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getReadBases()==SAMRecord.NULL_SEQUENCE) continue;
					if(rec.getFirstOfPairFlag()) {
						read1_idx = i;
					}
					else if(rec.getSecondOfPairFlag()) {
						read2_idx = i;
					}
				}
				
				
				if(read1_idx==-1 || read2_idx==-1 || read1_idx==read2_idx)  continue;
				final SAMRecord rec1a = buffer.get(read1_idx);
				final SAMRecord rec2a = buffer.get(read2_idx);
				if(!rec1a.overlaps(rec2a)) continue;
				final int chromStart = Math.max(rec1a.getStart(),rec2a.getStart());
				final int chromEnd  = Math.min(rec1a.getEnd(),rec2a.getEnd());
				if(chromStart > chromEnd)  continue;
				final Integer sampleid = rgid2idx.get(rec1a.getReadGroup().getId());
				if(sampleid==null) continue;
				
				if(genome==null || !genome.contigsMatch(rec1a)) {
					genome = new GenomicSequence(this.indexedFastaRef, rec1a.getContig());
				}
				final Set<Call> calls1 = new HashSet<>();
				final Set<Call> calls2 = new HashSet<>();
				
				for(int side=0;side<2;side++) {
					final SAMRecord rec = (side==0?rec1a:rec2a);
					final Set<Call> calls = (side==0?calls1:calls2);
					final byte[] bases = rec.getReadBases();
					for(AlignmentBlock ab:rec.getAlignmentBlocks()) {
						for(int n=0;n< ab.getLength();++n) {
							final int ref= ab.getReferenceStart() + n;
							if(ref <chromStart) continue;
							if(ref >chromEnd) break;
							if(ref<0 || ref>= genome.length()) continue;
							final byte refBase = (byte)Character.toUpperCase(genome.charAt(ref-1));//0 based
							if(!AcidNucleics.isATGC(refBase)) continue;
							final byte readBase = (byte)Character.toUpperCase(bases[(ab.getReadStart()-1/* 1 based */)+n]);//0 based
							if(readBase==refBase) continue;
							final Call call = new Call();
							call.sampleid = sampleid;
							call.tid = rec1a.getReferenceIndex().intValue();
							call.ref= refBase;
							call.alt = readBase;
							call.pos = ref;
							calls.add(call);
							}
						}
					}
				calls1.retainAll(calls2);
				if(calls1.isEmpty()) continue;
				for(final Call c:calls1) {
					sorting.add(c);
		
					}
				}
			sorting.doneAdding();
			out = this.writingVariantsDelegate.dictionary(dict).open(this.outputFile);
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFStandardHeaderLines.addStandardInfoLines(metaData,true,VCFConstants.DEPTH_KEY);
			VCFStandardHeaderLines.addStandardInfoLines(metaData,true,VCFConstants.ALLELE_COUNT_KEY);
			VCFStandardHeaderLines.addStandardInfoLines(metaData,true,VCFConstants.ALLELE_NUMBER_KEY);
			VCFStandardHeaderLines.addStandardInfoLines(metaData,true,VCFConstants.ALLELE_FREQUENCY_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData,true,VCFConstants.GENOTYPE_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData,true,VCFConstants.DEPTH_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData,true,VCFConstants.GENOTYPE_ALLELE_DEPTHS);
			final VCFHeader header2 = new VCFHeader(metaData, samples);
			header2.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header2);
			out.writeHeader(header2);
			try(CloseableIterator<Call> iter1 = sorting.iterator()) {
				try(EqualRangeIterator<Call> iter2 = new EqualRangeIterator<>(iter1, (A,B)->A.compare1(B))) {
					while(iter2.hasNext()) {
						final List<Call> calls = iter2.next();
						if(calls.isEmpty()) continue;
						final Call first = calls.get(0);
						final Set<Allele> altAllelesSet = calls.stream().
							map(A->Allele.create(A.alt,false)).
							collect(Collectors.toSet());
						final List<Allele> altAllelesList = new ArrayList<>(altAllelesSet);
						final List<Allele> vcAlleles = new ArrayList<>(altAllelesList.size()+1);
						vcAlleles.add(Allele.create(first.ref,true));
						vcAlleles.addAll(altAllelesList);
	
						
						final List<Genotype> genotypes = new ArrayList<>(samples.size());
						int DP=0;
						
						for(int i=0;i<samples.size();i++) {
							final String sn = samples.get(i);
							final Counter<Allele> allele2count = new Counter<>();
							for(Call c:calls) {
								if(c.sampleid!=i) continue;
								allele2count.incr(Allele.create(c.alt,false));
								}
							Genotype gt;
							if(allele2count.isEmpty()) {
								gt=GenotypeBuilder.createMissing(sn,this.ploidy);
								}
							else
								{
								final int[] array = new int[vcAlleles.size()];
								for(int j=0;j< vcAlleles.size();j++) {
									array[j] = (int)allele2count.count(vcAlleles.get(j));
									}
								final GenotypeBuilder gb=new GenotypeBuilder(sn, new ArrayList<>(allele2count.keySet()));
								gb.AD(array);
								gb.DP((int)allele2count.getTotal());
								DP+=(int)allele2count.getTotal();
								gt= gb.make();
								}
							genotypes.add(gt);
							}
						
						final VariantContextBuilder vcb = new VariantContextBuilder(null,
								dict.getSequence(first.tid).getContig(),
								first.pos, first.pos,
								vcAlleles
								);
						vcb.attribute(VCFConstants.DEPTH_KEY, DP);
						vcb.genotypes(genotypes);
						
						out.add(vcb.make());
						}
					}
				}
			}
		sorting.cleanup();
		out.close();
		out=null;
		out=null;
		this.indexedFastaRef.close();
		this.indexedFastaRef = null;
		return 0;
	} catch(final Throwable err) {
		getLogger().error(err);
		return -1;
	} finally {
		CloserUtil.close(out);
		CloserUtil.close(this.indexedFastaRef);;
	}
}

public static void main(final String[] args) {
	new Biostar489074().instanceMainWithExit(args);
}
}
