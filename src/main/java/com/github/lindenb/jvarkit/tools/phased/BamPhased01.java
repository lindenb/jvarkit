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
package com.github.lindenb.jvarkit.tools.phased;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

END_DOC
*/
@Program(name="bamphased01",
description="Extract Reads supporting at least two variants.",
keywords={"vcf","phased","genotypes","svg","bam"},
creationDate="20210218",
modificationDate="20210218"
)
public class BamPhased01 extends OnePassBamLauncher {
	private static final Logger LOG=Logger.build(BamPhased01.class).make();

	
	private static class PosToCheck extends SimplePosition {
		final Set<Byte> alts;
		PosToCheck(final String ctg,int pos,final Set<Byte> alts) {
			super(ctg,pos);
			this.alts=alts;
		}
		@Override
		public String toString() {
			return super.toString()+" "+ alts;
			}
		}

	
	@Parameter(names={"-V","--vcf"},description="Indexed VCf file",required=true)
	protected Path vcfFile=null;
	@Parameter(names={"--buffer-size"},description=BufferedVCFReader.OPT_BUFFER_DESC,splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int buffSizeInBp = 1_000;
	@Parameter(names={"--tag"},description="Tag in metadata of read containing the position. Ignored if empty")
	private String XTAG = "XP";

	
	private	VCFReader vcfReader = null;
	private BufferedVCFReader bufferedVCFReader = null;
	private Set<String> samplesInBam = new HashSet<>();
	
	private VariantContext simplify(VariantContext vc) {
		if(!vc.isSNP()) return null;
		if(!vc.isBiallelic()) return null;
		return vc;
	}
	
	@Override
	protected int beforeSam() {
		if(!(this.XTAG.length()==0 || this.XTAG.length()==2)) {
			LOG.error("tag should be empty of length==2 but got "+this.XTAG);
			return -1;
			}
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
	protected SAMFileHeader createOutputHeader(final SAMFileHeader headerIn) {
		final SAMFileHeader outHeader= super.createOutputHeader(headerIn);
		this.samplesInBam = headerIn.getReadGroups().
				stream().
				map(RG->RG.getSample()).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toSet());
		
		this.samplesInBam.retainAll(this.vcfReader.getHeader().getSampleNamesInOrder());
		if(this.samplesInBam.isEmpty()) {
			throw new RuntimeIOException("No overlapping samples between the SAM read groups (@RG SM:xxxx) and the vcf file "+this.vcfFile);
			}
		
		return outHeader;
		}
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction() {
		return rec->{
			if(rec.getReadUnmappedFlag()) return Collections.emptyList();
			final SAMReadGroupRecord rg = rec.getReadGroup();
			if(rg==null) return Collections.emptyList();
			final byte[] bases = rec.getReadBases();
			if(bases==null || bases==SAMRecord.NULL_QUALS || bases.length==0) Collections.emptyList();
			final String sn = rg.getSample();
			if(StringUtils.isBlank(sn) || !this.samplesInBam.contains(sn)) return Collections.emptyList();
			
			
			final List<PosToCheck> candidates = new ArrayList<>();
			try(CloseableIterator<VariantContext> iter=this.bufferedVCFReader.query(rec)) {
				while(iter.hasNext()) {
					final VariantContext ctx = iter.next();
					final Genotype gt = ctx.getGenotype(sn);
					if(gt.isHomRef() || gt.isNoCall()) continue;
					final Set<Byte> alts = gt.getAlleles().stream().
						filter(A->A.isCalled() && !A.isReference() && !A.isSymbolic()).
						filter(A->AcidNucleics.isATGC(A)).
						map(A->(byte)A.getDisplayString().charAt(0)).
						collect(Collectors.toSet());
					if(alts.isEmpty()) continue;
					final PosToCheck pos = new PosToCheck(ctx.getContig(),ctx.getStart(),alts);
					candidates.add(pos);
					}
				}
			if(candidates.size()<2) return Collections.emptyList();
			
			final List<PosToCheck> supporting = new ArrayList<>(candidates.size());

			for(AlignmentBlock ab:rec.getAlignmentBlocks()) {
				final int readPos1 = ab.getReadStart();
				final int refPos1 = ab.getReferenceStart();
				for(int x=0;x< ab.getLength();++x) {
					for(PosToCheck pos:candidates) {
						if(pos.getPosition() != refPos1+x) continue;
						final byte readBase = bases [ (readPos1-1) + x ];
						if(pos.alts.contains(readBase)) {
							supporting.add(pos);
							break;
							}
					}
				}
			}
			
			if(supporting.size()<2) return Collections.emptyList();
			
			if(!StringUtils.isBlank(this.XTAG)) {
				rec.setAttribute(this.XTAG,
						supporting.stream().
							map(S->S.getContig()+"_"+S.getStart()+"_"+S.alts.stream().map(B->""+(char)B.byteValue()).collect(Collectors.joining("_"))).
							collect(Collectors.joining(";"))
						);
					}

			
			return Collections.singletonList(rec);
			};
		
		}
	
	
	public static void main(final String[] args) {
		new BamPhased01().instanceMainWithExit(args);

	}

}
