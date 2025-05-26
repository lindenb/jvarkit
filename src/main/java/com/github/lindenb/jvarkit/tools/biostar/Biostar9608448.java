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
package com.github.lindenb.jvarkit.tools.biostar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord.PlatformValue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CoordMath;


/**
BEGIN_DOC

Motivation:

convert single end long read to paired end reads.

Note that short single read, paired reads, unmapped read, secondary & supplementary alignements are discarded

Example:

```
$ java -jar dist/jvarkit.jar  biostar9608448 src/test/resources/FAB23716.nanopore.bam -L 10 | samtools view | head
44767a9a-a0b9-4d7e-a324-d0d3ea113d8c_Basecall_Alignment_template	67	chr1	17123	38	6M2D4M	=	31047	13925	GTGCGCCGCT	-2.1.,(')-	MC:Z:10M
44767a9a-a0b9-4d7e-a324-d0d3ea113d8c_Basecall_Alignment_template	147	chr1	31047	38	10M	=	17123	-13925	CACCTTGAAC	&#&$$$%$$%	MC:Z:6M2D4M
d324a4bc-aa2c-4ee8-be69-934cc58c0003_Basecall_Alignment_template	67	chr1	38469	1	10M	=	43735	5267	ATGCTGCCTG	2-,.314443	MC:Z:10M
d324a4bc-aa2c-4ee8-be69-934cc58c0003_Basecall_Alignment_template	147	chr1	43735	1	10M	=	38469	-5267	AGCAAACTTT	-',12()./.	MC:Z:10M
76862e2e-98eb-4ad3-a523-6a8709c0b56a_Basecall_Alignment_template	67	chr1	44403	0	10M	=	44481	79	TCAACAACAA	&&&%&)'''&	MC:Z:10M
76862e2e-98eb-4ad3-a523-6a8709c0b56a_Basecall_Alignment_template	147	chr1	44481	0	10M	=	44403	-79	GGTAGCCGAA	''&$%(((%*	MC:Z:10M
3330d9a6-d2a9-423b-accc-92a6d1fe646e_Basecall_Alignment_template	67	chr1	52105	1	8M2D2M	=	53738	1634	ATTCCTACGA	%).,.%+$))	MC:Z:10M
3330d9a6-d2a9-423b-accc-92a6d1fe646e_Basecall_Alignment_template	147	chr1	53738	1	10M	=	52105	-1634	ACTTAGGCAA	,)((%%''((	MC:Z:8M2D2M
c6055e6a-9a1c-4126-84ec-64549fd4d264_Basecall_Alignment_template	67	chr1	63945	5	10M	=	67887	3943	TCACCATGAT	*+'*-.111-	MC:Z:10M
c6055e6a-9a1c-4126-84ec-64549fd4d264_Basecall_Alignment_template	147	chr1	67887	5	10M	=	63945	-3943	AGTATTATCA	+$+*+/((*&	MC:Z:10M
```

END_DOC

 */

@Program(name="biostar9608448",
	keywords={"sam","bam","paired","long","illumina","nanopore"},
	description="Convert long reads to short paired reads",
	creationDate = "20250130",
	modificationDate = "20250130",
	biostars=9608448,
	jvarkit_amalgamion = true
	)
public class Biostar9608448 extends OnePassBamLauncher
	{
	private static final Logger LOG = Logger.of(Biostar9608448.class);
	@Parameter(names={"--read-length","-L"},description="short read length")
	private int readLength = 150;
	
	private SAMFileHeader outHeader;
	
	private static class Base {
		int ref;
		CigarOperator op;
		byte base;
		byte qual;
		@Override
		public String toString() {
			return "ref:"+ref+" op:"+op+"\n";
			}
	}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}

	@Override
	protected SAMFileHeader createOutputHeader(SAMFileHeader headerIn) {
		this.outHeader = super.createOutputHeader(headerIn);
		this.outHeader .setSortOrder(SAMFileHeader.SortOrder.unsorted);
		this.outHeader .setReadGroups(headerIn.getReadGroups().stream().
				map(RG->{
					RG.setPlatform(PlatformValue.ILLUMINA.name());
					return RG;
				}).collect(Collectors.toList()));
		return this.outHeader;
		}
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction() {
	return READ->{
		if(READ.isSecondaryOrSupplementary()) return Collections.emptyList();
		if(READ.getReadUnmappedFlag()) return Collections.emptyList();
		if(READ.getReadPairedFlag()) return Collections.singletonList(READ);
		if(READ.getReadLength()< this.readLength*2) return Collections.emptyList();
		if(READ.getReadBases()==SAMRecord.NULL_SEQUENCE) return Collections.emptyList();
		if(READ.getBaseQualities()==SAMRecord.NULL_QUALS) return Collections.emptyList();
		final Cigar cigar = READ.getCigar();
		if(cigar.isEmpty()) return Collections.emptyList();

		final List<Base>  array = new ArrayList<>(READ.getReadLength());
		byte[] bases = READ.getReadBases();
		byte[] quals = READ.getBaseQualities();

		int refpos=READ.getUnclippedStart();
		int readpos=0;
		for(CigarElement ce: cigar.getCigarElements()) {
			final CigarOperator op = ce.getOperator();
			if(op.equals(CigarOperator.P)) continue;
 			for(int x=0;x< ce.getLength();++x) {
				final Base base = new Base();
				base.op = op;
				base.ref = refpos;
				switch(op) {
					case D: case N:case H:
						refpos++;
						break;
					case I:
						base.base = bases[readpos];
						base.qual = quals[readpos];
						readpos++;
						break;
					case S: case M: case X: case EQ:
						base.base = bases[readpos];
						base.qual = quals[readpos];
						readpos++;
						refpos++;
						break;
					default:
						throw new IllegalStateException();
					}
				array.add(base);
				}
			}
		
		// remove clipping/non-align in 5' en 3'
		while(!array.isEmpty()) {
			if(!array.get(array.size()-1).op.isAlignment())
				{
				array.remove(array.size()-1);
				continue;
				}
			if(!array.get(0).op.isAlignment())
				{
				array.remove(0);
				continue;
				}
			break;
			}
		if(array.isEmpty()) return Collections.emptyList();
		
		final SAMRecord[] reads=new SAMRecord[2];
		for(int side=0;side<2;++side) {
			final SAMRecord rec = READ.deepCopy();
			
			rec.setHeader(outHeader);
			rec.setReadPairedFlag(true);
			rec.setProperPairFlag(true);
			rec.setMateUnmappedFlag(false);
			rec.setMappingQuality(60);
			rec.clearAttributes();
			
			final List<Base> L2=new ArrayList<>(this.readLength);
			int x = (side==0?0:array.size()-1);
			int new_read_len=0;
			
			while( new_read_len < this.readLength && x>=0 && x< array.size()) {
				final Base b = array.get(x);
				L2.add(b);
				if(b.op.isAlignment() || b.op.equals(CigarOperator.I)) {
					new_read_len++;
					}
				if(b.op.isClipping()) throw new IllegalStateException();
				x+=(side==0?1:-1);
				}
			//for 3' side , reverse the order
			if(side==1) {
				Collections.reverse(L2);
				}
			
			
			//cleanup 3' side of reads
			x=L2.size()-1;
			while(x>=0) {
				Base b = L2.get(x);
				if(b.op.isAlignment())break;
				b.op=CigarOperator.S;
				--x;
				}
			if(x==0) throw new IllegalStateException();

			rec.setCigar(Cigar.fromCigarOperators(L2.stream().map(B->B.op).collect(Collectors.toList())));
			rec.setReadString(L2.stream().filter(B->!(B.op.equals(CigarOperator.D) || B.op.equals(CigarOperator.N))).map(B->""+(char)B.base).collect(Collectors.joining()));
			rec.setBaseQualityString(L2.stream().filter(B->!(B.op.equals(CigarOperator.D) || B.op.equals(CigarOperator.N))).map(B->""+SAMUtils.phredToFastq(B.qual)).collect(Collectors.joining()));
			rec.setAlignmentStart(L2.get(0).ref);
			
			final Object rg= READ.getAttribute(SAMTag.RG);
			if(rg!=null) rec.setAttribute(SAMTag.RG,rg);

			
			switch(side) {
				case 0: 
					rec.setFirstOfPairFlag(true);
					rec.setSecondOfPairFlag(false);
					rec.setReadNegativeStrandFlag(false);
					break;
				case 1: 
					rec.setFirstOfPairFlag(false);
					rec.setSecondOfPairFlag(true);
					rec.setReadNegativeStrandFlag(true);
					break;
				}
			//restore order
			
			reads[side]=rec;
			}
		
		final int distance = CoordMath.getLength(reads[0].getStart(),reads[1].getStart());
		reads[0].setInferredInsertSize(distance);
		reads[1].setInferredInsertSize(-distance);
		
		reads[0].setMateReferenceIndex(reads[1].getReferenceIndex());
		reads[0].setMateAlignmentStart(reads[1].getAlignmentStart());
		reads[1].setMateReferenceIndex(reads[0].getReferenceIndex());
		reads[1].setMateAlignmentStart(reads[0].getAlignmentStart());
		reads[0].setAttribute(SAMTag.MC, reads[1].getCigarString());
		reads[1].setAttribute(SAMTag.MC, reads[0].getCigarString());
		
		
		return Arrays.asList(reads[0],reads[1]);
		};
	}
	
	public static void main(final String[] args) {
		new Biostar9608448().instanceMainWithExit(args);
	}
}
