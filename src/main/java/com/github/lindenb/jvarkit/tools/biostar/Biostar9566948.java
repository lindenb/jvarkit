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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;

/**

BEGIN_DOC

## Example

```bash
$ java -jar dist/jvarkit.jar biostar9566948 -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam   -S 
@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S1	SM:S1	LB:L1	CN:Nantes
@CO	biostar9566948. compilation:20230621100700 githash:8a9a881b htsjdk:3.0.4 date:20230621100755. cmd:-R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam -S
RF01_1_483_2:0:0_3:0:0_41	163	RF01	1	60	1=	=	414	483	G	2	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:60	XS:i:0
RF01_8_542_1:0:0_2:0:0_95	99	RF01	8	60	1=	=	473	535	A	2	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:69	XS:i:0
RF01_11_507_0:0:0_1:0:0_9e	99	RF01	11	60	1=	=	438	497	G	2	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:70	XS:i:0
RF01_12_501_0:0:0_2:0:0_62	99	RF01	12	60	1=	=	432	490	C	2	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:70	XS:i:0
RF01_27_590_3:0:0_1:0:0_68	163	RF01	27	60	1=	=	521	564	G	2	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:55	XS:i:0
RF01_44_622_1:0:0_1:0:0_3a	99	RF01	44	60	1=	=	553	579	C	2	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:65	XS:i:0
```

```
$ java -jar dist/jvarkit.jar biostar9566948 -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam 
@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S1	SM:S1	LB:L1	CN:Nantes
@CO	biostar9566948. compilation:20230621100700 githash:8a9a881b htsjdk:3.0.4 date:20230621100859. cmd:-R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam
RF01_1_483_2:0:0_3:0:0_41	163	RF01	1	60	1=69S	=	414	483	GGCTATTAAAGCTATACAATGGGGCCGTATAATCTAATCTTGTCAGAATATTTATCATTTATATATAACT	2222222222222222222222222222222222222222222222222222222222222222222222	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:60	XS:i:0
RF01_8_542_1:0:0_2:0:0_95	99	RF01	8	60	1=69S	=	473	535	AAAGCTATACAATGGGGAAGTATAATCTAATCTTGTCAGAATATTTATCATTTATATATAACTCACAATG	2222222222222222222222222222222222222222222222222222222222222222222222	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:69	XS:i:0
RF01_11_507_0:0:0_1:0:0_9e	99	RF01	11	60	1=69S	=	438	497	GCTATACAATGGGGAAGTATAATCTAATCTTGTCAGAATATTTATCATTTATATATAACTCACAATCCGC	2222222222222222222222222222222222222222222222222222222222222222222222	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:70	XS:i:0
RF01_12_501_0:0:0_2:0:0_62	99	RF01	12	60	1=69S	=	432	490	CTATACAATGGGGAAGTATAATCTAATCTTGTCAGAATATTTATCATTTATATATAACTCACAATCCGCA	2222222222222222222222222222222222222222222222222222222222222222222222	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:70	XS:i:0
RF01_27_590_3:0:0_1:0:0_68	163	RF01	27	60	1=69S	=	521	564	GTATCATCTAATCTTGTCATAATATTTATCATATATATATAACTCACAATCCGCAGTTCAAATTCCAATA	2222222222222222222222222222222222222222222222222222222222222222222222	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:55	XS:i:0
RF01_44_622_1:0:0_1:0:0_3a	99	RF01	44	60	1=69S	=	553	579	CAGAATATTTATCATTTATATATAACTCAGAATCCGCAGTTCAAATTCCAATATACTATTCTTCCAATAG	2222222222222222222222222222222222222222222222222222222222222222222222	MC:Z:M1	RG:Z:S1	NM:i:0	AS:i:65	XS:i:0
```

END_DOC

*/
@Program(name="biostar9566948",
	description="Trim Reads So Only First Base Remains",
	biostars={9566948},
	keywords={"bam","sam"},
	modificationDate="20230621",
	creationDate="20230621",
	jvarkit_amalgamion =  true,
	menu="BAM Manipulation"
	)
public class Biostar9566948 extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(Biostar9566948.class).make();
	private ReferenceSequenceFile referenceSequenceFile = null;
	private GenomicSequence genomicSequence =null;
	@Parameter(names={"--disable-soft-clipping","-S"},description="disable soft clipping. Remove bases.")
	private boolean no_soft_clipping = false;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam()
		{
		this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(getRequiredReferencePath());
		return super.beforeSam();
		}
	@Override
	protected void afterSam() {
		try { this.referenceSequenceFile.close();}
		catch(IOException err) {}
		super.afterSam();
		}
	
	@Override
	protected void scanIterator(
			final SAMFileHeader headerIn,
			CloseableIterator<SAMRecord> iter,
			final SAMFileWriter out)
		{
		while(iter.hasNext()) {
			SAMRecord rec = iter.next();
			if(rec.getReadUnmappedFlag()) {
				out.addAlignment(rec);
				continue;
				}
			final byte[] bases = rec.getReadBases();
			final byte[] quals  = rec.getBaseQualities();
			if(bases==SAMRecord.NULL_SEQUENCE) {
				out.addAlignment(rec);
				continue;
				}
			final int refpos1 = rec.getAlignmentStart();
			int readpos1 = SAMRecord.getReadPositionAtReferencePosition(rec, refpos1 , false);
			if(readpos1<=0) throw new IllegalStateException("cannot get read pos 1 in " + rec.getSAMString());
			if(genomicSequence==null || !genomicSequence.getContig().equals(rec.getContig())) {
				genomicSequence = new GenomicSequence(this.referenceSequenceFile, rec.getContig());
				}
			
			
			final char base = refpos1>0 && refpos1<= this.genomicSequence.length() ?
					Character.toUpperCase(this.genomicSequence.charAt(refpos1-1)):
					'N';
			
			final CigarElement cigarElement;
			if(bases[readpos1-1]==base) {
				cigarElement = new CigarElement(1,CigarOperator.EQ);
				rec.setAttribute(SAMTag.NM,0);
			} else {
				cigarElement = new CigarElement(1,CigarOperator.X);
				rec.setAttribute(SAMTag.NM,1);
			}
			
			final int readpos0 = readpos1 -1;
			List<CigarElement> cigarlist = new ArrayList<>(3);
			if(readpos0>0 && !no_soft_clipping) 	cigarlist.add(new CigarElement(readpos0,CigarOperator.S));
			cigarlist.add(cigarElement);
			if(readpos0+1 <= bases.length && !no_soft_clipping) {
				cigarlist.add(new CigarElement(bases.length - (readpos0+1),CigarOperator.S));
				}
			
			rec.setCigar(new Cigar(cigarlist));
			if(quals!=SAMRecord.NULL_QUALS && no_soft_clipping) {
				rec.setBaseQualities(new byte[] {quals[readpos0]});
				}
			
			if(no_soft_clipping) {
				rec.setReadBases(new byte[] {bases[readpos0]});
				}
			
			if(rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
				rec.setAttribute(SAMTag.MC, "M1");
				}
			out.addAlignment(rec);
			}
		}
	


	public static void main(final String[] args) throws IOException
		{
		new Biostar9566948().instanceMainWithExit(args);
		}
	}
