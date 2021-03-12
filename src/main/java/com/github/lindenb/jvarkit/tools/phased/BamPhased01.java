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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

for @isamtalves

Only diallelic SNV are supported.

## Example

```
$ java -jar dist/bamphased01.jar \
	-V src/test/resources/rotavirus_rf.vcf.gz \
	src/test/resources/S4.bam 


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
@RG	ID:S4	SM:S4	LB:L4	CN:Nantes
@PG	ID:0	CL:-V src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S4.bam	VN:04c54fe	PN:bamphased01
@CO	bamphased01. compilation:20210218183501 githash:04c54fe htsjdk:2.23.0 date:20210218183539. cmd:-V src/test/resources/rotavirus_rf.vcf.gz src/test/resources/S4.bam
RF10_109_650_2:2:0_2:0:0_13	99	RF10	109	60	70M	=	581	542	CGATACTCGAGGATCCAGGGATGGCGTATTATCCTTATCTAGCAACTGTCCTAACAGTTTTGTTCAGGTT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:4	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:51	XS:i:0
RF10_128_592_1:2:0_1:0:0_27	163	RF10	128	60	70M	=	523	465	GATGGCGTATTATCCTTAAATAGCATCTGTCCTAACAGTTTTGTTCAGGTTGCACAAAGCATCTATTCCA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:55	XS:i:0
RF10_133_733_1:2:0_1:0:0_f	99	RF10	133	60	70M	=	664	601	CGTATTATCCTTATATAGCAACTGTCCTAACAGTTTTGTTCAGGTTGCACAAAGCATCTATTCCAACAAT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:55	XS:i:0
RF10_137_727_1:2:0_2:0:0_8	163	RF10	137	60	70M	=	658	591	TTATCCTTATATAGCATCTGTCCTAACAGTTTTGTTCAGGTTGCACAAAGCATCTATTGCAACAATGAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:57	XS:i:0
RF10_138_562_1:2:0_1:0:0_1e	163	RF10	138	60	70M	=	493	425	TATCCTTATATAGCATCTGTCCTAACAGTTATGTTCAGGTTGCACAAAGCATCTATTCCAACAATGAAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S4	NM:i:3	XP:Z:RF10_139_T_A;RF10_175_C_G	AS:i:58	XS:i:0
```


paired mode:

```
$ samtools sort -n -T xx src/test/resources/S2.bam | \
	java -jar dist/bamphased01.jar -V src/test/resources/rotavirus_rf.vcf.gz   --min-supporting 2 --paired --ignore-discordant-rg

@HD	VN:1.6	SO:queryname
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
@RG	ID:S2	SM:S2	LB:L2	CN:Nantes
@PG	ID:0	PN:bamphased01	VN:5d263eb	CL:-V src/test/resources/rotavirus_rf.vcf.gz --min-supporting 2 --paired --ignore-discordant-rg
@CO	bamphased01. compilation:20210220193801 githash:5d263eb htsjdk:2.23.0 date:20210220193911. cmd:-V src/test/resources/rotavirus_rf.vcf.gz --min-supporting 2 --paired --ignore-discordant-rg
RF03_2142_2592_1:1:0_2:1:0_1c	99	RF03	2142	60	70M	=	2523	451	AACTATTTTATTGGAATTAAGTTCAAAAATATACCCTATGAATATGATGATAAAGTACCCCATCTTACAT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF03_2201_G_C;RF03_2573_A_G	AS:i:60	XS:i:0
RF03_2142_2592_1:1:0_2:1:0_1c	147	RF03	2523	60	70M	=	2142	-451	ATAAAAGGCGACACACTGTTAGATATGACTGAGTGAGCTAAAAACTTAACGCACTGGTCACATCGTGACC	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF03_2201_G_C;RF03_2573_A_G	AS:i:55	XS:i:0
RF05_813_1307_2:1:0_1:1:0_10	99	RF05	813	60	70M	=	1238	495	AAAGCGTGAACGCATATAGTTGATGCTAGAAATTATATTAGTATTATGAACTCATCGTATACTGAGAATT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:56	XS:i:0
RF05_813_1307_2:1:0_1:1:0_10	147	RF05	1238	60	70M	=	813	-495	CATAAATGGAACGGAGTATGTATTATTAGACTATGAAGTGAACTGGGAAGTGAGGGGACGAGTCAGGCAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:60	XS:i:0
RF05_829_1327_3:1:0_4:1:0_16	99	RF05	829	60	70M	=	1258	499	TCGTTGATGCTCGAAATTATATTAGTATTATGAAGTCATCGTATACTGAGAATTACAGTGTGTCACAAAG	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:4	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:53	XS:i:0
RF05_829_1327_3:1:0_4:1:0_16	147	RF05	1258	60	70M	=	829	-499	TATTATTAGACTATCAAGTGAACTGGGAAGTGAGGGGACGAGTCATGCATAACAGGGATGCGAAAGTACC	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:5	XP:Z:RF05_879_C_A;RF05_1297_T_G	AS:i:45	XS:i:0
RF05_843_1306_0:1:0_2:1:0_23	83	RF05	1237	60	70M	=	843	-464	ACATAAATCGAACCGAGTATGTATTATTAGACTATGAAGTGAACTGGGAAGTGAGGGGACGAGTCATGCA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:55	XS:i:0
RF05_843_1306_0:1:0_2:1:0_23	163	RF05	843	60	70M	=	1237	464	AATTATATTAGTATTATGAACTCATCGTATACTGAGAATTACAGTGTGTCACAAAGATGTAAATTGTTTA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:1	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:65	XS:i:0
RF05_857_1394_2:1:0_1:1:0_2f	99	RF05	857	60	70M	=	1325	538	TATGAACTCATCGTATACTGAGAATTACAGTGTGTCACAAAGATGTACATTGATTACTAAGTATAAATTT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:55	XS:i:0
RF05_857_1394_2:1:0_1:1:0_2f	147	RF05	1325	60	70M	=	857	-538	TCCAAGAATTTTGACTATGAATGATACAAAGAAGATACTGAGTGCAATGATATTTGACTGGTTTGACACA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:64	XS:i:0
RF05_863_1322_3:1:0_0:1:0_0	83	RF05	1253	60	70M	=	863	-460	GTATGTATTATTAGACTATGAAGTGAACTGGGAAGTGAGGGGACGAGTCATGCAAAACATGGATGGGAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:1	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:65	XS:i:0
RF05_863_1322_3:1:0_0:1:0_0	163	RF05	863	60	70M	=	1253	460	CTCATCGTATACTGAGAATTACAGTTTGTCACAAAGATGTAAAATGTTTACTAAGTATAAATATGGGATT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:4	XP:Z:RF05_1297_T_G;RF05_879_C_A	AS:i:50	XS:i:0
RF05_874_1375_3:1:0_1:1:0_4d	99	RF05	874	60	70M	=	1306	502	CTGAGAATTACAGTGTGTCACAAAGATTTAAATTGTTTACTAAGTATACATTTGGGATTGTATCAAGATA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:4	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:54	XS:i:0
RF05_874_1375_3:1:0_1:1:0_4d	147	RF05	1306	60	70M	=	874	-502	AAAACATGGATGGGAAAGTACCAAGAATTTTGACTATGAATGATACAAAGAAGATACTCAGTGCAATGAT	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:2	XP:Z:RF05_879_C_A;RF05_1339_A_C	AS:i:60	XS:i:0
RF05_916_1354_3:0:0_1:2:0_53	147	RF05	1285	60	70M	=	916	-439	AAGTGAGGGGAAGAGTCATGCAAAACATGGATGGGAAAGTACCAAGAATTTTGACTATGAATGATACAAA	2222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	RG:Z:S2	NM:i:3	XP:Z:RF05_1297_T_G;RF05_1339_A_C	AS:i:55	XS:i:0
```

END_DOC
*/
@Program(name="bamphased01",
description="Extract Reads from a SAM/BAM file supporting at least two variants in a VCF file.",
keywords={"vcf","phased","genotypes","bam"},
creationDate="20210218",
modificationDate="20210218"
)
public class BamPhased01 extends OnePassBamLauncher {
	private static final Logger LOG=Logger.build(BamPhased01.class).make();

	
	private static class PosToCheck extends SimplePosition {
		final String ref;
		final Set<Byte> alts;
		PosToCheck(final VariantContext vc, final Set<Byte> alts) {
			super(vc.getContig(),vc.getStart());
			this.ref = vc.getReference().getDisplayString();
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
	@Parameter(names={"--tag"},description="Tag in metadata of read containing the variants positions. Ignored if empty")
	private String XTAG = "XP";
	@Parameter(names={"--mapq0"},description="If set. Do not remove the reads failing the test but set the MAPQ to 0.")
	private boolean failing_mapq0 = false;
	@Parameter(names={"--min-supporting"},description="Min number of variants that should be supported by one read.")
	private int num_supporting_variants = 2;
	@Parameter(names={"--paired"},description="Activate Paired-end mode. Variant can be supported by the read or/and is mate. Input must be sorted on query name using for example 'samtools collate'.")
	private boolean paired_mode = false;
	@Parameter(names={"--ignore-discordant-rg"},description="In paired mode, ignore discordant read-groups RG-ID.")
	private boolean ignore_discordant_rg = false;

	
	private	VCFReader vcfReader = null;
	private BufferedVCFReader bufferedVCFReader = null;
	private Set<String> samplesInBam = new HashSet<>();
	private SAMProgramRecord samProgramRecord = null;
	
	private VariantContext simplify(final VariantContext vc) {
		if(!vc.isSNP()) return null;
		if(!vc.isBiallelic()) return null;
		if(!samplesInBam.stream().
				map(S->vc.getGenotype(S)).
				anyMatch(G->G!=null && !G.isHomRef() && !G.isNoCall()))
				return null;
		
		return new VariantContextBuilder(vc).noID().passFilters().
				log10PError(VariantContext.NO_LOG10_PERROR).
				attributes(Collections.emptyMap()).
				make();
		}
	
	@Override
	protected int beforeSam() {
		if(!(this.XTAG.length()==0 || this.XTAG.length()==2)) {
			LOG.error("tag should be empty of length==2 but got "+this.XTAG);
			return -1;
			}
		if(this.XTAG.length()==2 && !this.XTAG.startsWith("X")) {
			LOG.error("tag should start with 'X' but got "+this.XTAG);
			return -1;
			}
		if(this.num_supporting_variants<2) {
			LOG.error("Bad number of supporting variant (should be >=2) "+ this.num_supporting_variants);
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
		final SAMSequenceDictionary vcfDict = this.vcfReader.getHeader().getSequenceDictionary();
		if(vcfDict!=null) {
			SequenceUtil.assertSequenceDictionariesEqual(vcfDict, SequenceDictionaryUtils.extractRequired(headerIn));
			}
		
		if(this.paired_mode) {
			final SAMFileHeader.SortOrder order= headerIn.getSortOrder();
			switch(order) {
				case queryname: break;
				case unsorted: break;
				default:throw new SAMException("in paired mode input must be sorted on query name. But got " + order);
				}
			}
	
	
	
		final SAMFileHeader outHeader= super.createOutputHeader(headerIn);
		
		this.samProgramRecord = outHeader.createProgramRecord();
		this.samProgramRecord.setProgramName(this.getProgramName());
		this.samProgramRecord.setProgramVersion(this.getVersion());
		this.samProgramRecord.setCommandLine(this.getProgramCommandLine());
		if(this.paired_mode) outHeader.setSortOrder(SAMFileHeader.SortOrder.unknown);
		
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
	
	
	private void failingSAMRecord(final SAMRecord rec,final SAMFileWriter sfw) {
		if(!failing_mapq0) return;
		rec.setMappingQuality(0);
		rec.setAttribute("PG", this.samProgramRecord.getId());
		sfw.addAlignment(rec);
		}
	
	@Override
	protected void scanIterator(final SAMFileHeader headerIn,final CloseableIterator<SAMRecord> iter0,final SAMFileWriter sfw) {
		if(this.paired_mode) {
			try(EqualIterator<SAMRecord> iter = new EqualIterator<>(iter0,(A,B)->A.getReadName().compareTo(B.getReadName()))) {
			while(iter.hasNext()) {
					final LinkedList<SAMRecord> buffer = new LinkedList<>(iter.next());
					SAMRecord R1= null;
					SAMRecord R2= null;
					while(!buffer.isEmpty()) {
						final SAMRecord rec = buffer.pop();
						if(rec.getReadUnmappedFlag()){
							failingSAMRecord(rec,sfw);
							}
						else if(!rec.getReadPairedFlag() || rec.isSecondaryOrSupplementary()) {
							scanVariants(Collections.singletonList(rec),sfw);
							}
						else if(R1==null && rec.getFirstOfPairFlag()) {
							R1 = rec;
							}
						else if(R2==null && rec.getSecondOfPairFlag()) {
							R2 = rec;
							}
						else
							{
							failingSAMRecord(rec,sfw);
							}
						}
					if(R1!=null && R2!=null) {
						if(R1.contigsMatch(R2)) {
							scanVariants(Arrays.asList(R1,R2),sfw);
							}
						else
							{
							scanVariants(Collections.singletonList(R1),sfw);
							scanVariants(Collections.singletonList(R2),sfw);
							}
						}
					else if(R1!=null && R2==null) {
						scanVariants(Collections.singletonList(R1),sfw);
						}
					else if(R2!=null && R1==null) {
						scanVariants(Collections.singletonList(R2),sfw);
						}
					}
				}
			}
		else
			{
			while(iter0.hasNext()) {
				final SAMRecord rec = iter0.next();
				if(rec.getReadUnmappedFlag()){
					failingSAMRecord(rec,sfw);
					continue;
					}
				scanVariants(Collections.singletonList(rec),sfw);
				}
			}
		}

	
	private void scanVariants(final List<SAMRecord> buffer0,final SAMFileWriter sfw) {
			if(buffer0.isEmpty()) return;
			final List<SAMRecord> buffer = new ArrayList<>(buffer0);
			
			final SAMReadGroupRecord rg0 = buffer.get(0).getReadGroup();
			if(rg0==null) {
				for(SAMRecord rec:buffer) failingSAMRecord(rec, sfw);
				return;
				}
			final String rgid0 = buffer.get(0).getReadGroup().getId();
			for(int i=1;!this.ignore_discordant_rg && i< buffer.size();i++) {
				if(buffer.get(i).getReadGroup()==null || 
					!buffer.get(i).getReadGroup().getId().equals(rgid0)) 
					{
					/*
					throw new SAMException("Two paired read without the same RG-ID:\n"+ 
							buffer.get(0).getSAMString()+"\n"+
							buffer.get(i).getSAMString()
							);
					*/
					}
				}
			
			final List<PosToCheck> supporting = new ArrayList<>();

			
			int i=0;
			while(i<buffer.size()) {
				final SAMRecord rec = buffer.get(i);
				final byte[] bases = rec.getReadBases();
				if(bases==null || bases==SAMRecord.NULL_QUALS || bases.length==0) {
					failingSAMRecord(rec,sfw);
					buffer.remove(i);
					continue;
					}
				final String sn = rg0.getSample();
				if(StringUtils.isBlank(sn) || !this.samplesInBam.contains(sn)) {
					failingSAMRecord(rec,sfw);
					buffer.remove(i);
					continue;
					}
				
				final List<PosToCheck> candidates = new ArrayList<>();

				try(CloseableIterator<VariantContext> iter=this.bufferedVCFReader.query(rec)) {
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						final Genotype gt = ctx.getGenotype(sn);
						if(gt.isHomRef() || gt.isNoCall()) continue;
						final Set<Byte> alts = gt.getAlleles().stream().
							filter(A->A.isCalled() && !A.isReference() && !A.isSymbolic() && A.length()==1 ).
							filter(A->AcidNucleics.isATGC(A)).
							map(A->(byte)Character.toUpperCase(A.getDisplayString().charAt(0))).
							collect(Collectors.toSet());
						if(alts.isEmpty()) continue;
						final PosToCheck pos = new PosToCheck(ctx,alts);
						candidates.add(pos);
						}
					}
				if(candidates.isEmpty()) {
					failingSAMRecord(rec,sfw);
					buffer.remove(i);
					continue;
					}
				
				final List<PosToCheck> supporting0 = new ArrayList<>();

				for(AlignmentBlock ab:rec.getAlignmentBlocks()) {
					final int readPos1 = ab.getReadStart();
					final int refPos1 = ab.getReferenceStart();
					for(int x=0;x< ab.getLength();++x) {
						for(PosToCheck pos:candidates) {
							if(pos.getPosition() != refPos1+x) continue;
							final byte readBase = bases [ (readPos1-1) + x ];
							if(pos.alts.contains(readBase)) {
								supporting0.add(pos);
								break;
								}
							}
						}
					}
				if(supporting0.isEmpty()) {
					failingSAMRecord(rec,sfw);
					buffer.remove(i);
					continue;
					}
				supporting.addAll(supporting0);
				i++;
				}
			
			if(supporting.size() < this.num_supporting_variants) {
				for(SAMRecord rec:buffer) failingSAMRecord(rec, sfw);
				return;
				}
			for(final SAMRecord rec:buffer) {
				if(!StringUtils.isBlank(this.XTAG)) {
					rec.setAttribute(this.XTAG,
							supporting.stream().
								map(S->S.getContig()+"_"+S.getStart()+"_"+S.ref+"_"+S.alts.stream().map(B->""+(char)B.byteValue()).collect(Collectors.joining("_"))).
								collect(Collectors.joining(";"))
							);
						}
	
				rec.setAttribute("PG", this.samProgramRecord.getId());
				sfw.addAlignment(rec);
				}
				
			
		}
	
	
	public static void main(final String[] args) {
		new BamPhased01().instanceMainWithExit(args);

	}

}
