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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
/**
BEGIN_DOC

The program removes all the existing read group and create some new one from the 'position file'.
For now, only simple alleles are supported.
Reads group are affected if a specific variant is found in the 'position file'.
If two samples share the same group, the read group is AMBIGOUS.
If the read is unmapped, the read group is UNMAPPED.
If no sample is affected to a read, the read group will be UNAFFECTED;

## see also:

* [https://www.biostars.org/p/283969](https://www.biostars.org/p/283969)  " How to extract reads with a known variant form a bam file"


## Example

the positions file

```
$ cat positions.tsv
rotavirus       267     C       SAMPLE1
rotavirus       267     G       SAMPLE2
```

processing :

```
$ java -jar dist/biostar214299.jar -p positions.tsv input.bam

@HD     VN:1.5  SO:coordinate
@SQ     SN:rotavirus    LN:1074
@RG     ID:UNAFFECTED   SM:UNAFFECTED   LB:UNAFFECTED
@RG     ID:UNMAPPED     SM:UNMAPPED     LB:UNMAPPED
@RG     ID:SAMPLE1      SM:SAMPLE1      LB:SAMPLE1
@RG     ID:SAMPLE2      SM:SAMPLE2      LB:SAMPLE2
@RG     ID:AMBIGOUS     SM:AMBIGOUS     LB:AMBIGOUS
(...)
rotavirus_237_744_6:0:0_3:0:0_29c       163     rotavirus       237     60      70M     =       675     508     ATCCGGCGTTAAATGGAAAGTTTCGGTGATCTATTAGAAATAGAAATTGGATGACTGATTCAAAAACGGT  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      MD:Z:3A19A1C1C1G31T8    RG:Z:SAMPLE1    NM:i:6  AS:i:41 XS:i:0
rotavirus_234_692_6:0:1_4:0:0_3ac       163     rotavirus       237     60      6S30M5I1M5D28M  =       623     456     TTGGTAATCAGGCGTTAAATGGAAAGTTTAGCTCAGGACAACGAAATAGAAATTGGATGACTGATTCTAA  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      MD:Z:31^TATTA28 RG:Z:SAMPLE2    NM:i:10 AS:i:37 XS:i:0
rotavirus_237_777_6:0:0_7:0:0_216       99      rotavirus       237     60      70M     =       708     541     ATCAGGGGTTAAATTGAAAGTTTAGCTCAGCTCTTAGACATAGAAATTGGATGACTGATTGTACAACGGT  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      MD:Z:6C7G17A5A21C2A6    RG:Z:SAMPLE1    NM:i:6  AS:i:40 XS:i:0
rotavirus_237_699_3:0:0_8:0:0_22f       163     rotavirus       237     60      70M     =       650     463     ATGAGGCGTTAAATGGAAAGTTTATCTCAGCTATTAGAAATAGCAATTGGATGACTGATTCTAAAACGGT  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      MD:Z:2C21G18A26 RG:Z:SAMPLE1    NM:i:3  AS:i:57 XS:i:0
(...)
rotavirus_311_846_10:0:0_11:0:0_3d7     141     *       0       0       *       *       0       0       AACTTAGATGAAGACGATCAAAACCTTAGAATGACTTTATGTTCTAAATGGCTCGACCCAAAGATGAGAG  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      RG:Z:UNMAPPED   AS:i:0  XS:i:0
rotavirus_85_600_7:0:0_9:0:0_3e0        77      *       0       0       *       *       0       0       AGCTGCAGTTGTTTCTGCTCCTTCAACATTAGAATTACTGGGTATTGAATATGATTCCAATGAAGTCTAT  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      RG:Z:UNMAPPED   AS:i:0  XS:i:0
rotavirus_85_600_7:0:0_9:0:0_3e0        141     *       0       0       *       *       0       0       TATTTCTCCTTAAGCCTGTGTTTTATTGCATCAAATCTTTTTTCAAACTGCTCATAACGAGATTTCCACT  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      RG:Z:UNMAPPED   AS:i:0  XS:i:0
```

## Cited In

 * Anatomy, transcription dynamics and evolution of wheat ribosomal RNA loci deciphered by a multi-omics approach.  https://doi.org/10.1101/2020.08.29.273623
 * Reciprocal allopolyploid grasses (Festuca X Lolium) display stable patterns of genome dominance . Marek Glombik & al. 2021. The plant journal.  doi:10.1111/tpj.15375
 * Fine structure and transcription dynamics of bread wheat ribosomal DNA loci deciphered by a multi-omics approach. Z Tulpova & al.  2022.  The Plante Genome. https://doi.org/10.1002/tpg2.20191
 * Watson CM, Jackson L, Crinnion LA, et al. Long-read sequencing to resolve the parent of origin of a de novo pathogenic UBE3A variant. Journal of Medical Genetics Published Online First: 12 April 2022. doi: 10.1136/jmedgenet-2021-108314
 * Benjamin McClinton, Christopher M. Watson, Laura A. Crinnion, Martin McKibbin, Manir Ali, Chris F. Inglehearn, Carmel Toomes, Haplotyping Using Long-Range PCR and Nanopore Sequencing of Phase Variants; Lessons Learned From the ABCA4 Locus, Laboratory Investigation, 2023, 100160, ISSN 0023-6837
 * Shinde SS, Sharma A and Vijay N (2023) Decoding the fibromelanosis locus complex chromosomal rearrangement of black-bone chicken: genetic differentiation, selective sweeps and protein-coding changes in Kadaknath chicken. Front. Genet. 14:1180658. doi: 10.3389/fgene.2023.1180658
 * Hany U, Watson CM, Liu L, et al. Novel Ameloblastin Variants, Contrasting Amelogenesis Imperfecta Phenotypes. Journal of Dental Research. 2023;0(0). doi:10.1177/00220345231203694
 * Malena Daich Varela, Elena Schiff, Samantha Malka, Genevieve Wright, Omar A. Mahroo, Andrew R. Webster, Michel Michaelides, Gavin Arno; PHYH c.678+5G>T Leads to In-Frame Exon Skipping and Is Associated With Attenuated Refsum Disease. Invest. Ophthalmol. Vis. Sci. 2024;65(2):38. https://doi.org/10.1167/iovs.65.2.38

END_DOC
 */
@Program(name="biostar214299",
	description="Extract allele specific reads from bamfiles",
	biostars=214299,
	keywords={"sam","bam","variant","snp"},
	creationDate="20160930",
	modificationDate="20220420",
	jvarkit_amalgamion =  true,
	menu="Biostars"
	)
public class Biostar214299 extends OnePassBamLauncher {
	private static final Logger LOG = Logger.of(Biostar214299.class);

	@Parameter(names={"-p","--positions"},description="Position file. A Tab delimited file containing the following 4 column: (1)chrom (2)position (3) allele A/T/G/C (4) sample name.",required=true)
	private Path positionFile = null;
	
	private static class Position
		{
		//String contig;
		int refpos;
		final Map<Character,String> base2sample = new HashMap<>();
		@Override
		public String toString() {
			return "pos:"+this.refpos+" "+this.base2sample;
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	@Override
	protected int beforeSam() {
		if(this.positionFile==null) {
			LOG.error("position File is not defined.");
			return -1;
			}
		
		return super.beforeSam();
		}

	@Override
	protected int processInput(final SAMFileHeader header, final CloseableIterator<SAMRecord> iter) {
		final String UNAFFECTED_SAMPLE="UNAFFECTED";
		final String AMBIGOUS_SAMPLE="AMBIGOUS";
		final String UNMAPPED="UNMAPPED";
		final IntervalTreeMap<Position> positionsTreeMap = new IntervalTreeMap<>();
		final Set<String> samples = new HashSet<>();
		try
			{
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			
			try ( BufferedReader br = IOUtils.openPathForBufferedReading(this.positionFile)) {
				String line;
				final CharSplitter tab = CharSplitter.TAB;
				while((line=br.readLine())!=null) {
					if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
					final String tokens[]=tab.split(line);
					if(tokens.length<4) {
						LOG.error("Not enough columns in "+line);
						return -1;
						}
					final String contig = tokens[0];
					if(dict.getSequence(contig)==null) 
						{
						LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(contig, dict));
						return -1;
						}
					final int refpos = Integer.parseInt(tokens[1]);
					final Interval interval = new Interval(contig, refpos, refpos);
					Position position = positionsTreeMap.get(interval);
					if(position==null) {
						position = new Position();
						//position.contig = contig;
						position.refpos = refpos;
						 positionsTreeMap.put(interval, position);
					}
					
					final String bases = tokens[2].toUpperCase();
					if(bases.length()!=1 || !bases.matches("[ATGC]"))
						{
						LOG.error("in "+line+" bases should be one letter and ATGC");
						return -1;
						}
					if(position.base2sample.containsKey(bases.charAt(0))) {
						LOG.error("in "+line+" bases already defined for this position");
						return -1;
					}
					
					final String sampleName = tokens[3].trim();
					if(sampleName.isEmpty())
						{
						LOG.error("sample name cannot be empty");
						return -1;
						}
					samples.add(sampleName);
					position.base2sample.put(bases.charAt(0), sampleName);
					
					}
			} catch (final IOException err) {
				LOG.error(err);
				return -1;
			}
			
			if(samples.contains(UNAFFECTED_SAMPLE)) 
				{
				LOG.error("Sample cannot be named "+UNAFFECTED_SAMPLE);
				return -1;
				}
			if(samples.contains(AMBIGOUS_SAMPLE)) 
				{
				LOG.error("Sample cannot be named "+AMBIGOUS_SAMPLE);
				return -1;
				}
			if(samples.contains(UNMAPPED)) 
				{
				LOG.error("Sample cannot be named "+UNMAPPED);
				return -1;
				}

			samples.add(UNAFFECTED_SAMPLE);
			samples.add(AMBIGOUS_SAMPLE);
			samples.add(UNMAPPED);
			
			
		
		
			final SAMFileHeader newHeader = new SAMFileHeader();
			newHeader.setSortOrder(header.getSortOrder());
			newHeader.setSequenceDictionary(dict);
			/* create groups */
			for(final String sample: samples) {
				final SAMReadGroupRecord rg = new SAMReadGroupRecord(sample);
				rg.setSample(sample);
				rg.setLibrary(sample);
				newHeader.addReadGroup(rg);
			}
			
			
			try(SAMFileWriter sfw = super.openSamFileWriter(newHeader)) {
	
				while(iter.hasNext())
					{
					final SAMRecord rec = iter.next();
					rec.setAttribute("RG",null);
					if(rec.getReadUnmappedFlag()) {
						rec.setAttribute("RG",UNMAPPED);
						sfw.addAlignment(rec);
						continue;
						}
					
					final Cigar cigar = rec.getCigar();
					final Collection<Position> snps =  positionsTreeMap.getContained(new Interval(rec.getContig(),rec.getUnclippedStart(),rec.getUnclippedEnd()));
					if(snps== null || snps.isEmpty()) {
						rec.setAttribute("RG",UNAFFECTED_SAMPLE);
						sfw.addAlignment(rec);
						continue;
					}
					final Map<Integer,Position> index2pos= snps.stream().
							collect(Collectors.toMap(P->P.refpos,P->P));
					final Set<String> selectedSamples = new HashSet<>();
					
					final byte bases[] =rec.getReadBases();
					if(bases==null || bases.equals(SAMRecord.NULL_SEQUENCE)) {
						LOG.error("Bases missing in read "+rec);
						return -1;
						}
					
					int refPos1=rec.getUnclippedStart();
					int readPos0=0;
				
					for(final CigarElement ce:cigar.getCigarElements())
						{
						final CigarOperator op = ce.getOperator();
						final boolean consummeReadBaseOrSoftClip= 
								op.consumesReadBases() || 
								op.equals(CigarOperator.S);
						if(op.consumesReferenceBases() && consummeReadBaseOrSoftClip) 
							{
							for(int i=0;i< ce.getLength();++i){
								final int nowRefPos1 = (refPos1+i);
								final int nowReadPos0 = (readPos0+i);
								final Position position = index2pos.get(nowRefPos1);
								if(position==null) continue;
								if(nowReadPos0>= bases.length) continue;
								final char base = (char)Character.toUpperCase(bases[nowReadPos0]);
								
								final String sample = position.base2sample.get(base);
								if(sample==null) continue;
								selectedSamples.add(sample);
								
								index2pos.remove(nowRefPos1);
								if(index2pos.isEmpty()) break;
								}
							}
						if(op.consumesReferenceBases() || op.isClipping())
							{
							refPos1+=ce.getLength();
							}
						if(consummeReadBaseOrSoftClip) {
							readPos0+=ce.getLength();
							}
						}
					if(selectedSamples.isEmpty())  {
						rec.setAttribute("RG",UNAFFECTED_SAMPLE);
						}
					else if(selectedSamples.size()==1) {
						rec.setAttribute("RG",selectedSamples.iterator().next());
						}
					else {
						rec.setAttribute("RG",AMBIGOUS_SAMPLE);
						}
					
					sfw.addAlignment(rec);
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			}
		}

	public static void main(final String[] args) {
		new Biostar214299().instanceMainWithExit(args);
	}

}
