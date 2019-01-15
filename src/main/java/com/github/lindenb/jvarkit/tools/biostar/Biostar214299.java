/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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

## History:

* 20180808: fixing bug https://github.com/lindenb/jvarkit/issues/108
* 20180212: fixing bug https://github.com/lindenb/jvarkit/issues/95

END_DOC
 */
@Program(name="biostar214299",
	description="Extract allele specific reads from bamfiles",
	biostars=214299,
	keywords={"sam","bam","variant","snp"}
	)
public class Biostar214299 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar214299.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-p","--positions"},description="Position file. A Tab delimited file containing the following 4 column: (1)chrom (2)position (3) allele A/T/G/C (4) sample name.",required=true)
	private File positionFile = null;
	
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
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
	public int doWork(final List<String> args) {
		if(this.positionFile==null) {
			LOG.error("position File is not defined.");
			return -1;
			}
		final String UNAFFECTED_SAMPLE="UNAFFECTED";
		final String AMBIGOUS_SAMPLE="AMBIGOUS";
		final String UNMAPPED="UNMAPPED";
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		final IntervalTreeMap<Position> positionsTreeMap = new IntervalTreeMap<>();
		final Set<String> samples = new HashSet<>();
		try
			{
			sfr = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header=sfr.getFileHeader();
			final SAMSequenceDictionary dict =header.getSequenceDictionary();
			if(dict==null) 
				{
				LOG.error("Dictionary missing in input sam");
				return -1;
				}
			
			
			try ( BufferedReader br = IOUtils.openFileForBufferedReading(this.positionFile)) {
				String line;
				final CharSplitter tab = CharSplitter.TAB;
				while((line=br.readLine())!=null) {
					if(line.trim().isEmpty() || line.startsWith("#")) continue;
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
						LOG.error("in "+line+" bases should be one letter an ATGC");
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
			newHeader.addComment("generated with "+getProgramName()+" "+getVersion()+" Pierre Lindenbaum : "+getProgramCommandLine());
			/* create groups */
			for(final String sample: samples) {
				final SAMReadGroupRecord rg = new SAMReadGroupRecord(sample);
				rg.setSample(sample);
				rg.setLibrary(sample);
				newHeader.addReadGroup(rg);
			}
			
			
			sfw = this.writingBamArgs.openSAMFileWriter(this.outputFile,newHeader, true);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header).logger(LOG);
			final SAMRecordIterator iter = sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.watch(iter.next());
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
				} else if(selectedSamples.size()==1) {
					rec.setAttribute("RG",selectedSamples.iterator().next());
				} else
				{
					rec.setAttribute("RG",AMBIGOUS_SAMPLE);
				}
				
				sfw.addAlignment(rec);
				}
			progress.finish();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	public static void main(final String[] args) {
		new Biostar214299().instanceMainWithExit(args);

	}

}
