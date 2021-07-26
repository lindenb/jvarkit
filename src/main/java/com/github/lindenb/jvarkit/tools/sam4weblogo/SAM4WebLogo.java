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
package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.SAMRecordNaming;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.filter.SamRecordFilter;
/**
BEGIN_DOC

## Motivation

"Sequence logo ( http://weblogo.berkeley.edu/logo.cgi ) for different alleles or generated from SAM/BAM" http://www.biostars.org/p/73021

![ScreenShot](https://raw.github.com/lindenb/jvarkit/master/doc/sam2weblogo.png)

## History

On October 14th, I tried to implement the insertions. I haven't tested this feature in depth.

## Example

```bash
$ java -jar dist/sam4weblogo.jar -r seq1:80-110  sorted.bam  2> /dev/null | head -n 50
>B7_593:4:106:316:452/1
TGTTG--------------------------
>B7_593:4:106:316:452a/1
TGTTG--------------------------
>B7_593:4:106:316:452b/1
TGTTG--------------------------
>B7_589:8:113:968:19/2
TGGGG--------------------------
>B7_589:8:113:968:19a/2
TGGGG--------------------------
>B7_589:8:113:968:19b/2
TGGGG--------------------------
>EAS54_65:3:321:311:983/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983a/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983b/1
TGTGGG-------------------------
>B7_591:6:155:12:674/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674a/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674b/2
TGTGGGGG-----------------------
>EAS219_FC30151:7:51:1429:1043/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043a/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043b/2
TGTGGGGGGCGCCG-----------------
>B7_591:5:42:540:501/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501a/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410a/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501b/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410b/1
TGGGGGGGGCGCAGT----------------
```

### fastq-like output

```
$ java -jar dist/sam4weblogo.jar -r 'RF01:100-130' src/test/resources/S1.bam --format fastq -c

@RF01_44_622_1:0:0_1:0:0_3a/1
TATTCTTCCAATAG-----------------
+
22222222222222                 
@RF01_44_499_0:0:0_3:0:0_7b/2
TATTCTTCCAATAG-----------------
+
22222222222222                 
@RF01_67_565_0:0:0_2:0:0_67/2
TATTCTTCCAATAGTGAATTAGAGAATAGAT
+
2222222222222222222222222222222
@RF01_94_620_1:0:0_2:0:0_15/2
TATTCTTCCAATAGTGAATTAGAGAATAGAT
+
2222222222222222222222222222222
@RF01_102_665_1:0:0_1:0:0_71/1
--TTCTTCCAATAGTGAATTAGAGAATAGAT
+
  22222222222222222222222222222
@RF01_110_504_2:0:0_1:0:0_5d/2
----------ATAGTGAATTAGATAATAGAT
+
          222222222222222222222
@RF01_121_598_1:0:0_3:0:0_6e/2
---------------------GAGAATAGAT
+
                     2222222222
```
## tabular output

```
$ java -jar /home/lindenb/src/jvarkit-git/dist/sam4weblogo.jar  -r 'chr4:15648762-15648862'    src/test/resources/retrocopy01.bwa.bam --format tabular -c
CTTGAACCCAGGAGGGGGAGGTGCCAGGGAGCCGAGATCATGCCACTGCACCCCAGCCTGGGCAACAAACCAAACCTCCATCC-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1111:24870:51570/1
CTTGAACCCAGGAGGCAGAGGTGGCAGTGAGCAGAAAACAAGCCACTGCACCCCAGCCGGGGAAAAAAAACAAGACCCCATCTaA-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1107:19532:10468/1
CTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2121:23043:26536/2
ctggacCCCAGGAGGCAGAGGTTGCAGGTACCCGAAATAATGCCCCTGCCCCCCAGCCTGGGCAACAAAGCAAGACCCCATCT-A-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1212:10642:64966/2
-------------GGCAGAGGTTGCAGTGAGCCGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-CaAAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1119:27661:8095/2
--------------------------------CGAGATCATGCCACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2224:26189:58778/2
-------------------------------------------CACTGCACTCCAGCCTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2103:11972:45611/1
ctggacccaaggaggcagaggttgaagggagcagaaacaaggccacggccccCCAGCCTGGGAAACAAAGCAAACCCCAATCT-A-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2121:29701:6108/1
---------------------------------------------------------CTGGGCAACAAAGCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1114:11992:43044/1
---------------------------------------------------------CTGGGCAACAAAGCAAGACTCCATCT-CaAAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2109:9019:23477/2
----------------------------------------------------------------ACAAAGCAAGACTCCATCT-CaAAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2110:15706:60765/2
---------------------------------------------------------------------GCAAGACTCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1218:27600:65793/1
-------------------------------------------------------------------------GACNCCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1209:1326:20805/1
-----------------------------------------------------------------------------CCATCT-C-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2201:25895:8095/2
gggccgatggtgaaaagtggcagggcgcaacaaaaaatgcccgtcacaccaaccgggggcaaaaaaaaaaaaccccaaccaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1102:32319:57829/2
-----------------------------------------------------agcagggaaaacaaaccaaaacacaacata-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1106:29437:72033/1
--------------------------------------------------------------------------------------AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1109:27813:4104/1
------------------------------------------ccacgaccccaccacggaaaaaaaaaaaaaacacaaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1112:17137:65951/1
tttaaaccaagacgccagacgtccaagatatacgactacccgccactgcaaaaccacctgtgcacaaacacacaaataaaact-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1117:10693:25042/2
------------------------------------------------------caccagggaaaaaaaaaaaaactcaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1118:24962:18502/1
ttgacaccagaagcaaaagttaaaagtcaaccgaaatcagcaaaataatctacaacgtaggaaaaagccaaagtctcattcta-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1119:15432:30791/1
tgaaagacggggggcagaggttgcaagagcccaaaaacatgccccggcacccaaacgtgggcaacaaagaaaaacccaatcaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1123:1590:53891/1
----------------------------------------------------------------------------------a-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1215:26585:67586/1
--------------------------------------------------------------------------------------AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1216:18142:49900/1
---------------------------------------atgccactgaacccaagccagggcaaaaaagaaaaaaccaatcc-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1219:25043:11048/2
-------------------------------------------------------------------aaaaaaacacaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2113:14215:55086/1
aataaacccaaaagaaaaagtgcgcagtaggccaaaatcaagcaaatgaaatcaagccggggaaaaaaacaaagaccccacat-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2113:14418:23618/1
-------------------------------------------aaccgaaacccaacccagaaaaaaaaaaaaaacacaaaca-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2120:21420:20629/2
--------------------------------------aaagcacagacaaaacaaccagacaaacaaaaaaaaaacaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2212:28564:14301/1
tacgaaccagggaggctcgggctgaggaagccgaacttcattcctcccaccgcaagactggcaaataaatcaagtcccaattt-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2216:20324:64685/1
--------------------------------------------accaccccccacccgggccacaaaaaaaaaaccaaaaaa-a-AAAAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:2216:23429:26607/1
---------------------------------------------------------------------------------------AAAAAAAAAAAGAGAA X2:1:H7T57CCXY:8:2124:1306:6249/1
----------------------------------------------------------------------------------------AAAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1222:16579:28488/1
-----------------------------------------------------------------------------------------AAAAAAAAAAAAAA X2:1:H7T57CCXY:8:1222:17513:27890/1
-------------------------------------------------------------------------------------------AAAAAAAAAAAA X2:1:H7T57CCXY:8:2109:14316:7902/1
-------------------------------------------------------------------------------------------------AAAAAA X2:1:H7T57CCXY:8:2108:14357:62329/1
--------------------------------------------------------------------------------------------------AAAAA X2:1:H7T57CCXY:8:1124:18223:49162/1
```

## See also

* https://www.biostars.org/p/103052/
* http://www.sciencedirect.com/science/article/pii/S1874778715300210

END_DOC

 */
@Program(name="sam4weblogo",
	description="Sequence logo for different alleles or generated from SAM/BAM ",
	biostars= {73021,368200},
	creationDate="20130524",
	modificationDate="20191014",
	keywords={"sam","bam","visualization","logo"}
	)
public class SAM4WebLogo extends Launcher
	{
	private static final int NO_POS= -99999999; 
	private static final Logger LOG = Logger.build(SAM4WebLogo.class).make();
	
	private enum Format{fasta,fastq,tabular};
	
	@Parameter(names={"-c","--clipped","--clip"},description="Use Clipped Bases")
	private boolean useClip = false;
	
	@Parameter(names={"-r","--region","--interval"},description="Region to observe. "+IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class,required=true)
	private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"-readFilter","--readFilter"},description="[20171201](moved to jexl)"+SamRecordJEXLFilter.FILTER_DESCRIPTION)
	private SamRecordFilter SamRecordFilter = SamRecordJEXLFilter.buildAcceptAll();
	
	@Parameter(names={"--format","-F"},description="output format.")
	private Format output_format = Format.fasta;
	@Parameter(names={"-fqu","--fqu"},description="[20180813] fastq unknown quality character")
	private String fastq_quality_unknown_str = "!";
	@Parameter(names={"-fqp","--fqp"},description="[20180813] fastq padding quality character")
	private String fastq_quality_padding_str = "-";
	@Parameter(names={"-R","--reference"},description="For Reading CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"--no-insert"},description="Do not show insertions")
	private boolean hide_insertions = false;
	@Parameter(names={"--naming"},description=SAMRecordNaming.OPT_DESC,converter=SAMRecordNaming.StringConverter.class,splitter=NoSplitter.class)
	private SAMRecordNaming samRecordNaming = SAMRecordNaming.compile("%n (%f) %N");
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.fastq_quality_padding_str.length()!=1) {
			LOG.error("Bad fastq padding character (length!=1)");
			return -1;
		}
		
		if(this.fastq_quality_unknown_str.length()!=1) {
			LOG.error("Bad fastq unknown character (length!=1)");
			return -1;
		}
		
		PrintWriter out=null;
		SamReader samReader=null;
		SAMRecordIterator iter=null;
		try {
			

			
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) srf.referenceSequence(this.faidx);
			
			
			final List<Locatable> intervals = this.intervalListProvider.
					stream().
					collect(Collectors.toList());
			
			if(!this.output_format.equals(Format.tabular) && intervals.size()>1) {
				LOG.error("Only one interval is allowed for format "+this.output_format);
				return -1;
			}
			
	    	final List<Path> inputPath = IOUtils.unrollPaths(args);
			out = super.openPathOrStdoutAsPrintWriter(outputFile);

			
			

			
			for(int locatable_idx=0; locatable_idx < intervals.size();++locatable_idx) {
				if(locatable_idx>0) out.println();
				final Locatable interval  = intervals.get(locatable_idx);
				final List<SAMRecord> buffer = new ArrayList<>();
				final SortedMap<Integer,Integer> genomicpos2insertlen = new TreeMap<>();
					
				
				for(final Path samPath: inputPath) {
					samReader = srf.open(SamInputResource.of(samPath));
					
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samReader.getFileHeader());
					final ContigNameConverter ctgNameConverter =  ContigNameConverter.fromOneDictionary(dict);
					
					final SimpleInterval segment = ctgNameConverter.convertToSimpleInterval(interval).orElseThrow(()->new JvarkitException.ContigNotFoundInDictionary(interval.getContig(), dict));
					
					
					final int extend = useClip?1000:0;
					
					iter=samReader.query(
							segment.getContig(),
							Math.max(1,segment.getStart()-extend),
							segment.getEnd()+extend,
							false);
					while(iter.hasNext())
		                {
		                final SAMRecord rec=iter.next();
		                if(rec.getReadUnmappedFlag()) continue;
		                if(this.SamRecordFilter.filterOut(rec)) continue;
		                final Cigar cigar =rec.getCigar();
		                if(cigar==null || cigar.isEmpty()) continue;
		               
		                if(!rec.getReferenceName().equals(segment.getContig())) continue;
		               
		                if(useClip) 
		                	{
		                	if(rec.getUnclippedEnd() < segment.getStart() ) continue;
			                if(rec.getUnclippedStart() > segment.getEnd() ) continue;
		                	}
		                else
		                	{
		                	if(rec.getEnd() < segment.getStart() ) continue;
			                if(rec.getStart() > segment.getEnd() ) continue;
		                	}
		                
		                // free memory
		               if(!this.output_format.equals(Format.fastq)) {
		            	rec.setBaseQualities(SAMRecord.NULL_QUALS);   
		               	}
		               final List<SAMTagAndValue> atts = rec.getAttributes();
		               for(final SAMTagAndValue stv:atts) {
		            	   if(stv.tag.equals("RG")) continue;
		            	   rec.setAttribute(stv.tag, null);
		               }
		               
		               
	                    int refPos = rec.getStart();
	                    for(final CigarElement ce: cigar) {
	                    	final CigarOperator op = ce.getOperator();
	                    	if(op.equals(CigarOperator.I)) {
	                    		genomicpos2insertlen.put(refPos, Math.max(ce.getLength(), genomicpos2insertlen.getOrDefault(refPos, 0)));
	                    		}
	                    	if(op.consumesReferenceBases()) refPos+=ce.getLength();
	                    	if(refPos > interval.getEnd()) break;
	                    	}
	                    buffer.add(rec);
		                }
					iter.close();
					iter=null;
					samReader.close();
					samReader=null;
					}
				
				
				
				
				/* compute columns */
				if(interval!=null) {
				/** test if genomic position is in interval */
				final IntPredicate testInInterval = pos-> interval.getStart() <= pos && pos <= interval.getEnd();
	
					
				final ArrayList<Integer> column2genomic = new ArrayList<>(interval.getLengthOnReference());
				for(int x=interval.getStart();x<=interval.getEnd();++x)
					{
					column2genomic.add(x);
					}
				
				// insert insertions in reverse order
				for(int pos: genomicpos2insertlen.keySet().stream().
						sorted(Collections.reverseOrder()).
						collect(Collectors.toList()))
					{
					if(!testInInterval.test(pos)) continue;
					
					final int insert_len = genomicpos2insertlen.get(pos);
					for(int x=0;x<insert_len;++x)
						{
						column2genomic.add(pos-interval.getStart(),NO_POS);
						}
					}
					
				
				// ruler
				if(this.output_format.equals(Format.tabular))
				{
					int x=0;
					while(x< column2genomic.size()) {
						if(column2genomic.get(x)==NO_POS || column2genomic.get(x)%10!=0) {
							out.print(' ');
							++x;
							} 
						else
							{
							final String s = String.valueOf(column2genomic.get(x));
							if(s.length()<10 && x+s.length()<  column2genomic.size())
								{
								out.print(s);
								x+=s.length();
								}
							else
								{
								out.print(' ');
								++x;
								}
							}
					}
					out.println(" "+interval.toString());
				}
				
					
				for(final SAMRecord rec:buffer) {
			        final String recQualString = rec.getBaseQualityString();
				        final Function<Integer,Character> read2qual;
				        if(recQualString==null || SAMRecord.NULL_QUALS_STRING.equals(recQualString)) {
				        	read2qual = IDX -> '~';
				        	}
				        else
				        	{
				        	read2qual = IDX -> {
					        	if(IDX<0 || IDX>=recQualString.length()) return '~';
					        	return recQualString.charAt(IDX);
					        	};
				        	}
				        final byte rec_bases[] = rec.getReadBases();
				        final Function<Integer,Character> read2base;
				        if(SAMRecord.NULL_SEQUENCE.equals(rec_bases)) {
				        	read2base = IDX-> '?'; 
				        	}
				        else
				        	{
				        	read2base = 
						        IDX -> {
							        	if(IDX<0 || IDX>=rec_bases.length) return '?';
							        	return (char)rec_bases[IDX];
							        	};
				        	}    
	
						
						final StringBuilder seqBuilder = new StringBuilder(column2genomic.size());
						final StringBuilder qualBuilder = new StringBuilder(column2genomic.size());
						int x = 0;
						int readRef = rec.getUnclippedStart();
						int readpos = 0;
						
						final Cigar cigar  = rec.getCigar();
						for(final CigarElement ce:cigar)
							{
							final CigarOperator op = ce.getOperator();
							if(op.equals(CigarOperator.PADDING)) continue;
							/* IN INSERTION, only print if showInsertions is true */
							
							
							for(int cigarIdx =0; cigarIdx < ce.getLength();++cigarIdx)
								{
								//pad before base
								while( !op.equals(CigarOperator.INSERTION) &&
									x < column2genomic.size() &&
									(column2genomic.get(x)==NO_POS || column2genomic.get(x) < readRef))
									{
									seqBuilder.append('-');
									qualBuilder.append('-');
									++x;
									}
								switch(op) {
									case I:
										{
										if(!this.hide_insertions) {
											if(testInInterval.test(readRef)) {
												seqBuilder.append(Character.toLowerCase(read2base.apply(readpos)));
												qualBuilder.append(read2qual.apply(readpos));
												++x;
												}
											}
										++cigarIdx;
										++readpos;
											
										break;
										}
									case P: break;
									case H:
										{
										if(testInInterval.test(readRef)) {
											seqBuilder.append(this.useClip?'n':'-');
											qualBuilder.append(this.useClip?this.fastq_quality_unknown_str:"-");
											++x;
											}
										++readRef;
										break;
										}
									case S:
										{
										if(testInInterval.test(readRef)) {
											seqBuilder.append(this.useClip?Character.toLowerCase((read2base.apply(readpos))):'-');
											qualBuilder.append(this.useClip?(char)read2qual.apply(readpos):'-');
											++x;
											}
										++readpos;
										++readRef;
										break;
										}
									case D: case N:
										{
										if(testInInterval.test(readRef)) {
											seqBuilder.append('>');
											qualBuilder.append('>');
											++x;
											}
										++readRef;
										break;
										}
									case X: case EQ: case M:
										{
										if(testInInterval.test(readRef)) {
											seqBuilder.append(read2base.apply(readpos));
											qualBuilder.append(read2qual.apply(readpos));
											++x;
											}
										++readpos;
										++readRef;
										break;
										}
									default : throw new IllegalArgumentException(op.name());
									}
								} // end cigarIdx
							}// end while cigar
							
						/* right pad */
						while(x < column2genomic.size())
							{
							seqBuilder.append('-');
							qualBuilder.append('-');
							++x;
							}
						
						final String readName= this.samRecordNaming.apply(rec);
						switch(this.output_format)
							{
							case fastq:
								out.print(FastqConstants.SEQUENCE_HEADER);
						    	out.print(readName);
						    	out.println();
						    	out.print(seqBuilder);
						    	out.println();
						    	out.println(FastqConstants.QUALITY_HEADER);
						    	out.println(qualBuilder);
								break;
							case tabular:
						    	out.print(seqBuilder);
						    	out.print(' ');
						    	out.println(readName);
								break;
							default:
								out.print('>');
						    	out.print(readName);
						    	out.println();
						    	out.print(seqBuilder);
						    	out.println();
								break;
							}
						
						}	 //end loop over SAMRec
					} // end if interval !=null
				}
			out.flush();
	        out.close();out=null;
	        return 0;
			} 
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(out);
			}
		}

public static void main(final String[] args)
	{
	new SAM4WebLogo().instanceMainWithExit(args);
	}
}
