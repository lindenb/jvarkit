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
package com.github.lindenb.jvarkit.tools.splitbam;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.AbstractBamSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;


/**

BEGIN_DOC


Split a BAM by chromosome group. Create EMPTY bams if no reads was found for a given group.
![img](https://chart.googleapis.com/chart?chl=+digraph+G+%7B%0D%0ABWA+-%3E+SPLITBAM%5Blabel%3D%22stdout%22%5D%3B%0D%0A+++SPLITBAM-%3ECHR1_bam%3B%0D%0A+++SPLITBAM-%3ECHR2_bam%3B%0D%0A+++SPLITBAM-%3ECHR3_bam%3B+%0D%0A+++CHR1_bam+-%3E+CHR1_vcf%3B%0D%0A+++CHR2_bam+-%3E+CHR2_vcf%3B%0D%0A+++CHR3_bam+-%3E+CHR3_vcf%3B%0D%0A+++CHR1_vcf+-%3E+merged_vcf%3B%0D%0A+++CHR2_vcf+-%3E+merged_vcf%3B%0D%0A+++CHR3_vcf+-%3E+merged_vcf%3B%0D%0A+%7D%0D%0A++++++++&cht=gv)


### Example


the content of 'split_g1k_v37_01.txt'



```

CHROMS_01_09	1 2 3 4 5 6 7 8 9
CHROMS_10_0Y	10 11 12 13 14 15 16 17 18 19 20 21 22 X Y 
CHROMS_OTHER	MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 

```



split the output of bwa sampe on the fly:



```

bwa mem (...) | samtools sort (...) | \
java -jar dist/splitbam3.jar \
	-o TESTSPLITBAM/__GROUPID__.bam \
	-m \
	-g split_g1k_v37_01.txt 


[Fri Jul 26 13:25:56 CEST 2013] Executing as lindenb@master on Linux 2.6.32-358.6.2.el6.x86_64 amd64; OpenJDK 64-Bit Server VM 1.7.0_19-mockbuild_2013_04_17_19_18-b00; Picard version: null
INFO	2013-07-26 13:25:56	SplitBam	reading stdin
INFO	2013-07-26 13:25:56	SplitBam	opening TESTSPLITBAM/CHROMS_01_09.bam
INFO	2013-07-26 13:25:57	SplitBam	opening TESTSPLITBAM/CHROMS_10_0Y.bam
INFO	2013-07-26 13:25:58	SplitBam	opening TESTSPLITBAM/CHROMS_OTHER.bam
INFO	2013-07-26 13:35:58	SplitBam	closing group CHROMS_01_09
INFO	2013-07-26 13:35:59	SplitBam	closing group CHROMS_10_0Y
INFO	2013-07-26 13:35:59	SplitBam	closing group CHROMS_OTHER
INFO	2013-07-26 13:36:00	SplitBam	closing group Unmapped
Runtime.totalMemory()=1916600320

```







END_DOC
*/
@Program(name="splitbam3",
	description="Split a BAM by chromosome group",
	creationDate="20150317",
	modificationDate="20210402",
	keywords={"sam","bam","split"}
	)
public class SplitBam3 extends AbstractBamSplitter<String> {
	private static final Logger LOG = Logger.build(SplitBam3.class).make();
	private final static String OTHER_NAME="OTHERS";

	@Parameter(names={"-g","--groupfile"},description="Chromosome group file. Interval are 1 based. Otherwise split per chromosome.")
	private Path chromGroupFile = null;
	@Parameter(names={"--others"},description="save mapped reads that don't belong to any group into an extra group")
	private boolean other_flags = false;
	@Parameter(names={"--unmapped"},description="save unmapped reads into an extra group")
	private boolean unmapped_flags = false;


	private final IntervalTreeMap<String> interval2group = new IntervalTreeMap<>();
	
	private SplitBam3() {
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	@Override
	protected int beforeIterator(final SAMFileHeader header) {
		try
			{
			final SAMSequenceDictionary samSequenceDictionary= SequenceDictionaryUtils.extractRequired(header);	
		
			if(this.chromGroupFile!=null)
				{
				final Set<String> set = new HashSet<>();
				try(BufferedReader r=IOUtils.openPathForBufferedReading(this.chromGroupFile)) {
				String line;
				while((line=r.readLine())!=null)
					{
					if(StringUtils.isBlank(line)|| line.startsWith("#")) continue;
					final String tokens[] =line.split("[ \t,]+");
					final String groupName=tokens[0].trim();
					if(StringUtils.isBlank(groupName)) throw new IllegalArgumentException("Empty group name in "+line);
					if(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(groupName))  throw new IOException("Group cannot be named "+SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
					if(OTHER_NAME.equals(groupName))  throw new IOException("Group cannot be named " + OTHER_NAME);
					if(set.contains(groupName))  throw new IOException("Group defined twice "+groupName);
					for(int i=1;i< tokens.length;i++)
						{
						String sequence;
						int start;
						int end;
						String segment = tokens[i].trim();
						
						if(StringUtils.isBlank(segment)) continue;
						
						int colon= segment.indexOf(':');
						if(colon==-1)
							{
							final SAMSequenceRecord ssr=samSequenceDictionary.getSequence(segment);
							if(ssr==null)
								{
								throw new JvarkitException.ContigNotFoundInDictionary(segment,samSequenceDictionary);
								}
							sequence = segment;
							start = 1;
							end = ssr.getSequenceLength();
							}
						else
							{
							int hyphen  = segment.indexOf('-',colon);
							if(hyphen==-1)  throw new IOException("Bad segment:"+segment);
							sequence = segment.substring(0,colon);
							if(samSequenceDictionary.getSequence(sequence)==null) {
								throw new JvarkitException.ContigNotFoundInDictionary(segment,samSequenceDictionary);
								}
							
							//+1 because interval are 1-based
							start = Integer.parseInt(segment.substring(colon+1,hyphen));
							end = Integer.parseInt(segment.substring(hyphen+1));
							}
						
						final Interval interval = new Interval(sequence, start, end);
						if(interval2group.containsKey(interval) && !
								interval2group.get(interval).equals(groupName)) {
							LOG.error("interval defined twice:"+interval);
							}
						this.interval2group.put(interval, groupName);
						set.add(groupName);
						}
					}
				}
				}
			else
				{
				// creating default split interval
				for(final SAMSequenceRecord seq:samSequenceDictionary.getSequences())
					{
					final String groupName=seq.getSequenceName();
					final Interval interval= new Interval(groupName, 1, seq.getSequenceLength());
					this.interval2group.put(interval, groupName);
					}
				}
			return 0;
			}
		catch(final Throwable err ) {
			getLogger().error(err);
			return -1;
			}
		}
	
	@Override
	protected Set<String> createKeys(final SAMRecord record) {
		final Interval interval;
		if( record.getReadUnmappedFlag() )
			{
			if(record.getReadPairedFlag() && !record.getMateUnmappedFlag())
				{
				interval= new Interval(
						record.getMateReferenceName(),
						record.getMateAlignmentStart(),
						record.getMateAlignmentStart()
						);
				}
			else if(this.unmapped_flags) {
				return Collections.singleton(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				}
			else
				{
				return Collections.emptySet();	
				}
			}
		else
			{
			interval= new Interval(record);
			}
		
		final Collection<String> groupIds = this.interval2group.getOverlapping(interval);
		if(groupIds.isEmpty()) {
			if(this.other_flags) return Collections.singleton(OTHER_NAME);
			return Collections.emptySet();
			}
		return new HashSet<>(groupIds);
		}
			
	
	public static void main(final String[] args) throws Exception
		{
		new SplitBam3().instanceMainWithExit(args);
		}
	
	}
