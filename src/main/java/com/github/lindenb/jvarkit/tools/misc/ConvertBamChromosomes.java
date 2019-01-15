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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**
BEGIN_DOC


## Example

```bash
$ cat samtools-0.1.19/examples/toy.sam

@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
r001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*
r003	0	ref	9	30	5H6M	*	0	0	AGCTAA	*
r004	0	ref	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*
r003	16	ref	29	30	6H5M	*	0	0	TAGGC	*
r001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*
x1	0	ref2	1	30	20M	*	0	0	aggttttataaaacaaataa	????????????????????
x2	0	ref2	2	30	21M	*	0	0	ggttttataaaacaaataatt	?????????????????????
x3	0	ref2	6	30	9M4I13M	*	0	0	ttataaaacAAATaattaagtctaca	??????????????????????????
x4	0	ref2	10	30	25M	*	0	0	CaaaTaattaagtctacagagcaac	?????????????????????????
x5	0	ref2	12	30	24M	*	0	0	aaTaattaagtctacagagcaact	????????????????????????
x6	0	ref2	14	30	23M	*	0	0	Taattaagtctacagagcaacta	???????????????????????

java -jar dist/bamrenamechr.jar \
   -f <(echo -e "ref\tCHROM1\nref2\tCHROM2")
   samtools-0.1.19/examples/toy.sam

@HD	VN:1.4	SO:unsorted
@SQ	SN:CHROM1	LN:45
@SQ	SN:CHROM2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.misc.ConvertBamChromosomes	VN:dfab75cb8c06e47e9989e59df62ec8f3242934c4	CL:-f /dev/fd/63 /commun/data/packages/samtools-0.1.19/examples/toy.sam
r001	163	CHROM1	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112
r002	0	CHROM1	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*
r003	0	CHROM1	9	30	5H6M	*	0	0	AGCTAA	*
r004	0	CHROM1	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*
r003	16	CHROM1	29	30	6H5M	*	0	0	TAGGC	*
r001	83	CHROM1	37	30	9M	=	7	-39	CAGCGCCAT	*
x1	0	CHROM2	1	30	20M	*	0	0	AGGTTTTATAAAACAAATAA	????????????????????
x2	0	CHROM2	2	30	21M	*	0	0	GGTTTTATAAAACAAATAATT	?????????????????????
x3	0	CHROM2	6	30	9M4I13M	*	0	0	TTATAAAACAAATAATTAAGTCTACA	??????????????????????????
x4	0	CHROM2	10	30	25M	*	0	0	CAAATAATTAAGTCTACAGAGCAAC	?????????????????????????
x5	0	CHROM2	12	30	24M	*	0	0	AATAATTAAGTCTACAGAGCAACT	????????????????????????
x6	0	CHROM2	14	30	23M	*	0	0	TAATTAAGTCTACAGAGCAACTA	???????????????????????

```

## See also

  * https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/g1kv37_to_hg19.tsv
  * https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/hg19_to_g1kv37.tsv
  * [[VcfRenameChromosomes]]
  * http://plindenbaum.blogspot.fr/2013/07/g1kv37-vs-hg19.html

## history

  * 20180612 : rewrote it, using a output Dict, handle the SA tag


END_DOC

 */
@Program(
	name="bamrenamechr",
	description="Convert the names of the chromosomes in a BAM file",
	keywords={"sam","bam","chromosome","contig"}
	)
public class ConvertBamChromosomes
	extends Launcher
	{
	private static final Logger LOG = Logger.build(ConvertBamChromosomes.class).make();
	
	
	@Parameter(names={"-i","-ignore"},description="If the tool cannot convert a contig, skip the read ")
	private boolean skip_on_not_found = false;
	@Parameter(names={"-f","--mapping","-m"},description="load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+",required=false)
	private File mappingFile=null;
	@Parameter(names={"-R","--reference","-r"},description="Use this reference file. "+INDEXED_FASTA_REFERENCE_DESCRIPTION,required=false)
	private File referenceFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File output= null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	private ContigNameConverter mapper  = null;
	private final Set<String> unmappedChromosomes = new HashSet<>();
	
	
	private String convert(final String contig) {
		final String newctg = this.mapper.apply(contig);
		if(newctg==null)
			{
			this.unmappedChromosomes.add(contig);
			if(skip_on_not_found) return null;
			throw new JvarkitException.ContigNotFound("Cannot convert contig "+contig);
			}
		return newctg;
		}
	
	private Map<String,String> readMappingFile(
			final SAMSequenceDictionary dictIn,
			final SAMSequenceDictionary dictOut
			) throws IOException {
		final Map<String,String> mapper= new LinkedHashMap<>();
		BufferedReader in=null;
		try {
			in=IOUtils.openFileForBufferedReading(this.mappingFile);
			String line;
			while((line=in.readLine())!=null)
				{
				if(StringUtil.isBlank(line) || line.startsWith("#")) continue;
				final String tokens[]=line.split("[\t]");
				if(tokens.length!=2
						|| tokens[0].trim().isEmpty()
						|| tokens[1].trim().isEmpty()
						) {
					in.close();in=null;
					throw new IOException("Bad mapping line: \""+line+"\"");
					}
				tokens[0]=tokens[0].trim();
				tokens[1]=tokens[1].trim();
				if(tokens[0].equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME) ||
				   tokens[1].equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
					in.close();
					throw new IOException("Bad name for: "+line);
					}
				if(dictIn.getSequence(tokens[0])==null) {
					LOG.warn("Input dictionary doesn't contain "+tokens[0]+", skipping");
					continue;
					}
				if(dictOut!=null && dictOut.getSequence(tokens[1])==null) {
					LOG.warn("Output dictionary doesn't contain "+tokens[1]+", skipping");
					continue;
					}
				
				if(mapper.containsKey(tokens[0]))
					{
					in.close();
					throw new IOException("Mapping defined twice for src \""+tokens[0]+"\"");
					}
				if(mapper.values().contains(tokens[1]))
					{
					throw new IOException("Mapping defined twice for dest \""+tokens[1]+"\"");
					}
				mapper.put(tokens[0], tokens[1]);
				}
			return mapper;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	@SuppressWarnings("deprecation")
	@Override
	public int doWork(final List<String> args) {
		final CharSplitter SEMICOLON_PAT = CharSplitter.SEMICOLON;
	   	final CharSplitter COMMA_PAT = CharSplitter.COMMA;

		if(this.referenceFile==null && this.mappingFile==null) {
			LOG.error("reference undefined and mapping file is undefined");
			return -1;
			}
		
		
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			sfr = super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader inHeader = sfr.getFileHeader();
			if(inHeader==null)
				{
				LOG.error("File header missing");
				return -1;
				}
			final SAMSequenceDictionary dictIn = SequenceDictionaryUtils.extractRequired(inHeader);
			
			if(dictIn.isEmpty())
				{
				LOG.error("input Sequence dict is empty");
				return -1;
				}
			final SAMFileHeader outHeader =inHeader.clone();
			
			//create new sequence dict
			SAMSequenceDictionary dictOut = null;
			
			Map<String,String> contigMap = null;
			
			if(this.referenceFile!=null)
				{
				dictOut  = SAMSequenceDictionaryExtractor.extractDictionary(this.referenceFile);
				}
			
			if(this.mappingFile!=null)
				{
				contigMap = this.readMappingFile(dictIn, dictOut);
				}
			
			if(dictOut==null && contigMap!=null)
				{
				final List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>(dictIn.size());
				for(final String ci : contigMap.keySet())
					{
					final SAMSequenceRecord ssr1=dictIn.getSequence(ci);
					if(ssr1==null) continue;
					final SAMSequenceRecord ssr2=new SAMSequenceRecord(contigMap.get(ci), ssr1.getSequenceLength());
					for(final Map.Entry<String,String> atts:ssr1.getAttributes())
						{	
						ssr2.setAttribute(atts.getKey(), atts.getValue());
						}
					ssrs.add(ssr2);
					}
				dictOut = new SAMSequenceDictionary(ssrs);
				this.mapper=ContigNameConverter.fromMap(contigMap);
				}
			else if(dictOut!=null && contigMap==null)
				{
				this.mapper = ContigNameConverter.fromDictionaries(dictIn, dictOut); 
				}
			if(this.mapper==null || dictOut==null)
				{
				LOG.error("Illegal state mapper==null or dict-out=null");
				return -1;
				}
			
			outHeader.setSequenceDictionary(dictOut);
			
			// check order of dictionaries, do we need to set the sort order ?
			int prev_out_index=-1;
			boolean output_is_sorted = true;
			for(int i=0;i< dictIn.size();i++)
				{
				final String si = dictIn.getSequence(i).getSequenceName();
				final String so = this.mapper.apply(si);
				if(so==null) continue;
				final int o_index = dictOut.getSequenceIndex(so);
				if(o_index < prev_out_index  ) {
					output_is_sorted = false;
					LOG.info("output will NOT be sorted");
					break;
					}
				prev_out_index= o_index;
				}
			
			if(!output_is_sorted)
				{
				outHeader.setSortOrder(SortOrder.unsorted);
				}
			
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
					newInstance().dictionary(inHeader).logger(LOG).build();
			
			sfw = this.writingBamArgs.openSAMFileWriter(this.output, outHeader, true);
			
			
			long num_ignored=0L;
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec1=progress.apply(iter.next());
				String newName1=null;
				String newName2=null;
				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					newName1= convert(rec1.getReferenceName());
					}
				if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
					{
					newName2= convert(rec1.getMateReferenceName());
					}
				rec1.setHeader(outHeader);

				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					if(newName1==null)
						{
						++num_ignored;
						continue;
						}
					else
						{
						rec1.setReferenceName(newName1);
						}
					}
				if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
					{
					if(newName2==null)
						{
						++num_ignored;
						continue;
						}
					rec1.setMateReferenceName(newName2);
					}
				if(rec1.hasAttribute(SAMTag.SA.name()))
					{
					final Object sa = rec1.getStringAttribute(SAMTag.SA.name());
					if(sa!=null && (sa instanceof String)) {
						final String tokens[]=SEMICOLON_PAT.split(String.class.cast(sa));
						final List<String> L = new ArrayList<>(tokens.length);
						for(final String semiColonStr:tokens) {
							final String commaStrs[] = COMMA_PAT.split(semiColonStr);
							final String newctg = convert(commaStrs[0]);
							if(newctg==null) continue;
							commaStrs[0]=newctg;
							L.add(String.join(",", commaStrs));
							}
						if(L.isEmpty())
							{
							rec1.setAttribute(SAMTag.SA.name(),null);
							}
						else
							{
							rec1.setAttribute(SAMTag.SA.name(),String.join(";",L));
							}
						}
					}
				
				sfw.addAlignment(rec1);
				}
			if(!this.unmappedChromosomes.isEmpty())
				{
				LOG.warning("Unmapped chromosomes: "+unmappedChromosomes);
				}
			LOG.warning("num ignored read:"+num_ignored);
			progress.close();
			return 0;
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
	

	public static void main(final String[] args)
		{
		new ConvertBamChromosomes().instanceMainWithExit(args);
		}
	}
