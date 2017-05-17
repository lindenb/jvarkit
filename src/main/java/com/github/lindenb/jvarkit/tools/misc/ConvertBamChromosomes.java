package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

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


END_DOC

 */
@Program(
	name="bamrenamechr",description="Convert the names of the chromosomes in a BAM file",
	keywords={"sam","bam","chromosome","contig"}
	)
public class ConvertBamChromosomes
	extends Launcher
	{
	private static final Logger LOG = Logger.build(ConvertBamChromosomes.class).make();
	
	
	@Parameter(names={"-c","-convert"},description="What should I do when  a converstion is not found")
	private ContigNameConverter.OnNotFound onNotFound=ContigNameConverter.OnNotFound.RAISE_EXCEPTION;
	@Parameter(names={"-f","--mapping","-m"},description="load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+",required=true)
	private File mappingFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File output= null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	private ConvertBamChromosomes()
		{
		
		}
	
	@Override
	public int doWork(List<String> args) {
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		final Set<String> unmappedChromosomes = new HashSet<>();
		try
			{
			final ContigNameConverter customMapping = ContigNameConverter.fromFile(mappingFile);
			customMapping.setOnNotFound(this.onNotFound);

			
			sfr = super.openSamReader(oneFileOrNull(args));
			
			SAMFileHeader header1=sfr.getFileHeader();
			if(header1==null)
				{
				LOG.error("File header missing");
				return -1;
				}
			
			
			SAMFileHeader header2=header1.clone();
			
			//create new sequence dict
			final SAMSequenceDictionary dict1=header1.getSequenceDictionary();
			if(dict1==null)
				{
				LOG.error("Sequence dict missing");
				return -1;
				}
			final List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>(dict1.size());
			for(int i=0;i< dict1.size();++i)
				{
				SAMSequenceRecord ssr=dict1.getSequence(i);
				String newName=customMapping.apply(ssr.getSequenceName());
				if(newName==null)
					{
					//skip unknown chromosomes
					continue;
					}
				ssr=new SAMSequenceRecord(newName, ssr.getSequenceLength());
				ssrs.add(ssr);
				}
			header2.setSequenceDictionary(new SAMSequenceDictionary(ssrs));
			
			SAMSequenceDictionary dict2=new SAMSequenceDictionary(ssrs);
			header2.setSequenceDictionary(dict2);
			
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1);

			sfw = this.writingBamArgs.openSAMFileWriter(output, header2, true);
			
			
			long num_ignored=0L;
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec1=iter.next();
				progress.watch(rec1);
				String newName1=null;
				String newName2=null;
				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					newName1=customMapping.apply(rec1.getReferenceName());
					}
				if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
					{
					newName2=customMapping.apply(rec1.getMateReferenceName());
					}
				rec1.setHeader(header2);

				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					if(newName1==null)
						{
						++num_ignored;
						continue;
						}
					rec1.setReferenceName(newName1);
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
				sfw.addAlignment(rec1);
				}
			if(!unmappedChromosomes.isEmpty())
				{
				LOG.warning("Unmapped chromosomes: "+unmappedChromosomes);
				}
			LOG.warning("num ignored read:"+num_ignored);
			return 0;
			}
		catch(Exception err)
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
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new ConvertBamChromosomes().instanceMainWithExit(args);
		}
	}
