package com.github.lindenb.jvarkit.tools.structvar;

import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.DiskBasedBAMFileIndex;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC

## About

exploring  Brent Pedersen 's idea of use of the BAI indexes. 

see [https://github.com/brentp/goleft/tree/master/indexcov](https://github.com/brentp/goleft/tree/master/indexcov)

> Quickly estimate coverage from a whole-genome bam or cram index. A bam index has 16KB resolution so that's what this gives, but it provides what appears to be a high-quality coverage estimate in seconds per genome.


END_DOC
 */
@Program(
		name="baicov",
		description="exploring  Brent Pedersen 's idea of use of the BAI indexes",
		keywords={"vcf","cnv","bam"},
		creationDate="20200529",
		modificationDate="20200529",
		generate_doc=false
		)
public class BaiCov extends Launcher {
	private static final Logger LOG=Logger.build(BaiCov.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	
	private static int virtualOffset(long offset) {
		return BlockCompressedFilePointerUtil.getBlockOffset(offset);
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(this.faidx);
			for(Path path: IOUtils.unrollPaths(args)) {
				String sampleName = null;
				if(path.getFileName().toString().endsWith(FileExtensions.BAM)) {
					final Path path2  = SamFiles.findIndex(path);
					if(path2==null) {
						LOG.error("Cannot find "+ FileExtensions.BAI_INDEX+" associated to "+path);
						return -1;
						}
					try(SamReader sr = SamReaderFactory.makeDefault().referenceSequence(this.faidx).open(path)) {
						final SAMFileHeader header= sr.getFileHeader();
						SequenceUtil.assertSequenceDictionariesEqual(dict, header.getSequenceDictionary());
						sampleName = header.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse(null);
						}
					path = path2;
					}
				if(!path.getFileName().toString().endsWith(FileExtensions.BAI_INDEX)) {
					LOG.error("Not a path with suffix "+ FileExtensions.BAI_INDEX+":" + path);
					return -1;
					}
				if(StringUtils.isBlank(sampleName)) {
					sampleName = IOUtils.getFilenameWithoutCommonSuffixes(path);
					}
				try(DiskBasedBAMFileIndex index = new DiskBasedBAMFileIndex(path.toFile(),dict)) {
					for(int tid=0;tid< dict.size();tid++) {
						final SAMSequenceRecord ssr= dict.getSequence(tid);
						final BAMIndexMetaData meta = index.getMetaData(tid);
						final BAMFileSpan span = index.getSpanOverlapping(tid,1, ssr.getEnd());
						System.err.print(path+" "+sampleName+" "+ ssr.getSequenceName()+"\t"+meta.getAlignedRecordCount()+"\t"+meta.getUnalignedRecordCount());
						for(Chunk chunk: span.getChunks()) {
							System.err.print(" "+virtualOffset(chunk.getChunkStart())+","+virtualOffset(chunk.getChunkEnd()));
							}
						System.err.println();
						}
					}
				}
				        
			return 0;
			} 
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	
	
	public static void main(String[] args) {
		new BaiCov().instanceMainWithExit(args);
	}
	
}
