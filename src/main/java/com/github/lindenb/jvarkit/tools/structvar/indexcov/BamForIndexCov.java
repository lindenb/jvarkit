package com.github.lindenb.jvarkit.tools.structvar.indexcov;

import java.io.StringWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileSource;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StopWatch;

public class BamForIndexCov extends Launcher {
	private static final Logger LOG = Logger.build(BamForIndexCov.class).make();
	@Parameter(names={"-o","--output"},description="Output directory",required=true)
	private Path outputDir = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-p","--prefix"},description="File prefix")
	private String prefix = "bam4indexcov.";
	@Parameter(names={"-Q","--mapq"},description="min mapq")
	private int mapq = 10;
	@Parameter(names={"-x","--exclude-chromosomes"},description="regex of chromosomes to skip")
	private String excludeChrRegex ="(hs37d5|chrM|chrMT|MT)";

	
	private class BamWriterForIndexCov implements AutoCloseable {
	    private final BinaryCodec outputBinaryCodec;
	    private BAMRecordCodec bamRecordCodec = null;
	    private final BlockCompressedOutputStream blockCompressedOutputStream;
	    private BAMIndexer bamIndexer = null;
	    private final SAMFileHeader header;
	    private final Path bamPath;
	    BamWriterForIndexCov(final Path bamPath,final Path baiPath,final SAMFileHeader header) {
	        final String bamFilename = bamPath.toString();
	        this.header = header;
	        this.bamPath = bamPath;
	    	blockCompressedOutputStream = new BlockCompressedOutputStream(new NullOuputStream(),bamPath,0);
	        outputBinaryCodec = new BinaryCodec(blockCompressedOutputStream);
	        outputBinaryCodec.setOutputFileName(bamFilename);
	        bamRecordCodec = new BAMRecordCodec(header);
            bamRecordCodec.setOutputStream(outputBinaryCodec.getOutputStream(),bamFilename);
           
            bamIndexer =  new BAMIndexer(baiPath,header);
            
            //write header
            outputBinaryCodec.writeBytes("BAI\1".getBytes());

            // calculate and write the length of the SAM file header text and the header text
            final Writer stringWriter = new StringWriter();
            new SAMTextHeaderCodec().encode(stringWriter, header, true);

            outputBinaryCodec.writeString(stringWriter.toString(), true, false);

            // write the sequences binarily.  This is redundant with the text header
            outputBinaryCodec.writeInt(header.getSequenceDictionary().size());
            for (final SAMSequenceRecord sequenceRecord: header.getSequenceDictionary().getSequences()) {
                outputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
                outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
            	}

	    	}

	    
	    void add(final SAMRecord alignment) {
            final long startOffset = blockCompressedOutputStream.getFilePointer();
            bamRecordCodec.encode(alignment);
            final long stopOffset = blockCompressedOutputStream.getFilePointer();
            // set the alignment's SourceInfo and then prepare its index information
            alignment.setFileSource(new SAMFileSource(null, new BAMFileSpan(new Chunk(startOffset, stopOffset))));
            bamIndexer.processAlignment(alignment);
    		}

	    @Override
	    public void close() {
	        outputBinaryCodec.close();
	        
			final SAMFileWriterFactory sfw = new SAMFileWriterFactory().
					setCreateIndex(false).
					setCreateMd5File(false).
					setCompressionLevel(0);
			try(SAMFileWriter sw=sfw.makeBAMWriter(header, true, this.bamPath)) {
				//do nothing, just write header
				}
	        bamIndexer.finish();
	    	}
	    }



	
	
	@Override
	public int doWork(List<String> args)
		{
		try {
			final List<String> bams = IOUtils.unrollStrings2018(args);
			if(bams.isEmpty()) {
				LOG.error("input is missing");
				return -1;
				}
			IOUtil.assertFileIsReadable(this.faidx);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidx);
			final Set<Integer> tid_to_skips =  dict.getSequences().stream().
					filter(S->S.getSequenceName().matches(this.excludeChrRegex)).
					map(SSR->SSR.getSequenceIndex()).
					collect(Collectors.toSet());
				
			
			
			
			IOUtil.assertDirectoryIsWritable(this.outputDir);
			final SamReaderFactory srf = super.createSamReaderFactory().
					referenceSequence(this.faidx);
			for(final String bamPathStr:bams) {
				final StopWatch stopWatch = new  StopWatch();
				stopWatch.start();
				LOG.error("Processing:" + bamPathStr);
				try(SamReader sr = srf.open(SamInputResource.of(bamPathStr))) {
					final SAMFileHeader header = sr.getFileHeader();
					if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
						LOG.error("Error, file "+bamPathStr+" is not sorted on coordinate.");
						return -1;
						}	
					final String sn = header.
						getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtils.isBlank(S)).
						findFirst().
						orElseThrow(()->new SAMException("Cannot find sample name in "+bamPathStr));
					
					final Path outBamPath = this.outputDir.resolve(this.prefix+sn+FileExtensions.BAM);
					final Path outBaiPath = this.outputDir.resolve(this.prefix+sn+FileExtensions.BAM+FileExtensions.BAI_INDEX);
					if(Files.exists(outBamPath)) {
						LOG.error("Error, file "+outBamPath+". Already exists.");
						return -1;
						}
					if(Files.exists(outBaiPath)) {
						LOG.error("Error, file "+outBaiPath+". Already exists.");
						return -1;
						}
					
					
					
					try(BamWriterForIndexCov sw= new BamWriterForIndexCov(outBamPath,outBaiPath,header)) {
						try(CloseableIterator<SAMRecord> iter = sr.iterator()) {
							while(iter.hasNext()) {
								final SAMRecord rec = iter.next();
								if(!SAMRecordDefaultFilter.accept(rec,this.mapq)) continue;
								if(tid_to_skips.contains(rec.getReferenceIndex())) continue;
								sw.add(rec);
								}
							}
						}
						
					}
				stopWatch.stop();
				LOG.error("That took:" + StringUtils.niceDuration(stopWatch.getElapsedTime()) );
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(String[] args)
		{
		new BamForIndexCov().instanceMainWithExit(args);
		}

	}
