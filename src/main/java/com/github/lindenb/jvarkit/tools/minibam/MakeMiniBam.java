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
package com.github.lindenb.jvarkit.tools.minibam;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.vcf.VCFFileReader;
/**
 BEGIN_DOC
 
# Motivation

Bams are too bigs and my users often ask to visualize a small region of a set of bam
 
# Example
 
```
$  find src/test/resources/ -name "S*.bam" > bams.list
$   java -jar dist/mkminibam.jar -p "RF01:100" -o out.zip bams.list 
[INFO][MakeMiniBam]src/test/resources/S5.bam
[INFO][MakeMiniBam]src/test/resources/S2.bam
[INFO][MakeMiniBam]src/test/resources/S4.bam
[INFO][MakeMiniBam]src/test/resources/S3.bam
[INFO][MakeMiniBam]src/test/resources/S1.bam

$ unzip -t out.zip 

Archive:  out.zip
    testing: miniBam.S5.bam           OK
    testing: miniBam.S5.bai           OK
    testing: miniBam.S2.bam           OK
    testing: miniBam.S2.bai           OK
    testing: miniBam.S4.bam           OK
    testing: miniBam.S4.bai           OK
    testing: miniBam.S3.bam           OK
    testing: miniBam.S3.bai           OK
    testing: miniBam.S1.bam           OK
    testing: miniBam.S1.bai           OK
No errors detected in compressed data of out.zip.
```
 
 
 END_DOC
 
 */
@Program(
		name="mkminibam",
		description="Creates an archive of small bams with only a few region",
		keywords={"bam","sam"},
		creationDate="20190410",
		modificationDate="20190410"
		)
public class MakeMiniBam extends Launcher {
	private static final Logger LOG = Logger.build(MakeMiniBam.class).make();

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-x","--extend"},description="Extend the positions by 'x' bases. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extend = 5_000;
	@Parameter(names={"-p","--pos"},description="Add this position 'chrom:pos'")
	private Set<String> posStrSet = new HashSet<>();
	@Parameter(names={"-T","--tmp"},description="Tmp working directory")
	private Path tmpDir = IOUtils.getDefaultTempDir();
	@Parameter(names={"--prefix"},description="File prefix in the archive")
	private String filePrefix="miniBam.";
	@Parameter(names={"-V","--variant"},description="Use the intervals from this VCF file.")
	private Path vcfInput=null;
	@Parameter(names={"-B","--bed"},description="Use the intervals from this BED file.")
	private Path bedInput=null;

	
	
	@Override
	public int doWork(final List<String> args) {
		ArchiveFactory archive = null;
		int id_generator=0;
		final Set<String> outputFileNames = new HashSet<>();
		try {
			IOUtil.assertDirectoryIsWritable(tmpDir);
			final List<Path> bamFiles = IOUtils.unrollPaths(args);
			if(bamFiles.isEmpty()) {
				LOG.error("no bam file defined");
				return -1;
				}
			if(posStrSet.isEmpty() && this.vcfInput==null && this.bedInput==null) {
				LOG.error("no position defined");
				return -1;
				}
			final SAMFileWriterFactory swf = new SAMFileWriterFactory();
			swf.setCompressionLevel(9);
			swf.setCreateIndex(true);
			
			final SamReaderFactory srf = super.createSamReaderFactory();
			
			archive = ArchiveFactory.open(this.outputFile);
			for(final Path bamFile:bamFiles) {
				LOG.info(bamFile.toString());
				final SamReader sr = srf.open(bamFile);
				if(!sr.hasIndex()) {
					sr.close();
					LOG.error("file "+bamFile+" is not indexed.");
					return -1;
					}
				final SAMFileHeader header = sr.getFileHeader();
				if(!header.getSortOrder().equals(SortOrder.coordinate)) {
					sr.close();
					LOG.error("file "+bamFile+" is not sorted on coordinate.");
					return -1;
					}
				
				final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(header);
				final ContigNameConverter ctgConvert = ContigNameConverter.fromOneDictionary(dict);
				
				final List<QueryInterval> queryIntervals = new ArrayList<>(this.posStrSet.size());
				
				final Consumer<Locatable> locatableConsummer = (LOC)->{
					final String contig  = ctgConvert.apply(LOC.getContig());
					if(StringUtils.isBlank(contig)) {
						LOG.warn("Cannot find "+LOC.getContig()+" in "+bamFile);
						return;
						}
					final SAMSequenceRecord ssr = dict.getSequence(contig);
					
					if(LOC.getStart()> ssr.getSequenceLength()) {
						LOG.warn("pos "+LOC+" is greater than chromosome size "+ssr.getSequenceLength()+" in "+bamFile);
						return;
						}
					final QueryInterval queryInterval = new QueryInterval(
							ssr.getSequenceIndex(), 
							Math.max(1, LOC.getStart()-extend),
							Math.min(ssr.getSequenceLength(), LOC.getEnd()+extend)
							);
					queryIntervals.add(queryInterval);
					};
				
				for(final String posStr:this.posStrSet) {
					int colon = posStr.indexOf(":");
					if(colon <=0 ) throw new IllegalArgumentException("cannot find colon in "+posStr);
					final String contig1 = posStr.substring(0,colon);
					int pos = Integer.parseInt(posStr.substring(colon+1));
					locatableConsummer.accept(new Interval(contig1,pos,pos));
					}
				if(this.vcfInput!=null) {
					try(VCFFileReader vr = new VCFFileReader(this.vcfInput,false)) {
						vr.forEach(locatableConsummer);
						}
					}
				if(this.bedInput!=null) {
					final BedLineCodec codec= new BedLineCodec();
					try(BufferedReader br = IOUtils.openPathForBufferedReading(this.bedInput)) {
						br.lines().
							filter(L->!BedLine.isBedHeader(L)).
							map(L->codec.decode(L)).
							filter(B->B!=null).
							map(B->B.toInterval()).
							forEach(locatableConsummer);
						}
					}
				
				
				
				final QueryInterval array[]= QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));
				
				
				final Path tmpBam = Files.createTempFile(this.tmpDir,"tmp.", ".bam");
				
				final SAMFileHeader header2= header.clone();
				final SAMProgramRecord prg = header2.createProgramRecord();
				prg.setProgramName(this.getProgramName());
				prg.setProgramVersion(this.getGitHash());
				prg.setCommandLine(this.getProgramCommandLine());
				
				JVarkitVersion.getInstance().addMetaData(this, header2);
				header2.addComment("MiniBam : Bam was "+bamFile);
				
				try(SAMFileWriter sfw = swf.makeBAMWriter(header2, true,tmpBam)) {
					
				if(array.length>0) {
					try(CloseableIterator<SAMRecord> ssr= sr.query(array, false)) {
						while(ssr.hasNext())
							{
							final SAMRecord rec = ssr.next();
							rec.setAttribute(SAMTag.PG.name(), prg.getId());
							sfw.addAlignment(rec);
							}
						}
					}
				}
				
	            final Path tmpBai = SamFiles.findIndex(tmpBam);
				if(!Files.exists(tmpBai))
					{
					LOG.error("Cannot find tmp bai Index for "+bamFile+" "+tmpBam);
					return -1;
					}
	            
				
				final String sampleName = header.getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtils.isBlank(S)).
						findFirst().
						orElseThrow(()->new IllegalArgumentException("No Sample found in "+bamFile));
				String filename=this.filePrefix+sampleName;
				while(outputFileNames.contains(filename)) 	{
					filename=this.filePrefix+sampleName+"."+ (id_generator++);
					}
				outputFileNames.add(filename);
				archive.copyTo(tmpBam, filename + BamFileIoUtils.BAM_FILE_EXTENSION);
				archive.copyTo(tmpBai, filename + IOUtils.getFileSuffix(tmpBai));

				
				Files.deleteIfExists(tmpBam);
				Files.deleteIfExists(tmpBai);
				}
			
			archive.close();
			archive=null;
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(archive);
			}
		
		}
	
	public static void main(final String[] args) {
		new MakeMiniBam().instanceMainWithExit(args);
	}

}
