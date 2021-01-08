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
package com.github.lindenb.jvarkit.tools.minibam;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.zip.Deflater;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

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
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StopWatch;

/**
 BEGIN_DOC
 

# Motivation

Bams are too bigs and my users often ask to visualize a small region of a set of bam

# Input

 input is a set of bam files or a file with the suffix '.list' containing one path to a bam per line.
 
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
		description="Creates an archive of small bams with only a few regions.",
		keywords={"bam","sam"},
		creationDate="20190410",
		modificationDate="20191010"
		)
public class MakeMiniBam extends Launcher {
	private static final Logger LOG = Logger.build(MakeMiniBam.class).make();

	
	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path outputFile = null;
	@Parameter(names={"-x","--extend"},description="Extend the positions by 'x' bases. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extend = 5_000;
	@Parameter(names={"-T","--tmp"},description="Tmp working directory")
	private Path tmpDir = IOUtils.getDefaultTempDir();
	@Parameter(names={"--prefix"},description="File prefix in the archive. Special value 'now' or empty string will be replaced by the current date")
	private String filePrefix="miniBam.";
	@Parameter(names={"-B","--bed","-p","--pos","-V","--variant","--vcf"},description=IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,splitter=NoSplitter.class,required=true)
	private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();
	@Parameter(names={"-R","--reference"},description="Optional Reference file for CRAM files. "+ INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path referencePath = null;
	@Parameter(names={"--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter samRecordFilter = SamRecordFilterFactory.ACCEPT_ALL;
	@Parameter(names={"-b","--bounds","--edge"},description="[20190427] If `b` is greater than 0 and the user interval has a length greater than `b` then consider the edges of the object as two positions. the idea is to just save the boundaries of a large deletion. " + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int bound_edge = -1;
	@Parameter(names={"-C","--comment"},description="[20190427]Add a file '*.md' with this comment.")
	private String userComment="";
	@Parameter(names={"--bnd"},description="[20190427]When reading VCF file, don't get the mate position for the structural BND variants.")
	private boolean disable_sv_bnd = false;
	@Parameter(names={"--no-samples"},description="[20191129]Allow no sample/ no read group : use fileame")
	private boolean allow_no_sample = false;

	
	@Override
	public int doWork(final List<String> args) {
		ArchiveFactory archive = null;
		int id_generator=0;
		final Set<String> outputFileNames = new HashSet<>();
		try {
			
			
			if (StringUtils.isBlank(this.filePrefix) || this.filePrefix.equals("now")) {
				final SimpleDateFormat simpleDateFormat = new SimpleDateFormat("yyyyMMdd");
				this.filePrefix = simpleDateFormat.format(new Date()) + ".";
			}
			if(!this.filePrefix.endsWith(".")) {
				this.filePrefix+=".";
				}
			
			IOUtil.assertDirectoryIsWritable(tmpDir);
			final List<Path> bamFiles = IOUtils.unrollPaths(args);
			if(bamFiles.isEmpty()) {
				LOG.error("no bam file defined");
				return -1;
				}
			
			final List<Locatable> locatables = this.intervalListProvider.
					enableBreakEndInterval(!disable_sv_bnd).
					enableSinglePoint().
					stream().
					collect(Collectors.toList())
					;
			
			if(locatables.isEmpty()) {
				LOG.error("no position defined");
				return -1;
				}
			
			final SAMFileWriterFactory swf = new SAMFileWriterFactory();
			swf.setCompressionLevel(9);
			swf.setCreateIndex(true);
			
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.referencePath!=null) srf.referenceSequence(this.referencePath);
			
			archive = ArchiveFactory.open(this.outputFile);
			archive.setCompressionLevel(Deflater.NO_COMPRESSION);
			for(final Path bamFile:bamFiles) {
				LOG.info(bamFile.toString());
				final StopWatch stopWatch = new StopWatch();
				stopWatch.start();
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
				
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				final Optional<String> dictLabel = SequenceDictionaryUtils.getBuildName(dict);
				final String labelSuffix = 
						(dictLabel.isPresent()? "." + dictLabel.get():"") + 
						(locatables.size()==1?"."+locatables.get(0).getContig()+"_"+locatables.get(0).getStart()+(locatables.get(0).getStart()==locatables.get(0).getEnd()?"":"_"+locatables.get(0).getEnd()):"")
						;
						
				final ContigNameConverter ctgConvert = ContigNameConverter.fromOneDictionary(dict);
				
				final QueryInterval array[]= locatables.
						stream().
						flatMap(loc->{
							 if(this.bound_edge<1 || loc.getLengthOnReference() <= this.bound_edge) {
								 return Collections.singletonList(loc).stream();
								 }
							return Arrays.asList(
									(Locatable)new SimpleInterval(loc.getContig(),loc.getStart(),loc.getStart()),
									(Locatable)new SimpleInterval(loc.getContig(),loc.getEnd(),loc.getEnd())
							        ).stream();
							}).
						map(LOC->{
							final String contig  = ctgConvert.apply(LOC.getContig());
							if(StringUtils.isBlank(contig)) {
								LOG.warn("Cannot find "+LOC.getContig()+" in "+bamFile);
								return null ;
								}
							final SAMSequenceRecord ssr = dict.getSequence(contig);
							
							if(LOC.getStart()> ssr.getSequenceLength()) {
								LOG.warn("pos "+LOC+" is greater than chromosome size "+ssr.getSequenceLength()+" in "+bamFile);
								return null;
								}
							return new QueryInterval(
									ssr.getSequenceIndex(), 
									Math.max(1, LOC.getStart()-extend),
									Math.min(ssr.getSequenceLength(), LOC.getEnd()+extend)
									);
							}).
						filter(Q->Q!=null).
						collect(HtsCollectors.optimizedQueryIntervals());
				
				
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
							if(this.samRecordFilter.filterOut(rec)) continue;
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
	            
				
				final String sampleName1 = header.getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtils.isBlank(S)).
						findFirst().
						orElse(null);
				final String sampleName;
				if(!StringUtils.isBlank(sampleName1))
					{
					sampleName = sampleName1;
					}
				else if(this.allow_no_sample) {
					sampleName = IOUtils.getFilenameWithoutSuffix(bamFile,1);
					LOG.warn("No Read group in "+bamFile+" using filename : "+sampleName);
					}
				else
					{
					throw new IllegalArgumentException("No Sample found in "+bamFile+". Use --no-samples option ?");
					}
				
				String filename=this.filePrefix + sampleName + labelSuffix;
				while(outputFileNames.contains(filename)) 	{
					filename=this.filePrefix+sampleName+"."+ (id_generator++) + labelSuffix;
					}
				outputFileNames.add(filename);
				archive.copyTo(tmpBam, filename + FileExtensions.BAM);
				archive.copyTo(tmpBai, filename + IOUtils.getFileSuffix(tmpBai));
				stopWatch.stop();
				LOG.info("Added "+ StringUtils.niceFileSize(Files.size(tmpBam))+"(bam) and "+StringUtils.niceFileSize(Files.size(tmpBai))+" (Bai). "+StringUtils.niceDuration(stopWatch.getElapsedTime()));
				
				Files.deleteIfExists(tmpBam);
				Files.deleteIfExists(tmpBai);
				}
			
			if(!StringUtils.isBlank(this.userComment)) {
				try(final PrintWriter pw = archive.openWriter(this.filePrefix+(this.filePrefix.endsWith(".")?"":".")+"README.md")) {
					pw.append(this.userComment);
					pw.println();
					pw.println("## BAMS");
					pw.println();
					for(final Path bamFile:bamFiles) pw.println("  * "+bamFile);
					pw.println();
					pw.println("## Date");
					pw.println();
					pw.println( new SimpleDateFormat("yyyyMMdd").format(new Date()));
					pw.flush();
					}
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
