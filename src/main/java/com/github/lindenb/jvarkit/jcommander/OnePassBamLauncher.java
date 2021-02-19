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
package com.github.lindenb.jvarkit.jcommander;

import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.IntervalFilter;
import htsjdk.samtools.util.AbstractProgressLogger;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.ProgressLoggerInterface;

public abstract class OnePassBamLauncher extends Launcher {
private static final Logger LOG = Logger.build(OnePassBamLauncher.class).make();
@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
protected Path outputFile=null;
@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
protected Path faidxPath =null;
@Parameter(names={"--validation-stringency"},description="SAM Reader Validation Stringency")
protected ValidationStringency validationStringency = ValidationStringency.LENIENT;
@ParametersDelegate
protected WritingBamArgs writingBamArgs = new WritingBamArgs();
@Parameter(names={"--regions"},description="Limit analysis to this interval. "+ IntervalListProvider.OPT_DESC,splitter=NoSplitter.class,converter=IntervalListProvider.StringConverter.class)
protected IntervalListProvider regionFiles = null;


protected Logger getLogger() {
	return null;
	}

private void deleteOutputOnError() {
	if(this.outputFile==null) return;
	try {
		if(Files.deleteIfExists(this.outputFile)) {
			LOG.warning("The following file was deleted: "+this.outputFile);
			}
	} catch(final IOException err) {
		//ignore
	}
}
/** initialize things before opening the BAM */
protected int beforeSam() {
	return 0;
	}

/** initialize things after closing the SAM */
protected void afterSam() {
	}


/** create a new SamReaderFactory */
@Override
protected SamReaderFactory createSamReaderFactory()
	{
	return  SamReaderFactory.
			makeDefault().
			validationStringency(this.validationStringency).
			referenceSequence(this.faidxPath);
	}

protected SAMFileHeader createOutputHeader(final SAMFileHeader headerIn) {
	SAMFileHeader headerOut = headerIn.clone();
	JVarkitVersion.getInstance().addMetaData(this, headerOut);
	return headerOut;
	}

protected boolean isPresorted() {
	return true;
}

protected SAMFileWriter openSamFileWriter(final SAMFileHeader headerIn) {
	final SAMFileHeader headerOut = createOutputHeader(headerIn);
	return this.writingBamArgs.
			setReferencePath(this.faidxPath).
			openSamWriter(this.outputFile, headerOut, isPresorted());
}

/** the very basic filter , for one read return all the reads that should be produced by one record. Outut can be empty */
protected Function<SAMRecord,List<SAMRecord>> createSAMRecordFunction() {
	return R->Collections.singletonList(R);
}

/** create input SAMRecord iterator */
protected CloseableIterator<SAMRecord> openSamIterator(final SamReader sr) {
	
	if(this.regionFiles!=null) {
		this.regionFiles.dictionary(SequenceDictionaryUtils.extractRequired(sr.getFileHeader()));

		if(!sr.hasIndex()) {
			final List<Interval> L = this.regionFiles.stream().map(x->new Interval(x)).collect(Collectors.toList());
			final SAMRecordIterator st0 = sr.iterator();
			return new FilteringSamIterator(st0,new IntervalFilter(L, sr.getFileHeader()));
			}
		else
			{
			return sr.query(this.regionFiles.optimizedQueryIntervals(), false);
			}
		}
	return sr.iterator();
	}

/** create a progress logger for the BAM writer. Result may be null */
protected ProgressLoggerInterface createProgressLogger() {
	final Logger log = getLogger();
	if(log==null) return null;
	return new AbstractProgressLogger(getProgramName(),"Writer",1_000_000)
		{
		@Override
		protected void log(String... message)
			{
			log.info(String.join(" ", message));
			}
		};
}

@Override
public int doWork(final List<String> args) {
	SamReader in = null;
	final String input = super.oneFileOrNull(args);
	if(input!=null &&
		this.outputFile!=null &&
		!IOUtil.isUrl(input) && Paths.get(input).equals(this.outputFile)) {
		LOG.error("Input == output : "+ input);
		return -1;
		}
	
	try {
		if(beforeSam()!=0) {
			LOG.error("initialization failed");
			return -1;
			}
		}
	catch (final Throwable err) {
		LOG.error(err);
		return -1;
		}
	
	
	
	try {
		final SamReaderFactory srf = createSamReaderFactory();
		if(input==null) {
			in = srf.open(SamInputResource.of(stdin()));
			}
		else if(IOUtil.isUrl(input)){
			in = srf.open(SamInputResource.of(new URL(input)));
			}
		else
			{
			in = srf.open(Paths.get(input));
			}
		
		int err= 0;
		final Function<SAMRecord,List<SAMRecord>> modifier = createSAMRecordFunction();
		final SAMFileHeader headerIn = in.getFileHeader();
		try(SAMFileWriter sfw = openSamFileWriter(headerIn)) {
			final ProgressLoggerInterface progress = createProgressLogger();
			if(progress!=null) sfw.setProgressLogger( progress);
			try(CloseableIterator<SAMRecord> iter= openSamIterator(in)) {
				while(iter.hasNext()) {
				final SAMRecord rec = iter.next();
				for(final SAMRecord R :modifier.apply(rec)) {
					sfw.addAlignment(R);
					}
				}
			}
		}

		
		in.close();
		in = null;
		if(err!=0) deleteOutputOnError();
		return err;
		}
	catch (final Throwable err) {
		LOG.error(err);
		deleteOutputOnError();
		return -1;
		}
	finally
		{
		try {
			afterSam();
			}
		catch (final Throwable err) {
			LOG.error(err);
			}	
		try {if(in!=null) in.close(); } catch (final Throwable err) {LOG.error(err);}
		}
	}
}
