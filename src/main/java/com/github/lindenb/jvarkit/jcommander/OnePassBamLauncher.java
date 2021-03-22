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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.AbstractProgressLogger;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;

public abstract class OnePassBamLauncher extends MultiBamLauncher {
private static final Logger LOG = Logger.build(OnePassBamLauncher.class).make();
@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
protected Path outputFile=null;
@ParametersDelegate
protected WritingBamArgs writingBamArgs = new WritingBamArgs();

@Override
protected Logger getLogger() {
	return LOG;
	}

@Override
protected void deleteOutputOnError() {
	if(this.outputFile==null) return;
	try {
		if(Files.deleteIfExists(this.outputFile)) {
			LOG.warning("The following file was deleted: "+this.outputFile);
			}
	} catch(final IOException err) {
		//ignore
	}
}

@Override
/** check input != output */
protected int validateInputsPath(final List<String> inputs) {
	if(!inputs.isEmpty() &&
		this.outputFile!=null &&
		inputs.stream().
			filter(F->!IOUtil.isUrl(F)).
			map(F->Paths.get(F)).
			anyMatch(F->F.equals(this.outputFile))) {
		LOG.error("Input file == output file : "+ inputs);
		return -1;
		}
	return 0;
	}



/** create SAM output header */
protected SAMFileHeader createOutputHeader(final SAMFileHeader headerIn) {
	final SAMFileHeader headerOut = headerIn.clone();
	JVarkitVersion.getInstance().addMetaData(this, headerOut);
	return headerOut;
	}

/** is output presorted */
protected boolean isPresorted() {
	return true;
}

protected SAMFileWriter openSamFileWriter(final SAMFileHeader headerIn) {
	final SAMFileHeader headerOut = createOutputHeader(headerIn);
	return this.writingBamArgs.
			setReferencePath(super.faidxPath).
			openSamWriter(this.outputFile, headerOut, isPresorted());
}

/** the very basic filter , for one read return all the reads that should be produced by one record. Outut can be empty */
protected Function<SAMRecord,List<SAMRecord>> createSAMRecordFunction() {
	return R->Collections.singletonList(R);
}

/** create a new program SAM record in the header */
protected SAMProgramRecord createProgramRecord(final SAMFileHeader header2) {
	final SAMProgramRecord samProgramRecord = header2.createProgramRecord();
	samProgramRecord.setProgramName(this.getProgramName());
	samProgramRecord.setProgramVersion(this.getGitHash());
	samProgramRecord.setCommandLine(getProgramCommandLine().replace('\t', ' '));
	return samProgramRecord;
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

protected void scanIterator(final SAMFileHeader headerIn,final CloseableIterator<SAMRecord> iter,final SAMFileWriter sfw) {
final Function<SAMRecord,List<SAMRecord>> modifier = createSAMRecordFunction();
while(iter.hasNext()) {
	final SAMRecord rec = iter.next();
	for(final SAMRecord R :modifier.apply(rec)) {
		sfw.addAlignment(R);
		}
	}
}

@Override
protected int processInput(final SAMFileHeader headerIn, final CloseableIterator<SAMRecord> iter) {
	try(SAMFileWriter sfw = openSamFileWriter(headerIn)) {
		final ProgressLoggerInterface progress = createProgressLogger();
		if(progress!=null) sfw.setProgressLogger( progress);
		scanIterator(headerIn,iter,sfw);	
		}
	return 0;
	}

}
