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
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;

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

protected Function<SAMRecord,List<SAMRecord>> createSAMRecordFunction() {
	return R->Collections.singletonList(R);
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
		if(getLogger()!=null) {
			}
		int err= 0;
		final Function<SAMRecord,List<SAMRecord>> modifier = createSAMRecordFunction();
		final SAMFileHeader headerIn = in.getFileHeader();
		try(SAMFileWriter sfw = openSamFileWriter(headerIn)) {
			try(CloseableIterator<SAMRecord> iter= in.iterator()) {
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
