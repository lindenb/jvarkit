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
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.IntervalFilter;
import htsjdk.samtools.util.AbstractProgressLogger;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.SequenceUtil;

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

/** assert REF was declarated by user */
protected Path getRequiredReferencePath() {
	if(this.faidxPath==null) {
		throw new IllegalStateException("Reference was not specified. " + INDEXED_FASTA_REFERENCE_DESCRIPTION);
		}
	return this.faidxPath;
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

/** create input SAMRecord iterator for one sample*/
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

protected void scanIterator(final CloseableIterator<SAMRecord> iter,final SAMFileWriter sfw) {
final Function<SAMRecord,List<SAMRecord>> modifier = createSAMRecordFunction();
while(iter.hasNext()) {
	final SAMRecord rec = iter.next();
	for(final SAMRecord R :modifier.apply(rec)) {
		sfw.addAlignment(R);
		}
	}
}



@Override
public int doWork(final List<String> args0) {
	final Map<SamReader,CloseableIterator<SAMRecord>> sam2iterator = new HashMap<>();
	final CloseableIterator<SAMRecord> mainIterator;
	final SAMFileHeader mainHeader;
	final List<String> inputs = new ArrayList<>(args0.size());
		
		/* parse input */
		try {
			// unroll list
			if(args0.size()==1 && args0.get(0).endsWith(".list")) {
				final Path p1 = Paths.get(args0.get(0));
				IOUtil.assertFileIsReadable(p1);
				inputs.addAll(Files.
						lines(p1).
						filter(S->!S.startsWith("#")).
						filter(S->!StringUtils.isBlank(S)).
						collect(Collectors.toList())
						);
				if(inputs.isEmpty()) {
					LOG.error("List is empty " + p1);
					return -1;
					}
				}
			else
				{
				inputs.addAll(args0);
				}
		} catch(final IOException err) {
			LOG.error(err);
			return -1;
			}
		
		/* check input is not output */
		if(!inputs.isEmpty() &&
			this.outputFile!=null &&
			inputs.stream().
				filter(F->!IOUtil.isUrl(F)).
				map(F->Paths.get(F)).
				anyMatch(F->F.equals(this.outputFile))) {
			LOG.error("Input file == output file : "+ inputs);
			return -1;
			}
		
		/* before SAM */
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
			if(inputs.isEmpty() || inputs.size()==1) {
				final SamReader in;
				if(inputs.isEmpty()) {
					in = srf.open(SamInputResource.of(stdin()));
					}
				else if(IOUtil.isUrl(inputs.get(0))){
					in = srf.open(SamInputResource.of(new URL(inputs.get(0))));
					}
				else 
					{
					in = srf.open(Paths.get(inputs.get(0)));
					}
				final CloseableIterator<SAMRecord> iter= openSamIterator(in);
				sam2iterator.put(in, iter);
				mainHeader = in.getFileHeader();
				mainIterator = iter;
				}
			else
				{
				for(final String bamPath: inputs) {
					final SamReader  in = srf.open(Paths.get(bamPath));
					final CloseableIterator<SAMRecord> iter= openSamIterator(in);
					sam2iterator.put(in, iter);
					}
				final Set<SAMFileHeader.SortOrder> all_sort_orders = sam2iterator.
						keySet().
						stream().
						map(SR->SR.getFileHeader()).
						map(H->H.getSortOrder()).
						collect(Collectors.toSet());
				
				if(all_sort_orders.size()!=1) {
					LOG.error("Heterogenous sort order in input bams : " + all_sort_orders);
					return -1;
					}
				final SAMFileHeader.SortOrder sortOrder = all_sort_orders.iterator().next();
				
				final SamFileHeaderMerger headerMerger  = new SamFileHeaderMerger(
						sortOrder,
						sam2iterator.keySet().stream().map(SR->SR.getFileHeader()).collect(Collectors.toList()),
						false);
				mainHeader  = headerMerger.getMergedHeader();
				@SuppressWarnings("resource")
				final MergingSamRecordIterator iter = new MergingSamRecordIterator( headerMerger,sam2iterator, false);
				mainIterator = iter;
				}
			if(this.faidxPath!=null) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidxPath);
				SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(mainHeader));
			}
			
			
			try(SAMFileWriter sfw = openSamFileWriter(mainHeader)) {
				final ProgressLoggerInterface progress = createProgressLogger();
				if(progress!=null) sfw.setProgressLogger( progress);
				scanIterator(mainIterator,sfw);	
				}
	
			mainIterator.close();
			for(final SamReader sr: sam2iterator.keySet()) {
				sam2iterator.get(sr).close();
				sr.close();
				}
			sam2iterator.clear();
			return 0;
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
			
			for(final SamReader sr: sam2iterator.keySet()) {
				CloserUtil.close(sam2iterator.get(sr));
				CloserUtil.close(sr);
				}
			}
		}
	}
