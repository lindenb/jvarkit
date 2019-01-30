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
package com.github.lindenb.jvarkit.tools.onekgenomes;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
/**
BEGIN_DOC

## Example

```
$ $ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/current.tree"  | cut -f 1 | grep 'bam$' | \
	grep low_coverage   | grep -F '.mapped' | sed 's%^%http://ftp.1000genomes.ebi.ac.uk/vol1/%' | head -n 3 > jeter.list



```

END_DOC
*/
@Program(name="onethousandbams",
description="Scan BAMs from 1000G. Remote fail safe (I hope)",
keywords={"sam","bam","1000genomes"},
creationDate="2019-01-29",
generate_doc=false
)
public class OneThousandBam extends Launcher{
	private static final Logger LOG = Logger.build(OneThousandBam.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-cache","-C"},description="bai cache directory.",required=true)
	private Path baiCacheDir = null;
	@Parameter(names={"-tmp","-T"},description="tmp directory to store the tmp bams")
	private Path tmxpDir = IOUtils.getDefaultTempDir();
	@Parameter(names={"-bed","-B"},description="bed file.")
	private Path bedFile = null;
	@Parameter(names={"--transient"},description="delete BAI Files, in the cache, on exit.")
	private boolean deleteBaiFilesOnExit = false;
	@Parameter(names={"--cache-and-exit"},description="Exit with success once the bai have been downloaded in the cache. (for workflow managers working in parallel)")
	private boolean exitAfterBaiDownload = false;

	private class Sample1KG
		{
		final String url;
		String sampleName=null;
		SAMReadGroupRecord readGroupRecord = null;
		Sample1KG(final String url) {
			this.url = url;
			if(!IOUtil.isUrl(url)) throw new IllegalArgumentException("Not a url : "+url);
			}
		
		public URL getBaiUrl() {
			try {
				return new URL(this.url+".bai");
			} catch (MalformedURLException e) {
				throw new IllegalArgumentException("Not a url : ",e);
				}
			}
		
		private Path getBaiCached() {
			return baiCacheDir.resolve(StringUtils.md5(this.url)+".bai");
			}
		
		private void downloadBaiToCache() {
			final Path baiCached = getBaiCached();
			if(Files.exists(baiCached)) {
				return;
				}
			final int max_try = 10;
			for(int i=0;i<max_try;i++) {
				LOG.info("download "+this.url+".bai try : "+(i+1)+"/"+max_try);
				try(final InputStream in= getBaiUrl().openStream()) {
					IOUtils.copyTo(in,baiCached);

					if(deleteBaiFilesOnExit) {
						Runtime.getRuntime().addShutdownHook(new Thread() {
							  public void run() {
							    try { Files.delete(baiCached);} catch(Exception err) {}
							  }
							});
						}
					return;
					}
				catch(final IOException err) {
					LOG.error(err);
					IOUtil.deletePaths(baiCached);
					}
				}
			throw new RuntimeIOException("Cannot download "+this.url+".bai");
			}
		SamReader open() throws IOException {
			final URL bamUrl = new URL(this.url);
			final SamInputResource sir = SamInputResource.of(bamUrl);
			sir.index(getBaiCached());
			return samReaderFactory.open(sir);
			}
		
		
		
		} 
	
	private final List<Sample1KG>  samples=  new ArrayList<>();
	private SAMSequenceDictionary dict=null;
	private final SamReaderFactory samReaderFactory = super.createSamReaderFactory();
	@Override
	public int doWork(final List<String> args) {
		SAMFileWriter sfw = null;
		QueryInterval queryIntervalArray[]=null;
		try
			{
			if(this.bedFile==null && !this.exitAfterBaiDownload) {
				LOG.error("bed file is missing");
				return -1;
				}
			
			IOUtil.assertDirectoryIsWritable(this.baiCacheDir);
			if(args.isEmpty()) {
				LOG.error("no urls");
				}
			else if(args.size()==1 && args.get(0).endsWith(".list")) {
				Files.lines(Paths.get(args.get(0))).
					filter(L->!StringUtils.isBlank(L)).
					filter(L->!L.startsWith("#")).
					map(L->new Sample1KG(L)).
					forEach(S->this.samples.add(S));
				}
			else
				{
				args.stream().
					map(L->new Sample1KG(L)).
					forEach(S->this.samples.add(S));
				}
			if(this.samples.isEmpty()) {
				LOG.error("no sample defined");
				return -1;
			}
			
			for(int x=0;x<this.samples.size();x++) {
				LOG.debug("bai "+(x+1)+"/"+this.samples.size());
				this.samples.get(x).downloadBaiToCache(); 
				}

			if(this.exitAfterBaiDownload) {
				LOG.info("Bai downloaded");				
				return 0;
				}
			for(final Sample1KG sample: this.samples) {
				LOG.debug("get sam file header for "+sample.url);
				try (final SamReader sr=sample.open()) {
					final SAMFileHeader header = sr.getFileHeader();
					if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
						throw new RuntimeIOException(sample.url+" is not sorted");
						}
					sample.sampleName = header.getReadGroups().stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().
							orElseThrow(()->new IllegalArgumentException("Cannot find sample in "+sample.url));
					final SAMSequenceDictionary dict = header.getSequenceDictionary();
					if(this.dict==null) {
						LOG.debug("read bed file "+this.bedFile);
						final List<QueryInterval> queryIntervals= new ArrayList<>();
						this.dict = dict;
						//load bed here
						ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(this.dict);
						final BedLineCodec codec = new BedLineCodec();
						BufferedReader br = IOUtils.openPathForBufferedReading(this.bedFile);
						String line;
						while((line=br.readLine())!=null) {
							final BedLine bed = codec.decode(line);
							if(bed==null) continue;
							String newCtg=ctgConverter.apply(bed.getContig());
							if(StringUtils.isBlank(newCtg)) {
								throw new JvarkitException.ContigNotFoundInDictionary(bed.getContig(), this.dict);
								}
							final int tid = this.dict.getSequenceIndex(newCtg);
							final QueryInterval queryInterval = new QueryInterval(tid, bed.getStart(), bed.getEnd());
							queryIntervals.add(queryInterval);
							}
						br.close();
						queryIntervalArray = QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));
						}
					else if(!SequenceUtil.areSequenceDictionariesEqual(dict, this.dict)) {
						throw new JvarkitException.DictionariesAreNotTheSame(dict, this.dict);
						}
					}
				}
			if(queryIntervalArray==null || queryIntervalArray.length==0) {
				LOG.error("no query interval defined");
				return -1;
			}
			
			final SAMRecordComparator samRecordComparator = new SAMRecordCoordinateComparator();
			
			final SAMFileHeader ouSamFileHeader= new SAMFileHeader(this.dict);
			ouSamFileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
			for(final Sample1KG sample: this.samples) {
				sample.readGroupRecord = new SAMReadGroupRecord(sample.sampleName);
				sample.readGroupRecord.setSample(sample.sampleName);
				sample.readGroupRecord.setLibrary(sample.sampleName);
				sample.readGroupRecord.setAttribute("URL", sample.url);
				ouSamFileHeader.addReadGroup(sample.readGroupRecord);
			}
			
			
			final SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
			sfwf.setCreateIndex(false);
			sfwf.setCreateMd5File(false);
			if(this.outputFile==null) {
				sfw = sfwf.makeBAMWriter(ouSamFileHeader, true,stdout());
				}
			else
				{
				sfw = sfwf.makeSAMOrBAMWriter(ouSamFileHeader, true, this.outputFile);
				}
		
			int interval_index = 0;
			while(interval_index <queryIntervalArray.length) {
				LOG.debug("interval_index= "+interval_index);
				/* current interval fetched by bam, can be updated */
				QueryInterval interval = queryIntervalArray[interval_index];
				LOG.debug("interval= "+interval);
				final List<Path> tmpFiles = new ArrayList<>(this.samples.size());
				boolean done=false;
				while(!done)  {
					done= true;
					for(final Path p:tmpFiles) Files.delete(p);
					tmpFiles.clear();
					
					for(final Sample1KG sample: this.samples) {
						if(!done) break;
						//LOG.debug("sample= "+sample.sampleName);
						//try forever until it works !
						for(;;)
							{ 
							//try to download data for the current interval
							final Path tmpBam = Files.createTempFile(this.baiCacheDir,"tmp.", ".bam");
							SamReader srf = null;
							SAMFileWriter rgnWriter = null;
							SAMRecordIterator samIter = null;
							try {
								srf = sample.open();
								//LOG.debug("query= "+sample.sampleName+" "+interval);
								samIter = srf.query(new QueryInterval[] {interval},false);
								rgnWriter =  sfwf.makeBAMWriter(ouSamFileHeader, true,tmpBam);
								while(samIter.hasNext()) {
									final SAMRecord rec=samIter.next();
									if(rec.isSecondaryOrSupplementary()) continue;
									if(rec.getReadUnmappedFlag()) continue;
									if(rec.getAlignmentEnd()>interval.end && /* current read goes beyond current interval */
										interval_index+1<queryIntervalArray.length && /* there is a new interval */
										interval.referenceIndex == queryIntervalArray[interval_index+1 ].referenceIndex && /* on same chromosome */
										queryIntervalArray[interval_index+1].start <= rec.getAlignmentEnd()
										) {
										// update current interval, merge with next
										interval = new QueryInterval(interval.referenceIndex,interval.start, queryIntervalArray[interval_index+1].end);
										interval_index++;
										LOG.info("extending interval to "+interval);
										done=false;
										break;
										}
									
									rec.setAttribute(SAMTag.RG.name(),sample.readGroupRecord.getId());
									rgnWriter.addAlignment(rec);
									}
								samIter.close();samIter=null;
								srf.close();srf=null;
								rgnWriter.close();rgnWriter=null;
								//LOG.debug("donequery= "+sample.url);
								}
							catch(final Exception err) {
								LOG.error(err);
								Files.delete(tmpBam);
								continue;
								}
							finally
								{
								CloserUtil.close(samIter);
								CloserUtil.close(srf);
								CloserUtil.close(rgnWriter);
								}
							
							tmpFiles.add(tmpBam);
							break;
							}
						if(!done) break;
						}//end of download each sample
					if(!done) continue;
					}//end of while(!done)
				//LOG.info("merging "+interval);
				//merge everything
				final List<SamReader> chunkReaders = new ArrayList<>(samples.size());
				final List<CloseableIterator<SAMRecord>> chunkIterators = new ArrayList<>(samples.size());
				for(final Path tmpPath:tmpFiles) {
					final SamReader tmpReader=samReaderFactory.open(tmpPath);
					chunkReaders.add(tmpReader);
					chunkIterators.add(tmpReader.iterator());
					}
				final MergingIterator<SAMRecord> merger = new MergingIterator<>(samRecordComparator, chunkIterators);
				while(merger.hasNext()) {
					sfw.addAlignment(merger.next());
					}
				
				merger.close();
				CloserUtil.close(chunkIterators);
				CloserUtil.close(chunkReaders);
				
				//cleanup
				for(final Path p:tmpFiles) Files.delete(p);
				interval_index++;
				}
			
			sfw.close();
			
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{

			}
		}
	
	public static void main(final String[] args) {
		new OneThousandBam().instanceMain(args);

	}

}
