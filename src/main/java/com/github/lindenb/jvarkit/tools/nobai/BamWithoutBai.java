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
package com.github.lindenb.jvarkit.tools.nobai;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.PushbackInputStream;
import java.net.URL;
import java.nio.file.Path;
import java.util.List;

import org.apache.http.impl.client.HttpClientBuilder;
import org.disq_bio.disq.impl.formats.bgzf.BgzfBlockGuesser;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.CustomSeekableStreamFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.BamRecordGuesser;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import org.apache.http.impl.client.CloseableHttpClient;

/**
BEGIN_DOC

## Motivation

this tool queries a remote bam without bai (e.g: Encode bams - last checked 2019-12-16)

The idea is to run a **binary search** on the remote BAM, scanning for the BGZF blocks and the Bam Record.

This tool uses a class from the https://github.com/disq-bio/disq project , written by Tom White,  originally from the Hadoop-BAM project.

This tool expect 'small' reads. Long reads may fail.

## Acknowledgement

Louis Bergelson and John Marshall for their useful suggestions

## Example

```
$ java  -jar dist/bamwithoutbai.jar  -r "chr3:38548061-38649667" \
   -o jeter.bam https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam

$ samtools view jeter.bam
D2FC08P1:268:C3NPCACXX:8:1111:11204:70951	99	chr3	38548254	255	1S99M1S	=	38548442	289	CCACCTCCCACCCCCAACCACACCGTCTGCAGCCAGCCCCAGGCACCTGTCTCAAAGCTCCCGGGCTGTCCACACACACAAAAACCACAGTCTCCTTCCGC	@@@FFFFFFHHGHJIGIIJJIJJIGIIJGHGGHIIHHIIJJJGIIJJDGCHGHHHFEDFFFDCBDDBBBCDDDDDDDDDDB?BBDDBD@?>CCCCC@CC##	NH:i:1	HI:i:1	AS:i:198	NM:i:0	MD:Z:99
D2FC08P1:268:C3NPCACXX:8:1216:18650:7210	99	chr3	38548254	255	1S100M	=	38548442	289	CCACCTCCCACCCCCAACCACACCGTCTGCAGCCAGCCCCAGGCACCTGTCTCAAAGCTCCCGGGCTGTCCACACACACAAAAACCACAGTCTCCTTCCGG	CCCFFFFFHHHDHIJJJJJJIJJJJHHIJIJIJIJIIJIGIIJIIJJJJHHGHFHFFFFFFDDDBDDBBDCDDDDDDBBDDBDDDDDDDDCDDDDDDDD##	NH:i:1	HI:i:1	AS:i:199	NM:i:0	MD:Z:100
(...)
D2FC08P1:268:C3NPCACXX:8:1303:15665:9279	147	chr3	38649573	255	101M	=	38633218	-16456	TTGGCGCGGACTCGGCTCGGCGCGGGGCTCGGGGCACTGGGCGCAGGCTCAGCGGCCCCGGGGGAGCGATCCCTGCATCCTACGGGCGCCGCCGCCGTCTC	<9DDBB@C<2BB@DBDDBDBDDDDDBDB9BDDDDDCDDDDDDDDDDDDDDDDDDDDDDBDDDDDBBBDDDDDDDDDDDDDDBDDDDFJHHHHHFFFFFCCC	NH:i:1	HI:i:1	AS:i:196	NM:i:0	MD:Z:101
D2FC08P1:268:C3NPCACXX:8:2312:16447:12679	147	chr3	38649596	255	23S78M	=	38633222	-16452	GGGGCGGGGCCCGGGGGGGGGGGGGGGCGGGGGGCGCGGGGCGCAGCCTCGGCGCCCCGGGGGGAGCGATCCCTGCATCCTACGGGCGCCGCCGCCGTCTC	###################################################################################B>6<6F?DFFDDDDD@@@	NH:i:1	HI:i:1	AS:i:155	NM:i:8	MD:Z:5T0C5A1T8G3A3G3C42
```

comparing wget+samtools vs bamwithoutbai : 


```
$ cat interval.bed
chr22	41697506	41756151

$ time wget -q -O - "https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam" |\
	samtools view -L interval.bed -O BAM -owget.bam - 

real	17m48,710s
user	8m7,346s
sys	0m51,719s


$ samtools view wget.bam | sha1sum 
5fc261d051592edd791309211727ebbd0e3de909  -

time java -jar dist/bamwithoutbai.jar  -o nobai.bam -r "chr22:41697507-41756151" \
	 'https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam' 

real	0m24,848s
user	0m3,860s
sys	0m0,320s

$ samtools view nobai.bam | sha1sum 
5fc261d051592edd791309211727ebbd0e3de909  -
```


## Working behind a proxy.

see  https://docs.oracle.com/javase/8/docs/technotes/guides/net/proxies.html

```
$ java -Dhttp.proxyHost=webcache.example.com  -Dhttp.proxyPort=1234 \
    -Dhttps.proxyHost=webcache.example.com  -Dhttps.proxyPort=1234 \
    -jar dist/bamwithoutbai.jar  -r "chr3:38548061-38649667" \
   -o jeter.bam https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam
```

## Screenshot

https://twitter.com/yokofakun/status/1207337424935936001

![https://twitter.com/yokofakun/status/1207337424935936001](https://pbs.twimg.com/media/EMFS-xGXsAAVR20?format=jpg&name=small)


## See also

 * lbergelson on github: https://github.com/samtools/htsjdk/issues/1445#issuecomment-565599459
 * https://twitter.com/yokofakun/status/1202681859051859969
 * https://twitter.com/jomarnz/status/1205532441353560066

END_DOC

**/

@Program(name="bamwithoutbai",
description="Query a Remote BAM without bai",
keywords={"bam","sam","bai","remote"},
creationDate="20191213",
modificationDate="20191217"
)
public class BamWithoutBai extends Launcher{
	private static final Logger LOG = Logger.build(BamWithoutBai.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","--region","--interval"},description=IntervalParserFactory.OPT_DESC,required=true)
	private String intervalStr ="";
	@Parameter(names={"--repeat"},description="Max dichotomy repeat to perform during binary search.",hidden=true)
	private int dichotomy_repeat=50;
	@Parameter(names={"--reference","-R"},description="For writing CRAM. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;
	@Parameter(names={"--debug"},description="Enable debugging information.")
	private boolean do_debug = false;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	
	private CloseableHttpClient httpClient = null;
	
	/** empty closeable iterator */
	private static final CloseableIterator<SAMRecord> EOF_ITER  = new CloseableIterator<SAMRecord>() {
		@Override
		public SAMRecord next() { throw new IllegalStateException(); }
		@Override
		public boolean hasNext() { return false; }
		@Override
		public void close() { }
	};
	

	/** an iterator  converting an uncompressed input stream to SAMRecord */
	private class MySamRecordIterator 
		extends AbstractIterator<SAMRecord>
		implements CloseableIterator<SAMRecord> {
		/** bam record codec converting input stream to sam record */
        private final BAMRecordCodec bamRecordCodec;
        /** cleanup on close */
        private Runnable onClose = null;
        MySamRecordIterator(final SAMFileHeader header,final InputStream in,final String url) {
        	this.bamRecordCodec = new BAMRecordCodec(header);
        	this.bamRecordCodec.setInputStream(in,url);
        }
        
		@Override
		protected SAMRecord advance() {
            final SAMRecord next = this.bamRecordCodec.decode();
            if(next==null) close();
            return next;
			}
		@Override
		public void close() {
			if(onClose!=null) onClose.run();
			}
		}
	
	/** limit delegate iterator to the user's interval */
	private class MySamFilterIterator 
		extends AbstractIterator<SAMRecord>
		implements CloseableIterator<SAMRecord> {
		private final CloseableIterator<SAMRecord> delegate;
		/** user interval */
		private final QueryInterval interval;
		
		MySamFilterIterator(final CloseableIterator<SAMRecord> delegate,final QueryInterval interval) {
	    	this.delegate = delegate;
			this.interval = interval;
	    }
	    
		
		@Override
		protected SAMRecord advance() {
			while(this.delegate.hasNext()) {
				final SAMRecord rec = delegate.next();
				final int rec_tid= rec.getReferenceIndex();
				// unmapped are at the end.
				if(rec_tid==SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
					close();
					return null;
					}
				else if(rec_tid < this.interval.referenceIndex) {
					/* we're before the interval */
					continue;
					}
				else if(rec_tid > this.interval.referenceIndex) {
					/* we're beyond the interval */
					close();
					return null;
					}
				else if(rec.getStart() > this.interval.end) {
					/* we're beyond the interval */
					close();
					return null;
					}
				else if(rec.getEnd() < this.interval.start) {
					/* we're before the interval */
					continue;
					}
				return rec;
				}
			return null;
			}
		@Override
		public void close() {
			this.delegate.close();
			}
		}

	
	/* a custom SAMReader returning a closeable iterator for the user's interval */
	private class CustomSamReader implements Closeable {
		/* BAM url */
		final URL url;
		/* BAM header */
		final SAMFileHeader samFileHeader;
		/* seekable stream to the remote bam */
		final SeekableStream seekableStream;
		/* the underlying samReader */
		private SamReader samReader ;
		/* class finding the next BGZF block */
		final BgzfBlockGuesser bgzfBlockGuesser;
		
		CustomSamReader(final URL url ) throws IOException {
				this.url = url;
				final SamReaderFactory srf = SamReaderFactory.makeDefault().
						validationStringency(ValidationStringency.LENIENT);
				
				final CustomSeekableStreamFactory customSeekableStreamFactory = new CustomSeekableStreamFactory();
				// get Header
				if(do_debug) LOG.debug("opening "+url);

				
				this.seekableStream = customSeekableStreamFactory.
						setUserAgent(IOUtils.getDefaultUserAgent()).
						setUsingHttpHead(false).
						getStreamFor(url);
				
				this.bgzfBlockGuesser= new BgzfBlockGuesser(this.seekableStream, this.url.toString());

				/* get the header */ 
				this.samReader = srf.open(SamInputResource.of(this.seekableStream)); 
				
				this.samFileHeader = this.samReader.getFileHeader();
				
				// end try with resource
				if(!this.samFileHeader.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
					throw new IOException("Bam "+url+" is not sorted on coordinate: "+this.samFileHeader.getSortOrder());
					}
				}
		
		public SAMFileHeader getFileHeader() {
			return samFileHeader;
		}

		
		/* create a SAMRecord Iterator at the given offset in the rem */
		private CloseableIterator<SAMRecord> queryAtOffset(final long start_offset) throws IOException {
			if(do_debug) LOG.debug("opening at offset="+start_offset);
			
			if(start_offset>=this.seekableStream.length()) return EOF_ITER;
			
			final BgzfBlockGuesser.BgzfBlock bgzfBlock = this.bgzfBlockGuesser.guessNextBGZFPos(start_offset, seekableStream.length());
			if(bgzfBlock==null) {
				if(do_debug) LOG.debug("No bgz block found at "+start_offset);
				return EOF_ITER;
			}
			
			this.seekableStream.seek(bgzfBlock.pos);

			final BlockCompressedInputStream bcis = new BlockCompressedInputStream(this.seekableStream);
			final PushbackInputStream bpi = new PushbackInputStream(bcis,BamRecordGuesser.BUFFER_SIZE);
			final BamRecordGuesser bamRecordGuesser = new BamRecordGuesser(this.samFileHeader);
			bamRecordGuesser.setDebug(do_debug);
			if(bamRecordGuesser.find(bpi)) {
					if(do_debug) LOG.debug("Got SAMRecord at offset = "+start_offset);
					return new MySamRecordIterator(this.samFileHeader,bpi,this.url.toString());
				}
			return EOF_ITER;
			}
		
		
		/** return a SAM record iterator for the given genomic region */
		public CloseableIterator<SAMRecord> query(final Locatable locatable) throws IOException {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samFileHeader);			
			final int tid = dict.getSequenceIndex(locatable.getContig());
			if(tid==-1) {
				LOG.warn("unknown contig "+locatable.getContig());
				return EOF_ITER;
			}
			final QueryInterval userQueryInterval = new QueryInterval(tid, locatable.getStart(),locatable.getEnd());
			long byte_start = 0L;
			long byte_end=this.seekableStream.length();
			int repeat_dichotomy = dichotomy_repeat;
		    long len = byte_end - byte_start;
		    long last_offset_before= 0L;
		    SAMRecord previousBeforeRecord = null;
		    /* run binary search */
		    while (len > 0L)
		            {
		    		repeat_dichotomy--;
		            final long half = len / 2;
		            final long middle = byte_start + half;
		            
		            
		            CloseableIterator<SAMRecord> iter1 = queryAtOffset(middle);
		            /* something wrong happened or we've run too many binary search, start from here please */
		            if(!iter1.hasNext() || repeat_dichotomy==0) {
		            	if(do_debug) LOG.debug("break repeat="+repeat_dichotomy);
		            	break;
		            	}
		            
	            	/* we've run too many binary search, start from here please */

	            	final SAMRecord rec = iter1.next();
	            	if(do_debug) LOG.debug("Got SAM record at tid:"+rec.getReferenceIndex()+":"+rec.getStart());
	            	if(
	            		(rec.getReferenceIndex()!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && rec.getReferenceIndex()< userQueryInterval.referenceIndex) ||
	            		(rec.getReferenceIndex()== userQueryInterval.referenceIndex && rec.getAlignmentStart() <= userQueryInterval.start)
	            		)
	            		{
	            		final long save_byte_start = middle;
	            		if(do_debug) LOG.debug("tid:"+rec.getReferenceIndex()+":"+rec.getStart()+" before "+userQueryInterval.referenceIndex+":"+userQueryInterval.start);
	            		byte_start = middle + 1;
	                    len = len - half - 1;
	                    // are we **strictly** before the interval  rec.getAlignmentStart() < userQueryInterval.start ?
		            	if(
			            		(rec.getReferenceIndex()!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && rec.getReferenceIndex()< userQueryInterval.referenceIndex) ||
			            		(rec.getReferenceIndex()== userQueryInterval.referenceIndex && rec.getAlignmentStart() < userQueryInterval.start)
			            		)
		            		{
		            		if(do_debug) LOG.debug("tid:"+rec.getReferenceIndex()+":"+rec.getStart()+" strictly before "+userQueryInterval.referenceIndex+":"+userQueryInterval.start);

		            		if(do_debug) LOG.debug("new 'best' is offset="+save_byte_start);
		            		if(previousBeforeRecord!=null &&
		            			rec.getReadName().equals(previousBeforeRecord.getReadName()) &&
		            			rec.getFlags() == previousBeforeRecord.getFlags() &&
		            			rec.getStart() == previousBeforeRecord.getStart() &&
		            			rec.getEnd() == previousBeforeRecord.getEnd() &&
		            			rec.getReferenceIndex().equals(previousBeforeRecord.getReferenceIndex())) {
		            			if(do_debug) LOG.debug("new 'best' is same than previously record: break the loop");
		            			break;
		            			}
		            		last_offset_before = save_byte_start;
		            		previousBeforeRecord = rec;
		            		}
	            		}
	            	else
	            		{
	            		if(do_debug) LOG.debug("tid:"+rec.getReferenceIndex()+":"+rec.getStart()+" after "+userQueryInterval.referenceIndex+":"+userQueryInterval.start);
	            		len = half;
	            		}
		            }
    		if(do_debug) LOG.debug("best offset="+last_offset_before);
		    final CloseableIterator<SAMRecord> iter1 = queryAtOffset(last_offset_before);
		    if(!iter1.hasNext()) {
		    	if(do_debug) LOG.debug("iter has no next");
		    	return EOF_ITER;
				}
		    return new MySamFilterIterator(iter1,userQueryInterval);
			}
		
		
		@Override
		public void close()  {
			CloserUtil.close(this.seekableStream);
			this.bgzfBlockGuesser.close();
			}
	}
	
	
	
	@Override
	public int doWork(final List<String> args) {
	
	try {
		if(dichotomy_repeat<2) {
			LOG.error("dichotomy_repeat is too low "+dichotomy_repeat);
			return -1;
		}
		
		final String urlStr = oneAndOnlyOneFile(args);
		if(!( urlStr.startsWith("http://") || urlStr.startsWith("https://"))) {
			LOG.error("No a remote http url: "+urlStr);
			return -1;
		}

		if(do_debug ) {
			// https://hc.apache.org/httpcomponents-client-4.5.x/logging.html
			System.setProperty("org.apache.commons.logging.Log","org.apache.commons.logging.impl.SimpleLog");
			System.setProperty("org.apache.commons.logging.simplelog.showdatetime","true");
			System.setProperty("org.apache.commons.logging.simplelog.log.org.apache.http","DEBUG");
			System.setProperty("org.apache.commons.logging.simplelog.log.org.apache.http.wire","ERROR");
		}
		

		/** create http client */
		this.httpClient = HttpClientBuilder.
				create().
				setUserAgent(IOUtils.getDefaultUserAgent()).
				build();
		
		/* open SAM reader */
		try( CustomSamReader sr =  new CustomSamReader(new URL(urlStr))) {
			/* extract sam header */
			final SAMFileHeader samFileHeader = sr.getFileHeader();
			/* extract the SAM dictionary */
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samFileHeader);
			/* parse the user's interval*/
			final Locatable userInterval = IntervalParserFactory.newInstance(dict).make().apply(this.intervalStr).orElse(null);
			if(userInterval==null) {
				LOG.error("cannot parse interval "+this.intervalStr+" for "+urlStr);
				return -1;
				}
			/* create output SAM header */
			final SAMFileHeader header2 =  samFileHeader.clone();
			JVarkitVersion.getInstance().addMetaData(this, header2);
			/* open output SAM */
			try(SAMFileWriter sfw=this.writingBamArgs.setReferencePath(this.faidx).openSamWriter(this.outputFile, header2, true)) {
				/* query the interval and write in the output */
				try(CloseableIterator<SAMRecord> iter = sr.query(userInterval)) {
					while(iter.hasNext()) {
						final SAMRecord rec=iter.next();
						sfw.addAlignment(rec);
						}
					}
				}
			}
		/* we're done */
		return 0;
	} catch(final Throwable err) {
		LOG.error(err);
		return -1;
	} finally {
		CloserUtil.close(this.httpClient);
	}
	
	}
	
public static void main(final String[] args) {
	new BamWithoutBai().instanceMainWithExit(args);
	}
}
