/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriterFactory;
import com.github.lindenb.jvarkit.fastq.FastqRecordCodec;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;

/**
BEGIN_DOC


END_DOC
 
*/

@Program(name="repairfastq",
description="Join single end reads to paired end",
keywords="fastq",
modificationDate="20240128",
creationDate="20240128"
)
public class RepairFastq extends Launcher {
	private static final Logger LOG = Logger.of(RepairFastq.class);
	@Parameter(names={"-o","--output","--R1"},description="output for R1 reads. " + OPT_OUPUT_FILE_OR_STDOUT)
	private File fileoutR1 = null;
	@Parameter(names={"--R2"},description="fastq output for R2 reads.")
	private File fileoutR2 = null;
	@Parameter(names={"--Rx"},description="fastq output for other reads.")
	private File fileoutOthers = null;
	@Parameter(names={"--md5"},description="create MD5 checksum.")
	private boolean create_md5 = false;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private static class OneRead {
		short side;
		FastqRecord record;
		}
	
	private static class ReadCodec extends AbstractDataCodec<OneRead>
		{
		private final FastqRecordCodec fqCodec = new FastqRecordCodec();
		@Override
		public OneRead decode(final DataInputStream dis) throws IOException
			{
			final OneRead r=new OneRead();
			try {
				r.side = dis.readShort();
			} catch (final EOFException e) {
				return null;
				}
			r.record = this.fqCodec.decode(dis);
			return r;
			}
		@Override
		public void encode(final DataOutputStream dos,final OneRead r)
				throws IOException {
			dos.writeShort(r.side);
			this.fqCodec.encode(dos,r.record);
			}
		@Override
		public AbstractDataCodec<OneRead> clone() {
			return new ReadCodec();
			}
		}

	private OneRead mate(OneRead r) {
		String name= r.record.getReadName();
		OneRead or = new OneRead();
		switch(r.side) {
			case 1: or.side=2;break;
			case 2: or.side=1;break;
			default:or.side=0;break;
			}
		or.record =  new FastqRecord(name,"N",r.record.getBaseQualityHeader(),"2");
		return or;
		}
	
	
	private List<OneRead> fixPair( List<OneRead> array) {
		if(array.size()>2) {
			Map<String,OneRead> hash = array.stream().
					collect(Collectors.toMap(R->R.record.getReadString(),R->R));
			array=hash.values().stream().collect(Collectors.toList());
			}
		
		if(array.size()==2 && 
				array.get(0).record.getReadName().endsWith("/2") &&
				array.get(1).record.getReadName().endsWith("/1")
				) {
			array =  Arrays.asList(array.get(1),array.get(0));
			}
		
		if(array.size()==1) {
			return Arrays.asList(array.get(0),mate(array.get(0)));
			}
		return array;
		}
	
	
	
	private String normalizeName(final String s) {
		if(s.endsWith("/1") || s.endsWith("/2")) return s.substring(0,s.length()-2);
		return s;
		}
	
	@Override
	public int doWork(final List<String> args) {
		final Comparator<OneRead> readNameCompare1 = (R1,R2)-> normalizeName(R1.record.getReadName()).compareTo(normalizeName(R2.record.getReadName()));
		final Comparator<OneRead> readNameCompare2 = (R1,R2)-> {
			int i = readNameCompare1.compare(R1, R2);
			if(i!=0) return i;
			return Short.compare(R1.side, R2.side);
			};
		try {
			FastqPairedWriter pairedWriter;
			final FastqPairedWriterFactory fqwf = new FastqPairedWriterFactory().
					setCreateMd5(this.create_md5);
			
			if(fileoutR1==null && fileoutR2==null) { 
				pairedWriter = fqwf.open(stdout());
				}
			else if(fileoutR1!=null && fileoutR2==null) { 
				pairedWriter = fqwf.open(fileoutR1);
				}
			else  if(fileoutR1!=null && fileoutR2!=null) { 
				pairedWriter = fqwf.open(fileoutR1,fileoutR2);
				}
			else
				{
				LOG.error("fileout R1 undefined but R2 defined.");
				return -1;
				}
			final FastqWriter otherWiter;
			if(this.fileoutOthers!=null) {
				otherWiter = new BasicFastqWriter(this.fileoutOthers,this.create_md5);
			} else {
				otherWiter = null;
			}
			
			final SortingCollection<OneRead> sorter = SortingCollection.newInstance(
					OneRead.class,
					new ReadCodec(),
					readNameCompare2,
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths()
					);

			
			switch(args.size()) {
				case 0://continue
				case 1:
					{
					try(BufferedReader br = (args.isEmpty()?
							IOUtils.openStreamForBufferedReader(stdin()):
								IOUtils.openPathForBufferedReading(Paths.get(args.get(0))))) {
						try(FastqReader fr= new FastqReader(null, br, false)) {
							while(fr.hasNext()) {
								final OneRead r = new OneRead();
								r.side = 0;
								r.record = fr.next();
								sorter.add(r);
								}
							}
						}
					catch(IOException err) {
						LOG.error(err);
						return -1;
						}
					sorter.doneAdding();
					try(CloseableIterator<OneRead> iter1= sorter.iterator()) {
						EqualIterator<OneRead> iter = new EqualIterator<>(iter1,readNameCompare1);
						while(iter.hasNext()) {
							final List<OneRead> reads = fixPair(iter.next());
							if(reads.size()==2) {
								pairedWriter.write(reads.get(0).record, reads.get(1).record);
								}
							else if(otherWiter!=null)
								{
								for(OneRead rec: reads) {
									otherWiter.write(rec.record);
									}
								}
							}
						iter.close();
						}
					break;
					}
				case 2:
					{
					for(short side=0;side<2;++side) {
						final File file = new File(args.get(side));
						try(FastqReader fr= new FastqReader(file,false)) {
							while(fr.hasNext()) {
								final OneRead r = new OneRead();
								r.side = (short)(side+1);
								r.record = fr.next();
								sorter.add(r);
								}
							}
						}
					sorter.doneAdding();
					try(CloseableIterator<OneRead> iter1= sorter.iterator()) {
						final EqualIterator<OneRead> iter = new EqualIterator<>(iter1,readNameCompare1);
						while(iter.hasNext()) {
							final List<OneRead> reads = fixPair(iter.next());
							if(reads.size()==2 && reads.get(0).side==1 &&  reads.get(1).side==2) {
								pairedWriter.write(reads.get(0).record, reads.get(1).record);
								}
							else if(otherWiter!=null)
								{
								for(OneRead rec: reads) {
									otherWiter.write(rec.record);
									}
								}
							}
						iter.close();
						}
					break;
					}
				default:
					LOG.error("illegal number of arguments. Expected 0, 1 or 2 files");
					return -1;
				}
			sorter.cleanup();
			pairedWriter.close();
			if(otherWiter!=null) otherWiter.close();
			return 0;
			}
		catch(final  Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new RepairFastq().instanceMainWithExit(args);
	}

}
