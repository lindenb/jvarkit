package com.github.lindenb.jvarkit.bgen;

import java.io.Closeable;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.util.LongList;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;

public class BGenIndex extends BGenUtils {
	private static final byte[] MAGIC="BGIX1".getBytes(ENCODING);
	private static final int RECORD_SIZE=Integer.BYTES /* position */ + Long.BYTES /* offset */;
	protected final List<Contig> contigs = new ArrayList<>();
	protected final Map<String,Contig> name2contig = new HashMap<>();
	protected BGenIndex() {
		
		}
	
	public static IndexCreator newIndexCreator(final Path path) throws IOException {
		return new IndexCreatorImpl(path);
		}
	
	
	private static class Contig {
		int tid;
		String name;
		long countVariants=0L;
		long startOffset=0L;
		long endOffset() {
			return startOffset+countVariants*RECORD_SIZE;
			}
		}
	
	public interface IndexCreator extends AutoCloseable {
		public void add(final String contig,int position,long offset);
	}
	
	private static class Record implements Comparable<Record>{
		int tid;
		int position;
		long offset=0L;
		@Override
		public int compareTo(Record o) {
			int i = Integer.compare(tid, o.tid);
			if(i!=0) return i;
			i = Integer.compare(position, o.position);
			return i;
			}
		}
	private static class RecordCodec extends BinaryCodec implements SortingCollection.Codec<Record> {

		@Override
		public void encode(Record val) {
			writeInt(val.tid);
			writeInt(val.position);
			writeLong(val.offset);
			}

		@Override
		public Record decode() {
			Record rec=new Record();
			try {
				rec.tid = readInt();
				}
			catch(Throwable err)
				{
				return null;
				}
			rec.position= readInt();
			rec.offset= readLong();
			return rec;
			}
		@Override
		public RecordCodec clone()  {
			return new RecordCodec();
			}
		}
	private static class IndexCreatorImpl extends BGenIndex implements IndexCreator{
		final SortingCollection<Record> sorter;
		final Path indexPath;
		IndexCreatorImpl(Path indexPath) throws IOException {
			this.indexPath=indexPath;
			Path parentDir=indexPath.getParent();
			if(parentDir==null) throw new IOException("no parent for "+indexPath);
			this.sorter = SortingCollection.newInstance(
				Record.class,
				new RecordCodec(),
				(A,B)->A.compareTo(B),
				100_000,
				parentDir
				);
			}
		private int tid(String contig) {
			Contig ctg = this.name2contig.get(contig);
			if(ctg==null) {
				ctg = new Contig();
				ctg.name = contig;
				ctg.tid = this.name2contig.size();
				this.name2contig.put(contig, ctg);
				contigs.add(ctg);
				}
			ctg.countVariants++;
			return ctg.tid;
			}
		public void add(final String contig,int position,long offset)  {
			Record rec = new Record();
			rec.tid = tid(contig);
			rec.position = position;
			rec.offset = offset;
			sorter.add(rec);
			}
		@Override
		public void close() throws IOException {
			try {
				sorter.doneAdding();
						
				try(OutputStream os = Files.newOutputStream(this.indexPath)) {
					BinaryCodec codec = new BinaryCodec(os);
					codec.writeBytes(MAGIC);
					codec.writeInt(this.contigs.size());
					for(Contig c : this.contigs) {
						writeStringUInt16(codec, c.toString());
						codec.writeUInt(c.countVariants);
						}
					try(CloseableIterator<Record> it=sorter.iterator()) {
						final Record rec =it.next();
						codec.writeUInt(rec.position);
						codec.writeLong(rec.offset);
						}
					os.flush();
					}
				
				}
			finally {
				sorter.cleanup();
				}
			}
		}
	
	private static class IndexReader extends BGenIndex implements LongList<Record>,AutoCloseable{
		private FileInputStream fileInputStream=null; 
		private long num_variants=0L;
		private final long start_offset;
		private final BinaryCodec codec;
		IndexReader(Path path) throws IOException {
			fileInputStream = new FileInputStream(path.toFile());
			codec = new BinaryCodec(fileInputStream);
			long offset = 0L;
			byte[] magic = super.readNBytes(codec, MAGIC.length);
			offset+=magic.length;
			if(Arrays.compare(magic, MAGIC)!=0) throw new IOException("Bad magic");
			int n_contigs = codec.readInt();
			for(int i=0;i< n_contigs;++i) {
				Contig c = new Contig();
				int len = codec.readUShort();
				offset+=Short.BYTES;
				c.name = new String(readNBytes(codec, len),ENCODING);
				offset+=len;
				c.countVariants = codec.readUInt();
				offset+= Integer.BYTES;
				num_variants+=c.countVariants;
				c.tid =i;
				contigs.add(c);
				name2contig.put(c.name, c);
				}
			this.start_offset=offset;
			for(int i=0;i< this.contigs.size();i++) {
				Contig c = this.contigs.get(i);
				c.startOffset = offset;
				offset+= c.countVariants*(RECORD_SIZE);
				}
			}
		public long size() {
			return num_variants;
			}
		@Override
		public Record get(long idx) {
			final long offset = idx*RECORD_SIZE;
			Contig c = this.contigs.stream().filter(C->C.startOffset>=offset && offset<C.endOffset()).findFirst().get();
			try {
				this.fileInputStream.getChannel().position(offset);
				Record rec= new Record();
				rec.tid= c.tid;
				rec.position =  longToUnsignedInt(this.codec.readUInt());
				rec.offset = codec.readLong();
				return rec;
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		public void close() throws Exception {
			fileInputStream.close();
			codec.close();
			}
		}
	
	
	public static BGenIndex load(final Path path) throws IOException {
		BGenIndex instance= new BGenIndex();
		instance.read(path);
		return instance;
		}
	private void read(Path path) throws IOException  {
		try(InputStream is = Files.newInputStream(path)) {}
	}

}
