package com.github.lindenb.jvarkit.util.picard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.regex.Pattern;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.Tribble;


public class TribbleFactory<T>
	extends AbstractIntervalMapFactory<T>
	{
	public class MyFeature<X>
		implements Feature
		{
		private Interval interval;
		private X value;
		public MyFeature(Interval interval,X value)
			{
			this.interval=interval;
			this.value=value;
			}
		@Override
		public String getChr() {
			return interval.getSequence();
			}
		@Override
		public String getContig() {
			return interval.getSequence();
		}
		@Override
		public int getStart() {
			return interval.getStart();
			}
		@Override
		public int getEnd() {
			return interval.getEnd();
			}
		public X getValue()
			{
			return value;
			}
		}
	
	@SuppressWarnings("rawtypes")
	private class MyCodec
		implements FeatureCodec<MyFeature,SRC>
			{
			private Pattern tab=Pattern.compile("[\t]");
			@Override
			public boolean canDecode(String f) {
				return f.endsWith(".tsv");
				}
			@Override
			public Class<MyFeature> getFeatureType()
				{
				return MyFeature.class;
				}
			private String readLine(PositionalBufferedStream stream) throws IOException
				{
				StringBuilder sb=null; 
				int c;
				while((c=stream.read())!=-1 && c!='\n')
					{
					if(sb==null) sb=new StringBuilder();
					sb.append((char)c);
					}	
				return sb==null?null:sb.toString();
				}
			@Override
			public MyFeature<T> decode(PositionalBufferedStream stream)
				throws IOException
				{
				for(;;)
					{
					String line=readLine(stream);
					if(line==null) return null;
					String tokens[]=tab.split(line);
					Interval interval=getKeyFunction().apply(tokens);
					if(interval==null) continue;
					T value=getValueFunction().apply(tokens);
					if(value==null) continue;
					return new MyFeature<T>(interval, value);
					}
				}
			@Override
			public Feature decodeLoc(PositionalBufferedStream stream)
				throws IOException {
				return decode(stream);
				}
			@Override
			public FeatureCodecHeader readHeader(PositionalBufferedStream stream)
				throws IOException {
				return FeatureCodecHeader.EMPTY_HEADER;
				}
			}
	
	private int maxRecordsInRAM=1000000;
	private static class StringCodec extends  AbstractDataCodec<String>
		{
		@Override
		public String decode(DataInputStream dis) throws IOException 
			{
			try
				{
				return dis.readUTF();
				}
			catch(IOException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, String object)
				throws IOException {
			dos.writeUTF(object);
			}
		@Override
		public StringCodec clone() {
			return new StringCodec();
			}
		}
	
	private class StringCmp implements Comparator<String>
		{
		@Override
		public int compare(String o1, String o2)
			{
			Pattern tab=Pattern.compile("[\t]");
			Interval i1=getKeyFunction().apply(tab.split(o1));
			Interval i2=getKeyFunction().apply(tab.split(o2));
			int i=compareChromosomes(i1.getSequence(),i2.getSequence());
			if(i!=0) return i;
			i=i1.getStart()-i2.getStart();
			if(i!=0) return i;
			i=i1.getEnd()-i2.getEnd();
			return i;
			}
		}
	
	public File createIndex(LineReader r,File file)  throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		SortingCollection<String> sorting=SortingCollection.newInstance(
				String.class,
				new StringCodec(),
				new StringCmp(),
				maxRecordsInRAM
				);
		sorting.setDestructiveIteration(true);
		String line;
		while((line=r.readLine())!=null)
	    	{
			 String tokens[]=tab.split(line);
	         Interval interval=getKeyFunction().apply(tokens);
	         T value=getValueFunction().apply(tokens);
	         if(value==null) continue;
	         if(interval==null) continue;
	         sorting.add(line);
	    	}
		sorting.doneAdding();
		
		//save to file
		
		CloseableIterator<String> iter=sorting.iterator();
		PrintWriter pw=new PrintWriter(file);
		while(iter.hasNext())
			{
			pw.println(iter.next());
			}
		pw.flush();
		pw.close();
		iter.close();
		sorting.cleanup();
		
		
		MyCodec codec=new MyCodec();
		File indexFile=Tribble.indexFile(file);
		Index index=IndexFactory.createDynamicIndex(file, codec);
		
		LittleEndianOutputStream out=new LittleEndianOutputStream(new FileOutputStream(indexFile));
		index.write(out);
		out.close();
		return indexFile;
		}
	
	public AbstractFeatureReader<MyFeature<T>> getTribble(File baseFile)
		{
		File indexFile=Tribble.indexFile(baseFile);
		Index index=IndexFactory.loadIndex(indexFile.toString());
		@SuppressWarnings("unchecked")
		AbstractFeatureReader<MyFeature<T>> reader= AbstractFeatureReader.getFeatureReader(
				baseFile.getAbsolutePath(),
				new MyCodec(),
				index
		         );
		return reader;
		}

	}
