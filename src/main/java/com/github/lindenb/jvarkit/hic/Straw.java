/*
  The MIT License (MIT)
 
  Copyright (c) 2011-2016 Broad Institute, Aiden Lab
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/

/* 
 copied from

https://github.com/aidenlab/straw/blob/24a50c4777e8992270402fa4465cd5f7d1dad8ca/C%2B%2B/straw.cpp

*/	
	
package com.github.lindenb.jvarkit.hic;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.zip.InflaterInputStream;

import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.util.LittleEndianInputStream;

// https://github.com/aidenlab/straw/blob/24a50c4777e8992270402fa4465cd5f7d1dad8ca/C%2B%2B/straw.cpp
	
public class Straw {
	
	static class XYValue
		{
		int x;
		int y;
		float value;
		}
	
	static class IndexEntry {
		  int size;
		  long position;
		}
	
	// sparse matrix entry
	static class ContactRecord {
	  int binX;
	  int binY;
	  float counts;
	  }
	
	/** an iterator dealing with IO exception */
	static private abstract class AbstractIOIterator<T>
		extends AbstractIterator<T>
		{
		@Override
		protected final T advance() {
			try {
				return take();
				}
			catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		/** the next element or null if the iterator is at the end */
		protected abstract T take() throws IOException;
		}
	
	// map of block numbers to pointers
	//final Map<Integer, IndexEntry> blockMap = new TreeMap<>();
	long total_bytes;
	// version number
	int version;
	/** genome version */
	private String buildVersion = "undefined";
	/** attributes */
	private final Map<String,String> attributes = new HashMap<>();
	/** dictionary */
	private SAMSequenceDictionary dictionary = null;
	
	private void debug(Object o) {
		System.err.println(o);
	}
	
	// reads the header, storing the positions of the normalization vectors and returning the master pointer
	long readHeader(LittleEndianInputStream fin) throws IOException {
	  final String magic = fin.readString();	
	  if (!magic.equals("HIC")) {
		  throw new IOException("Hi-C magic string is missing, does not appear to be a hic file "+ magic);
	  	}

	  this.version = fin.readInt();
	  if (this.version < 6) {
		  throw new IOException("Version " + version + " no longer supported");
	  	  }
	  final long master = fin.readLong();
	  this.buildVersion = fin.readString();
	  
	  // reading and ignoring attribute-value dictionary
	  final int nattributes = fin.readInt();
	  for (int i=0; i<nattributes; i++) {
		  final String key = fin.readString();//key
		  final String value = fin.readString();//value
		  this.attributes.put(key, value);
	  	  }
	  final int nChrs = fin.readInt();
	  final List<SAMSequenceRecord> ssrList = new ArrayList<SAMSequenceRecord>(nChrs);
	  // chromosome map for finding matrix
	  for (int i=0; i<nChrs; i++) {
		  final String name = fin.readString();
		  final int length = fin.readInt();
		  final SAMSequenceRecord ssr = new SAMSequenceRecord(name,length);
		  ssr.setAssembly(this.buildVersion);
		  ssrList.add(ssr);
	  	  }
	  this.dictionary = new SAMSequenceDictionary(ssrList);
	  return master;
	 }
	
	// reads the header, storing the positions of the normalization vectors and returning the master pointer
	private long readHeader(
			final LittleEndianInputStream fin,
			final String intervalStr1,
			final String intervalStr2,
			final QueryInterval queryIntervals[]
			) throws IOException {
		
	  final long master = readHeader(fin);	
	  Function<String, Optional<SimpleInterval>> intervalParser = IntervalParserFactory.newInstance().
			  dictionary(this.dictionary).
			  make();
	  
	  final SimpleInterval interval1 = intervalParser.apply(intervalStr1).orElseThrow(IntervalParserFactory.exception(intervalStr1));
	  final SimpleInterval interval2 = intervalParser.apply(intervalStr2).orElseThrow(IntervalParserFactory.exception(intervalStr2));

	  
	  queryIntervals[0] = new QueryInterval(this.dictionary.getSequenceIndex(interval1.getContig()),interval1.getStart(),interval1.getEnd());
	  queryIntervals[1] = new QueryInterval(this.dictionary.getSequenceIndex(interval2.getContig()),interval2.getStart(),interval2.getEnd());
	  return master;
	}
	
	void readFooter(LittleEndianInputStream fin, long master, int c1, int c2,
		  Normalization normx, 
		  Unit unitx, 
		  int resolution,
		  long myFilePos[], IndexEntry c1NormEntry, IndexEntry c2NormEntry) throws IOException {
		  int nBytes = fin.readInt();
		  String key =String.valueOf(c1)+"_"+c2;
		  
		  int nEntries = fin.readInt();
		  boolean found = false;
		  for (int i=0; i<nEntries; i++) {
		    String str = fin.readString();
		    //System.err.println("str: "+str);
		    long fpos = fin.readLong();
		    //System.err.println("fpos: "+fpos);
		    int sizeinbytes = fin.readInt();
		    if (str.equals(key)) {
		      myFilePos[0] = fpos;
		      System.err.println("Got myFilePos="+fpos);
		      found=true;
		    }
		  }
		  if (!found) {
		    throw new IOException("File doesn't have the given chr_chr map" );
		  	}
		  if (normx.equals(Normalization.NONE)) return; // no need to read norm vector index
		  // read in and ignore expected value maps; don't store; reading these to 
		  // get to norm vector index
		  int nExpectedValues = fin.readInt();
		  //System.err.println("nExpectedValues="+nExpectedValues);
		  
		  for (int i=0; i<nExpectedValues; i++) {
		    String str = fin.readString();

		    int binSize = fin.readInt();
		    //System.err.println("binSize="+binSize);


		    int nValues = fin.readInt();
		    //System.err.println("nValues="+nValues);

		    for (int j=0; j<nValues; j++) {
		      double v = fin.readDouble();
		     // System.err.println("v="+v);
		    }

		    int nNormalizationFactors = fin.readInt();
		    //System.err.println("nNormalizationFactors["+i+"]="+nNormalizationFactors);
		    for (int j=0; j<nNormalizationFactors; j++) {
		      int chrIdx= fin.readInt();
		      //System.err.println("chrIdx["+j+"]="+chrIdx);

		      double v= fin.readDouble();
		      //System.err.println("v["+j+"]="+v);
		    }
		  }
		 
		  
		  nExpectedValues = fin.readInt();
		  //System.err.println("nExpectedValues "+nExpectedValues);
		  
		  for (int i=0; i<nExpectedValues; i++) {
		    String ignore1 = fin.readString(); //typeString
		    String ignore2 = fin.readString(); //unit
		    int binSize = fin.readInt();
		    //System.err.println("nExpectedValues["+i+"] "+ignore1 +" "+ignore2+" "+binSize);
		    

		    int nValues = fin.readInt();
		    for (int j=0; j<nValues; j++) {
		      double v= fin.readDouble();
		    }
		    int nNormalizationFactors = fin.readInt();
		    for (int j=0; j<nNormalizationFactors; j++) {
		      int chrIdx = fin.readInt();
		      double v= fin.readDouble();
		    }
		  }
		  // Index of normalization vectors
		  nEntries = fin.readInt();
		  boolean found1 = false;
		  boolean found2 = false;
		  
		  final Set<String> available_rom_resolutions = new HashSet<>();
		  
		  for (int i = 0; i < nEntries; i++) {
		    final Normalization normtype = Normalization.valueOf(fin.readString());
		    int chrIdx = fin.readInt();
		    final Unit unit1 = Unit.valueOf(fin.readString());
		    final int resolution1 = fin.readInt();
		    
		    available_rom_resolutions.add(normtype+"("+resolution1+")");
		    
		    long filePosition = fin.readLong();
		    int sizeInBytes= fin.readInt();
		    if (chrIdx == c1 && normtype.equals(normx) && unit1.equals(unitx) && resolution1 == resolution) {
		      c1NormEntry.position=filePosition;
		      c1NormEntry.size=sizeInBytes;
		      found1 = true;
		    }
		    if (chrIdx == c2 && normtype.equals(normx) && unit1.equals(unitx) && resolution1 == resolution) {
		      c2NormEntry.position=filePosition;
		      c2NormEntry.size=sizeInBytes;
		      found2 = true;
		    }
		  }
		  if (!found1 || !found2) {
		    throw new IOException(
		    		"Normalization vectors not found for one or both chromosomes at " + resolution + " " + unitx+" available "+ String.join("\n", available_rom_resolutions));
		  }
		}

	
	private class ReadMatrixZoomDataResult
		{
		int blockBinCount;
		int blockColumnCount;
		final Map<Integer, IndexEntry> blockMap = new TreeMap<>();
		}
	
	// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count 
	private ReadMatrixZoomDataResult readMatrixZoomData(LittleEndianInputStream fin, Unit myunitx, int mybinsize) throws IOException 
	  {
	  final Unit unit = Unit.valueOf(fin.readString());
	  fin.readInt(); // Old "zoom" index -- not used
	  fin.readFloat(); // sumCounts
	  fin.readFloat(); // occupiedCellCount
	  fin.readFloat(); // stdDev
	  fin.readFloat(); // percent95
	  final int binSize = fin.readInt();
	  final int blockBinCount = fin.readInt();
	  final int blockColumnCount = fin.readInt();
	  debug("readMatrixZoomData "+unit +" "+binSize+" "+blockBinCount+" "+blockColumnCount);
	  final ReadMatrixZoomDataResult storeBlockData;
	  
	  if (myunitx.equals(unit) && mybinsize==binSize) {
		storeBlockData = new ReadMatrixZoomDataResult();
		storeBlockData.blockBinCount = blockBinCount;
		storeBlockData.blockColumnCount = blockColumnCount;
	    }
	  else
	  	{
		storeBlockData = null;
	  	}
	  
	  final int nBlocks  = fin.readInt();

	  for (int b = 0; b < nBlocks; b++) {
	    final int blockNumber = fin.readInt();
	    final long filePosition = fin.readLong();
	    final int blockSizeInBytes = fin.readInt();
	    final IndexEntry entry = new IndexEntry();
	    entry.size = blockSizeInBytes;
	    entry.position = filePosition;
	    if (storeBlockData!=null) {
	    	debug(" storeBlockData "+blockNumber);
	    	storeBlockData.blockMap.put(blockNumber, entry);
	    	}
	  	}
	  return storeBlockData;
	}	
	
	
	private static Set<Integer> getBlockNumbersForRegionFromBinPosition(
			int regionIndices[], 
			int blockBinCount, 
			int blockColumnCount,
			boolean intra
			) {
		   final int col1 = regionIndices[0] / blockBinCount;
		   final int col2 = (regionIndices[1] + 1) / blockBinCount;
		   final int row1 = regionIndices[2] / blockBinCount;
		   final int row2 = (regionIndices[3] + 1) / blockBinCount;
		   
		   final Set<Integer> blocksSet = new TreeSet<>();
		   // first check the upper triangular matrix
		   for (int r = row1; r <= row2; r++) {
		     for (int c = col1; c <= col2; c++) {
		       final int blockNumber = r * blockColumnCount + c;
		       blocksSet.add(blockNumber);
		     }
		   }
		   // check region part that overlaps with lower left triangle
		   // but only if intrachromosomal
		   if (intra) {
		     for (int r = col1; r <= col2; r++) {
		       for (int c = row1; c <= row2; c++) {
			     final int blockNumber = r * blockColumnCount + c;
			     blocksSet.add(blockNumber);
		       }
		     }
		   }

		return blocksSet;
		}

	private ReadMatrixZoomDataResult readMatrix(
		final SeekableStream fin,
		long myFilePosition,
		final Unit unitx,
		int resolution
		) throws IOException {
		debug("readMatrix at "+ myFilePosition);  

		fin.seek(myFilePosition);  
		  LittleEndianInputStream in = new LittleEndianInputStream(fin);
		  in.skip(Integer.BYTES * 2);//ignore c1 + c2
		  final int nRes = in.readInt();
		  int i=0;
		  while (i<nRes) {
			debug("" + i + "/" + nRes);
		    final ReadMatrixZoomDataResult found = this.readMatrixZoomData(in, unitx, resolution);
		    if(found!=null) return found;
		    i++;
		    }
		throw new IOException("Error finding block data");
		}
	
	private class ReadBlockIteratorV6 extends AbstractIOIterator<ContactRecord>{
		private final LittleEndianInputStream in;
		private final int nRecords;
		private int index=0;
		ReadBlockIteratorV6(LittleEndianInputStream in,int nRecords) {
			this.in = in;
			this.nRecords = nRecords;
			}
		@Override
		protected ContactRecord take() throws IOException {
			if(this.index >= this.nRecords) return null;
			this.index++;
			final ContactRecord record=new ContactRecord();
			record.binX = in.readInt();
	        record.binY = in.readInt();
	        record.counts = in.readFloat();
			return record;
			}
		}
	
	private class ReadBlockIteratorV7 extends AbstractIOIterator<ContactRecord>{
		private final LittleEndianInputStream in;
		private final int binXOffset;
		private final int binYOffset;
		private final char useShort;
		private final char type;
		private final short rowCount;
		private short colCount = 0;
		private int i = 0;
		private int j = 0;
		private short y = 0;
		ReadBlockIteratorV7(final LittleEndianInputStream in) throws IOException{
			this.in = in;
			this.binXOffset = in.readInt();
			this.binYOffset = in.readInt();
			this.useShort = (char)in.readByte();
			this.type = (char)in.readByte();
			if(this.type==1) {
				this.rowCount = in.readShort();
				}
			else
				{
				this.rowCount = 0;
				}
			}
		@Override
		protected ContactRecord take() throws IOException {
			if(this.type!=1) return null;
		    
		    for(;;) {
		    	 if( this.i >= this.rowCount) return null;
		    	
		    	 if(this.j==0)
		    	 	{
		    		this.y = in.readShort() ;
		    		this.colCount = in.readShort();
		    	 	}

		    	 if(this.j >= this.colCount) {
		    		this.j = 0;
		    		this.i++;
		    		if( this.i >= this.rowCount) return null;
		    		continue;
		    	 	}
		    	 
		    	
		    	 final short x = in.readShort();
		    	 final ContactRecord record = new ContactRecord();

				  if (this.useShort == 0) { // yes this is opposite of usual
					record.counts = (float)this.in.readShort();
				    } 
				  else {
					record.counts = this.in.readFloat();
				  	}
		    	
				 record.binX = this.binXOffset + x;;
				 record.binY = this.y + this.binYOffset;
				 
				 this.j++;
				 
				 return record;
		    	}
		    }	
		}

	
	private Iterator<ContactRecord> readBlock(
			final SeekableStream fin,
			final int blockNumber,
			final Map<Integer, IndexEntry> blockMap) throws IOException {
		debug("readblock at "+fin.position()+" blck "+blockNumber);  
		IndexEntry idx = blockMap.get(blockNumber);
		if (idx==null) {
			debug("No block for " +blockNumber );
		    return Collections.emptyIterator();
		  }
		
			if (idx.size == 0) {
			    return Collections.emptyIterator();
			  }
		  debug("reading "+idx.size);
		  final byte compressedBytes[] = new byte[idx.size];
     
		  fin.seek(idx.position);
		  fin.read(compressedBytes);
		  final InflaterInputStream zipIn = new InflaterInputStream(new ByteArrayInputStream(compressedBytes));
		  
		  
		  // create stream from buffer for ease of use

		  @SuppressWarnings("resource")
		  final LittleEndianInputStream bufferin=new LittleEndianInputStream(zipIn);
		  final int nRecords = bufferin.readInt();
		  debug("NRecord "+nRecords);
		  
		  // different versions have different specific formats
		  if (this.version < 7) {
			return new ReadBlockIteratorV6(bufferin,nRecords);
		    }
		  else {
			return new ReadBlockIteratorV7(bufferin);
		    }
		  }
	
	
	// reads the normalization vector from the file at the specified location
	private static double[] readNormalizationVector(final ByteArrayInputStream is) throws IOException {
	  @SuppressWarnings("resource")
	  final LittleEndianInputStream in = new LittleEndianInputStream(is);
	  final int nValues = in.readInt();
	  double values[] = new double[nValues];
	  for (int i = 0; i < nValues; i++) {
	     values[i]=in.readDouble();
	  	 }
	  return values;
	  }
	
	static class Query {
		Locatable interval1;
		Locatable interval2;
		}
	
	Iterator<XYValue> straw(final Normalization norm, String fname, int binsize, String chr1loc, String chr2loc, Unit unit) throws IOException
	 {
	 final QueryInterval queryIntervals[]=new QueryInterval[] {null,null};
		

	  // read header into buffer; 100K should be sufficient
	 
	  SeekableStream seekIn = SeekableStreamFactory.getInstance().getStreamFor(fname);
	  LittleEndianInputStream fin = new LittleEndianInputStream(new BufferedInputStream(seekIn));
	  long master = this.readHeader(fin,chr1loc,chr2loc,queryIntervals);
	  debug("master="+master);

	  // from header have size of chromosomes, set region to read
	  int c1=Math.min(queryIntervals[0].referenceIndex,queryIntervals[1].referenceIndex);
	  int c2=Math.max(queryIntervals[0].referenceIndex,queryIntervals[1].referenceIndex);
	  int origRegionIndices[]=new int[4]; // as given by user
	  int regionIndices[]=new int[4]; // used to find the blocks we need to access
	  // reverse order if necessary
	  if (queryIntervals[0].referenceIndex > queryIntervals[1].referenceIndex) {
	    origRegionIndices[0] = queryIntervals[1].start;
	    origRegionIndices[1] = queryIntervals[1].end;
	    origRegionIndices[2] = queryIntervals[0].start;
	    origRegionIndices[3] = queryIntervals[0].end;
	  	}
	  else {
	    origRegionIndices[0] = queryIntervals[0].start;
	    origRegionIndices[1] = queryIntervals[0].end;
	    origRegionIndices[2] = queryIntervals[1].start;
	    origRegionIndices[3] = queryIntervals[1].end;
	    }
	  regionIndices[0] = origRegionIndices[0] / binsize;
	  regionIndices[1] = origRegionIndices[1] / binsize;
	  regionIndices[2] = origRegionIndices[2] / binsize;
	  regionIndices[3] = origRegionIndices[3] / binsize;
	  

	  IndexEntry c1NormEntry=new IndexEntry(), c2NormEntry=new IndexEntry();
	  long myFilePos[]=new long[] {0};

	  //long bytes_to_read = total_bytes - master;
	 debug("master:"+master);
	
	  seekIn.seek(master);
	  fin = new LittleEndianInputStream(seekIn);
	  readFooter(fin, master, c1, c2, norm, unit, binsize,
			  myFilePos,
			  c1NormEntry, c2NormEntry); 
	  
	  // readFooter will assign the above variables


	  double c1Norm[]=new double[0];
	  double c2Norm[]=new double[0];

	  if (!norm.equals(Normalization.NONE)) {
		
	    byte buffer3[] = new byte[c1NormEntry.size];
	    debug("seek1 "+c1NormEntry.position);
	    seekIn.seek(c1NormEntry.position);
	    seekIn.read(buffer3);
	    
	    ByteArrayInputStream bufferin = new ByteArrayInputStream(buffer3);
	    c1Norm = readNormalizationVector(bufferin);

	    byte buffer4[] = new byte[c2NormEntry.size];
	    debug("seek2 "+c2NormEntry.position);
	    seekIn.seek(c2NormEntry.position);
	    seekIn.read(buffer4);
	    
	  
	    
	    bufferin = new ByteArrayInputStream(buffer4);
	    c2Norm = readNormalizationVector(bufferin);
	  }

	  
	  // readMatrix will assign blockBinCount and blockColumnCount
	  final ReadMatrixZoomDataResult readMatrixResult = readMatrix(seekIn, myFilePos[0], unit, binsize); 
	  
	  final Set<Integer> blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, readMatrixResult.blockBinCount, readMatrixResult.blockColumnCount, c1==c2); 

	  // getBlockIndices
	  final List<XYValue> xyvalues = new ArrayList<>();
	  for (Integer it:blockNumbers) {
	    // get contacts in this block
		final  Iterator<ContactRecord> recordsIterator = readBlock(seekIn,it,readMatrixResult.blockMap);
	    while(recordsIterator.hasNext()) {     
	      final ContactRecord rec = recordsIterator.next();
	      int x = rec.binX * binsize;
	      int y = rec.binY * binsize;
	      float c = rec.counts;
	      if (!norm .equals(Normalization.NONE)) {
			  c = (float)(c / (c1Norm[rec.binX] * c2Norm[rec.binY]));
		      }

	      if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
		   y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
		  // or check regions that overlap with lower left
		  ((c1==c2) && y >= origRegionIndices[0] && y <= origRegionIndices[1] && x >= origRegionIndices[2] && x <= origRegionIndices[3])) {
		
	    	  XYValue xyv = new XYValue();
	    	  xyv.x=x;
	    	  xyv.y=y;
	    	  xyv.value=c;
	    	  xyvalues.add(xyv);
	      }
	    }
	  }
    return xyvalues.iterator();
	}
	
	
	
public static void main(String[] argv) {
	 if (argv.length != 6) {
		   System.err.println( "Not enough arguments");
		   System.err.println( "Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>");
		  System.exit(-1);
		  }
	 try {
		  Normalization norm=Normalization.valueOf(argv[0]);
		  String fname=argv[1];
		  String chr1loc=argv[2];
		  String chr2loc=argv[3];
		  Unit unit=Unit.valueOf(argv[4]);
		  String size=argv[5];
		  int binsize=Integer.parseInt(size);
		  Straw app=new Straw();
		  Iterator<XYValue> iter = app.straw(norm, fname, binsize, chr1loc, chr2loc, unit);
		  while(iter.hasNext()) {
		    XYValue xyv = iter.next();
			System.out.printf("%d\t%d\t%.14g\n", xyv.x, xyv.y, xyv.value);   
		  	}
	}	catch(Throwable err) {
		err.printStackTrace();
	}
}

}
