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
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.InflaterInputStream;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.util.LittleEndianInputStream;

// https://github.com/aidenlab/straw/blob/24a50c4777e8992270402fa4465cd5f7d1dad8ca/C%2B%2B/straw.cpp
	
public class Straw {
	enum  Normalization {NONE,VC,VC_SQRT,KR};
	enum Unit {BP,FRAG};
	
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
	
	// map of block numbers to pointers
	final Map<Integer, IndexEntry> blockMap = new TreeMap<>();
	long total_bytes;
	// version number
	int version;
	
	boolean readMagicString(LittleEndianInputStream fin) throws IOException {
		  String str = fin.readString();
		  return str.equals("HIC");
		}
	
	// reads the header, storing the positions of the normalization vectors and returning the master pointer
	long readHeader(LittleEndianInputStream fin, String chr1, String chr2, 
			int c1pos1[], int c1pos2[],
			int c2pos1[], int c2pos2[],
			int chr1ind[], int chr2ind[]) throws IOException {
	  if (!readMagicString(fin)) {
		  throw new IOException("Hi-C magic string is missing, does not appear to be a hic file");
	  	}

	  this.version = fin.readInt();
	  if (this.version < 6) {
		  throw new IOException("Version " + version + " no longer supported");
	  }
	  long master = fin.readLong();
	  final String genome = fin.readString();
	  int nattributes = fin.readInt();

	  // reading and ignoring attribute-value dictionary
	  for (int i=0; i<nattributes; i++) {
		  final String key = fin.readString();//key
		  final String value = fin.readString();//value
	  	  }
	  int nChrs = fin.readInt();
	  // chromosome map for finding matrix
	  boolean found1 = false;
	  boolean found2 = false;
	  for (int i=0; i<nChrs; i++) {
	    String name = fin.readString();
	    int length = fin.readInt();
	    if (name.equals(chr1)) {
	      found1=true;
	      chr1ind[0] = i;
	      if (c1pos1[0] == -100) {
			c1pos1[0] = 0;
			c1pos2[0] = length;
	      }
	    }
	    if (name.equals(chr2)) {
	      found2=true;
	      chr2ind[0] = i;
	      if (c2pos1[0] == -100) {
		c2pos1[0] = 0;
		c2pos2[0] = length;
	      }
	    }
	  }
	  if (!found1 || !found2) {
	   throw new IOException("One of the chromosomes wasn't found in the file. Check that the chromosome name matches the genome.");
	  }
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
			  //System.err.println("[i]="+i);
		    String str = fin.readString();
		    //System.err.println("str="+str);

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

	
	// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count 
	boolean readMatrixZoomData(LittleEndianInputStream fin, Unit myunitx, int mybinsize, int myBlockBinCount[], int myBlockColumnCount[]) throws IOException {
	  final Unit unit = Unit.valueOf(fin.readString());
	  fin.readInt(); // Old "zoom" index -- not used
	  fin.readFloat(); // sumCounts
	  fin.readFloat(); // occupiedCellCount
	  fin.readFloat(); // stdDev
	  fin.readFloat(); // percent95
	  final int binSize = fin.readInt();
	  final int blockBinCount = fin.readInt();
	  final int blockColumnCount = fin.readInt();
	  System.err.println("readMatrixZoomData "+unit +" "+binSize+" "+blockBinCount+" "+blockColumnCount);
	  boolean storeBlockData = false;
	  if (myunitx.equals(unit) && mybinsize==binSize) {
	    myBlockBinCount[0] = blockBinCount;
	    myBlockColumnCount[0] = blockColumnCount;
	    storeBlockData = true;
	  }
	  
	  final int nBlocks  = fin.readInt();

	  for (int b = 0; b < nBlocks; b++) {
	    final int blockNumber = fin.readInt();
	    final long filePosition = fin.readLong();
	    final int blockSizeInBytes = fin.readInt();
	    final IndexEntry entry = new IndexEntry();
	    entry.size = blockSizeInBytes;
	    entry.position = filePosition;
	    if (storeBlockData) {
	    	System.err.println(" storeBlockData "+blockNumber);
	    	this.blockMap.put(blockNumber, entry);
	    	}
	  	}
	  return storeBlockData;
	}	
	
	
	private Set<Integer> getBlockNumbersForRegionFromBinPosition(
			int regionIndices[], 
			int blockBinCount, 
			int blockColumnCount,
			boolean intra) {
		   int col1 = regionIndices[0] / blockBinCount;
		   int col2 = (regionIndices[1] + 1) / blockBinCount;
		   int row1 = regionIndices[2] / blockBinCount;
		   int row2 = (regionIndices[3] + 1) / blockBinCount;
		   
		   Set<Integer> blocksSet = new TreeSet<>();
		   // first check the upper triangular matrix
		   for (int r = row1; r <= row2; r++) {
		     for (int c = col1; c <= col2; c++) {
		       int blockNumber = r * blockColumnCount + c;
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

	private void readMatrix(final SeekableStream fin, long myFilePosition, Unit unitx, int resolution, int myBlockBinCount[], int myBlockColumnCount[]) throws IOException {
		System.err.println("readMatrix at "+ myFilePosition);  

		fin.seek(myFilePosition);  
		  LittleEndianInputStream in = new LittleEndianInputStream(fin);
		  in.skip(Integer.BYTES * 2);//ignore c1 + c2
		  int nRes = in.readInt();
		  int i=0;
		  boolean found=false;
		  while (i<nRes && !found) {
			 System.err.println("" + i + "/" + nRes);
		    found = readMatrixZoomData(in, unitx, resolution, myBlockBinCount, myBlockColumnCount);
		    i++;
		    }
		if (!found) throw new IOException("Error finding block data");
		}
	
	
	private List<ContactRecord> readBlock(SeekableStream fin, Object curl, boolean isHttp, int blockNumber) throws IOException {
		System.err.println("readblock at "+fin.position()+" blck "+blockNumber);  
		IndexEntry idx = blockMap.get(blockNumber);
		if (idx==null) {
			 System.err.println("No block for " +blockNumber );
		    return Collections.emptyList();
		  }
		
		if (idx.size == 0) {
		    return Collections.emptyList();
		  }
		  System.err.println("reading "+idx.size);
		  final byte compressedBytes[] = new byte[idx.size];
     
		
		  fin.seek(idx.position);
		  fin.read(compressedBytes);
		  
		 final InflaterInputStream zipIn = new InflaterInputStream(new ByteArrayInputStream(compressedBytes));
		 final ByteArrayOutputStream baos = new ByteArrayOutputStream();
		 byte[] buffer = new byte[1024];   
		 int count;
		 while ((count=zipIn.read(buffer))!=-1) {  
		    baos.write(buffer, 0, count);   
		    }  
		 
		 baos.close();  
		  
		  final byte uncompressedBytes[] = baos.toByteArray();
		  baos.close();
		  
		  System.err.println("Uncompressed size:"+uncompressedBytes.length);
		  

		  // create stream from buffer for ease of use

		  @SuppressWarnings("resource")
		LittleEndianInputStream bufferin=new LittleEndianInputStream(new ByteArrayInputStream(uncompressedBytes));
		  int nRecords = bufferin.readInt();
		  System.err.println("NRecord "+nRecords);
		  List<ContactRecord> v=new ArrayList<>(nRecords);
		  // different versions have different specific formats
		  if (version < 7) {
		    for (int i = 0; i < nRecords; i++) {
		      int binX= bufferin.readInt();
		      int binY= bufferin.readInt();
		      float counts= bufferin.readFloat();
		      ContactRecord record=new ContactRecord();
		      record.binX = binX;
		      record.binY = binY;
		      record.counts = counts;
		      v.add(record);
		    }
		  } 
		  else {
		    int binXOffset = bufferin.readInt();
		    int binYOffset = bufferin.readInt();
		    char useShort = (char)bufferin.readByte();
		    char type = (char)bufferin.readByte();
		    int index=0;
		    if (type == 1) {
		      // List-of-rows representation
		      short rowCount = bufferin.readShort();
		      for (int i = 0; i < rowCount; i++) {
			short y =bufferin.readShort() ;
			int binY = y + binYOffset;
			short colCount = bufferin.readShort();
			for (int j = 0; j < colCount; j++) {
			  short x = bufferin.readShort();
			  int binX = binXOffset + x;
			  float counts;
			  if (useShort == 0) { // yes this is opposite of usual
			    short c2 = bufferin.readShort();
			    counts = c2;
			  } 
			  else {
				  counts = bufferin.readFloat();
			  }
			  ContactRecord record = new ContactRecord();
			  record.binX = binX;
			  record.binY = binY;
			  record.counts = counts;
			  v.add(record);
			  index++;
			}
		      }
		    }
		  }
		  return v;
		}
	
	// reads the normalization vector from the file at the specified location
	double[] readNormalizationVector(final ByteArrayInputStream is) throws IOException {
	  LittleEndianInputStream in = new LittleEndianInputStream(is);
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
	
	List<XYValue> straw(final Normalization norm, String fname, int binsize, String chr1loc, String chr2loc, Unit unit) throws IOException
	 {
	  // parse chromosome positions
	  String ss[]= chr1loc.split("[\\:\\-]");
	  String chr1, chr2;
	  int c1pos1[]=new int[] {-100};
	  int c1pos2[]=new int[] {-100};
	  int c2pos1[]=new int[] {-100};
	  int c2pos2[]=new int[] {-100};
	  chr1=ss[0];
	    c1pos1[0] = Integer.parseInt(ss[1]);
	    c1pos2[0] = Integer.parseInt(ss[2]);
	  
	  ss= chr2loc.split("[\\:\\-]");
	  chr2=ss[0];
	  c2pos1[0] = Integer.parseInt(ss[1]);
	  c2pos2[0] = Integer.parseInt(ss[2]);
	    
	 
	  int[] chr1ind=new int[] {0};
	  int[] chr2ind=new int[] {0};

	  boolean isHttp = false;

	  // read header into buffer; 100K should be sufficient
	  long master;
	  SeekableStream seekIn = SeekableStreamFactory.getInstance().getStreamFor(fname);
	  LittleEndianInputStream fin = new LittleEndianInputStream(new BufferedInputStream(seekIn));
	  master = this.readHeader(fin, chr1, chr2, c1pos1, c1pos2, c2pos1, c2pos2, chr1ind, chr2ind);
	  System.err.println("master="+master);

	  // from header have size of chromosomes, set region to read
	  int c1=Math.min(chr1ind[0],chr2ind[0]);
	  int c2=Math.max(chr1ind[0],chr2ind[0]);
	  int origRegionIndices[]=new int[4]; // as given by user
	  int regionIndices[]=new int[4]; // used to find the blocks we need to access
	  // reverse order if necessary
	  if (chr1ind[0] > chr2ind[0]) {
	    origRegionIndices[0] = c2pos1[0];
	    origRegionIndices[1] = c2pos2[0];
	    origRegionIndices[2] = c1pos1[0];
	    origRegionIndices[3] = c1pos2[0];
	    regionIndices[0] = c2pos1[0] / binsize;
	    regionIndices[1] = c2pos2[0] / binsize;
	    regionIndices[2] = c1pos1[0] / binsize;
	    regionIndices[3] = c1pos2[0] / binsize;
	  }
	  else {
	    origRegionIndices[0] = c1pos1[0];
	    origRegionIndices[1] = c1pos2[0];
	    origRegionIndices[2] = c2pos1[0];
	    origRegionIndices[3] = c2pos2[0];
	    regionIndices[0] = c1pos1[0] / binsize;
	    regionIndices[1] = c1pos2[0] / binsize;
	    regionIndices[2] = c2pos1[0] / binsize;
	    regionIndices[3] = c2pos2[0] / binsize;
	  }

	  IndexEntry c1NormEntry=new IndexEntry(), c2NormEntry=new IndexEntry();
	  long myFilePos[]=new long[] {0};

	  //long bytes_to_read = total_bytes - master;
	  System.err.println("master:"+master);
	
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
	    System.err.println("seek1 "+c1NormEntry.position);
	    seekIn.seek(c1NormEntry.position);
	    seekIn.read(buffer3);
	    
	    ByteArrayInputStream bufferin = new ByteArrayInputStream(buffer3);
	    c1Norm = readNormalizationVector(bufferin);

	    byte buffer4[] = new byte[c2NormEntry.size];
	    System.err.println("seek2 "+c2NormEntry.position);
	    seekIn.seek(c2NormEntry.position);
	    seekIn.read(buffer4);
	    
	  
	    
	    bufferin = new ByteArrayInputStream(buffer4);
	    c2Norm = readNormalizationVector(bufferin);
	  }

	  int blockBinCount[]=new int[] {0};
	  int blockColumnCount[]=new int[] {0};
	  
	  // readMatrix will assign blockBinCount and blockColumnCount
	  readMatrix(seekIn, myFilePos[0], unit, binsize, blockBinCount, blockColumnCount); 
	  
	  Set<Integer> blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount[0], blockColumnCount[0], c1==c2); 

	  // getBlockIndices
	  List<XYValue> xyvalues = new ArrayList<>();
	  for (Integer it:blockNumbers) {
	    // get contacts in this block
		final  List<ContactRecord> records = readBlock(seekIn, null, isHttp, it);
	    for (ContactRecord rec:records) {	      
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
		//printf("%d\t%d\t%.14g\n", x, y, c);
	      }
	    }
	  }
	      //      free(chunk.memory);      
	      /* always cleanup */
	      // curl_easy_cleanup(curl);
	      //    curl_global_cleanup();
    return xyvalues;
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
		  List<XYValue> xyvalues = app.straw(norm, fname, binsize, chr1loc, chr2loc, unit);
		  for (XYValue xyv:xyvalues) {
		    System.out.printf("%d\t%d\t%.14g\n", xyv.x, xyv.y, xyv.value);   
		}
	}	catch(Throwable err) {
		err.printStackTrace();
	}
}

}
