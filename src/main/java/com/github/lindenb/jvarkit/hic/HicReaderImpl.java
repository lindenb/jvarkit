package com.github.lindenb.jvarkit.hic;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.InflaterInputStream;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.util.LittleEndianInputStream;

/**
HIC format: 
	https://github.com/aidenlab/Juicebox/blob/5a56089c63957cb15401ea7906ab77e242dfd755/HiCFormatV8.md
	https://github.com/aidenlab/straw/blob/24a50c4777e8992270402fa4465cd5f7d1dad8ca/C%2B%2B/straw.cpp
    https://www.encodeproject.org/files/ENCFF784GFP/@@download/ENCFF784GFP.hic <- doesn't work
    https://github.com/aidenlab/Juicebox/blob/5a56089c63957cb15401ea7906ab77e242dfd755/HiC_format_v8.docx
*/
public class HicReaderImpl implements HicReader {
	private static final Logger LOG = Logger.build(HicReaderImpl.class).make();
	private static final String MAGIC="HIC";
	private static final int DEFAULT_VERSION = 8;
	final Object source;
	private final SeekableStream seekableStream;
	/** Version number */
	private final int version;
	/** File position of master index **/
	private final long masterIndexPosition;
	/** Genome identifier **/
	private final String genomeId;
	/** Genome identifier **/
	private final Map<String,String> attributes;
	/** dict **/
	private final SAMSequenceDictionary dictionary;
	/**  base pair resolutions */
	private final Set<Integer> basePairResolutions;
	/** fragment resolutions */
	private final Set<Integer> fragmentResolutions;
	
	static class XYValue
		{
		int x;
		int y;
		float value;
		}
	
	 static class ContactRecord
		{
		int binX;
		int binY;
		float counts;
		}
	
	private LittleEndianInputStream streamToEndian() throws IOException  {
		return new LittleEndianInputStream(new BufferedInputStream(this.seekableStream));
		}
	
	HicReaderImpl(final Path path) throws IOException  {
		this(path.toString());
	}
	
	HicReaderImpl(final String path) throws IOException  {
		this(path,SeekableStreamFactory.getInstance().getStreamFor(path));
		}
	
	HicReaderImpl(final Object source,final SeekableStream seekableStream) throws IOException {
		this.source = source;
		this.seekableStream = seekableStream;
		
		@SuppressWarnings("resource")
		LittleEndianInputStream lis = this.streamToEndian();
		final String magic = lis.readString();
		if(!magic.equals(MAGIC)) {
			this.seekableStream.close();
			throw new IOException("Bad magic for " + source);
			}
		
		this.version = lis.readInt();
		
		if(this.version != DEFAULT_VERSION) {
			this.seekableStream.close();
			throw new IOException("Version "+this.version+" not supported. for " + source+" expected "+DEFAULT_VERSION);
			}
		
		this.masterIndexPosition = lis.readLong();
		
		this.genomeId = lis.readString();
		
		
		/* attributes */
		final int n_attributes = lis.readInt();
		
		final Map<String,String> atts = new HashMap<String, String>(n_attributes);
		for(int i=0;i< n_attributes;i++) {
			final String key = lis.readString();
			final String value = lis.readString();			
			atts.put(key, value);
			}
		this.attributes = Collections.unmodifiableMap(atts);
		
		final int n_chromosomes = lis.readInt();
		
		final List<SAMSequenceRecord> ssrList = new ArrayList<>(n_chromosomes);
		for(int i=0;i< n_chromosomes;i++) {
			final String chromName = lis.readString();
			
			final int chromLen = lis.readInt();
			
			final SAMSequenceRecord ssr = new SAMSequenceRecord(chromName, chromLen);
			ssr.setAssembly(this.genomeId);
			ssrList.add(ssr);
			}
		this.dictionary = new SAMSequenceDictionary(ssrList);
		
		/* number of base pair resolutions */
		final int nBpResolution = lis.readInt();
		final Set<Integer> resBPSet = new HashSet<>(nBpResolution);
		for(int i=0;i< nBpResolution;i++) {
			final int resBP = lis.readInt();
			resBPSet.add(resBP);
			}
		this.basePairResolutions = Collections.unmodifiableSet(resBPSet);
		
		/* number of fragment resolutions */
		final int nFragResolution = lis.readInt();
		final Set<Integer> resFragSet = new HashSet<>(nFragResolution);
		for(int i=0;i< nFragResolution;i++) {
			final int resF = lis.readInt();
			resFragSet.add(resF);
			}
		this.fragmentResolutions = Collections.unmodifiableSet(resFragSet);
		
		}
		
	
	@Override
	public Set<Integer> getBasePairResolutions() {
		return this.basePairResolutions;
		}
	
	@Override
	public Set<Integer> getFragmentResolutions() {
		return this.fragmentResolutions;
		}
	
	@Override
	public void close() {
		CloserUtil.close(this.seekableStream);
		}
	
	
	@Override
	public SAMSequenceDictionary getDictionary() {
		return this.dictionary;
		}
	
	@Override
	public String getBuild() {
		return this.genomeId;
		}
	
	@Override
	public Object getSource() {
		return this.source;
		}
	
	@Override
	public Map<String, String> getAttributes() {
		return this.attributes;
		}
	
	@Override
	public int getVersion() {
		return this.version;
		}
	
	
	
	@Override
	public Iterator<QueryResult> query(
			final Locatable interval1,
			final Locatable interval2,
			final Normalization norm,
			final int binsize, 
			final Unit unit
			)
		{
		try {
			final Query q= new Query();
			q.interval1 = interval1;
			q.interval2 = interval2;
			q.normalization = norm;
			q.unit = unit;
			q.binsize = binsize;
			
			SAMSequenceRecord ssr = getDictionary().getSequence(interval1.getContig());
			if(ssr==null) {
				q.errorMessage = "unknown contig in "+interval1;
				return Collections.emptyIterator();
				}
			q.qInterval1 = new QueryInterval(ssr.getSequenceIndex(),interval1.getStart(),interval1.getEnd());
			
			ssr = getDictionary().getSequence(interval2.getContig());
			if(ssr==null) {
				q.errorMessage = "unknown contig in "+interval2;
				return Collections.emptyIterator();
				}
			q.qInterval2 = new QueryInterval(ssr.getSequenceIndex(),interval2.getStart(),interval2.getEnd());

			/* swap if needed */
			if(q.qInterval1.referenceIndex < q.qInterval2.referenceIndex) {
				Locatable tmp1 = q.interval2;
				q.interval2 = q.interval1;
				q.interval1 = tmp1;
				
				QueryInterval tmp2 = q.qInterval2;
				q.qInterval2 = q.qInterval1;
				q.qInterval1 = tmp2;

				}
			
			
			q.scanFooter();
			if(q.chr_chri_fpos<0L) {
				return Collections.emptyIterator();
				}
			final double c1Norm[];
			final double c2Norm[];
			
			if (!q.normalization.equals(Normalization.NONE)) {
			    c1Norm = readNormalizationVector(q.normEntry1);
			    c2Norm = readNormalizationVector(q.normEntry2);
			    }
			else
				{
				c1Norm = null;
				c2Norm = null;
				}
			
		   q.readMatrix(); 
			
		  final Set<Integer> blockNumbers = q.getBlockNumbersForRegionFromBinPosition(); 
			
			
			
		  // getBlockIndices
		  final List<XYValue> xyvalues = new ArrayList<>();
		  for (Integer it:blockNumbers) {
		    // get contacts in this block
		    for(ContactRecord rec:q.readBlock(it)) {     
		      int x = rec.binX * binsize;
		      int y = rec.binY * binsize;
		      float c = rec.counts;
		      if (!norm .equals(Normalization.NONE)) {
				  c = (float)(c / (c1Norm[rec.binX] * c2Norm[rec.binY]));
			      }

		      if ((x >= q.qInterval1.start && x <= q.qInterval1.end &&
			   y >= q.qInterval2.start && y <= q.qInterval2.end) ||
			  // or check regions that overlap with lower left
			  ((q.qInterval1.referenceIndex==q.qInterval2.referenceIndex) && 
				 y >= q.qInterval1.start && y <= q.qInterval1.end && 
				 x >= q.qInterval2.start && x <= q.qInterval2.end)) {
			
		    	  XYValue xyv = new XYValue();
		    	  xyv.x=x;
		    	  xyv.y=y;
		    	  xyv.value=c;
		    	  xyvalues.add(xyv);
		      }
		    }
		  }

			return Collections.emptyIterator();
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}

	private static class IndexEntry {
		  final int size;
		  final long position;
		  IndexEntry(final int size,final long position)  {
			  this.size = size;
			  this.position = position;
		  }
		}
	
	// reads the normalization vector from the file at the specified location
	private double[] readNormalizationVector(final IndexEntry entry) throws IOException {
	      final byte buf[] = new byte[entry.size];
		  seekableStream.seek(entry.position);
		  seekableStream.readFully(buf);

	 	  @SuppressWarnings("resource")
		  final LittleEndianInputStream in = new LittleEndianInputStream(new ByteArrayInputStream(buf));
		  final int nValues = in.readInt();
		  double values[] = new double[nValues];
		  for (int i = 0; i < nValues; i++) {
		     values[i]=in.readDouble();
		  	 }
		  return values;
		  }
	
	
	
	
	private class Query {
		/** user query */
		Normalization normalization;
		Unit unit;
		Locatable interval1;
		Locatable interval2;
		int binsize;//resolution
		
		/** normalized intervals from query */
		QueryInterval qInterval1;
		QueryInterval qInterval2;
		
		/** error message */
		String errorMessage = null;
		
		
		/** set in readFooter */
		IndexEntry normEntry1 = null;
		/** set in readFooter */
		IndexEntry normEntry2 = null;
		
		/** chr-chr filepos found in scanHeader */
		long chr_chri_fpos = -1L;
		
		/** zoom data */
		int blockBinCount;
		int blockColumnCount;
		final Map<Integer, IndexEntry> blockMap = new TreeMap<>();


		private void readMatrix() throws IOException {
			seekableStream.seek(this.chr_chri_fpos);  
			LittleEndianInputStream in = streamToEndian();
			  in.skip(Integer.BYTES * 2);//ignore c1 + c2
			  final int nRes = in.readInt();
			  int i=0;
			  while (i<nRes) {
			    final boolean found = this.readMatrixZoomData(in);
			    if(found) return;
			    i++;
			    }
			throw new IOException("Error finding block data");
			}
		
		
	
		// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count 
		private boolean readMatrixZoomData(final LittleEndianInputStream fin) throws IOException 
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
		  //debug("readMatrixZoomData "+unit +" "+binSize+" "+blockBinCount+" "+blockColumnCount);
		  
		  boolean storeData=false;
		  if (this.unit.equals(unit) && this.binsize == binSize) {
			this.blockBinCount = blockBinCount;
			this.blockColumnCount = blockColumnCount;
			storeData = true;
		    }
		  
		  
		  final int nBlocks  = fin.readInt();
	
		  for (int b = 0; b < nBlocks; b++) {
		    final int blockNumber = fin.readInt();
		    final long filePosition = fin.readLong();
		    final int blockSizeInBytes = fin.readInt();
		    final IndexEntry entry = new IndexEntry(blockSizeInBytes,filePosition);
		    if (storeData) {
		    	this.blockMap.put(blockNumber, entry);
		    	}
		  	}
		  return storeData;
		}	
		
		
		
		private Set<Integer> getBlockNumbersForRegionFromBinPosition( ) {
			   final int regionIndices[]= {
					  this.qInterval1.start/this.binsize,
					  this.qInterval1.end/this.binsize,
					  this.qInterval2.start/this.binsize,
					  this.qInterval2.end/this.binsize
				};
			
			
		       
			   final int col1 = regionIndices[0] / this.blockBinCount;
			   final int col2 = (regionIndices[1] + 1) / this.blockBinCount;
			   final int row1 = regionIndices[2] / this.blockBinCount;
			   final int row2 = (regionIndices[3] + 1) / this.blockBinCount;
			   
			   final Set<Integer> blocksSet = new TreeSet<>();
			   // first check the upper triangular matrix
			   for (int r = row1; r <= row2; r++) {
			     for (int c = col1; c <= col2; c++) {
			       final int blockNumber = r * this.blockColumnCount + c;
			       blocksSet.add(blockNumber);
			     }
			   }
			   // check region part that overlaps with lower left triangle
			   // but only if intrachromosomal
			   if (this.qInterval1.referenceIndex==this.qInterval2.referenceIndex) {
			     for (int r = col1; r <= col2; r++) {
			       for (int c = row1; c <= row2; c++) {
				     final int blockNumber = r * this.blockColumnCount + c;
				     blocksSet.add(blockNumber);
			       }
			     }
			   }

			return blocksSet;
			}
		
		private void scanFooter() throws IOException {
		
		final String key = getDictionary().getSequence(this.qInterval1.referenceIndex).getSequenceName()+
			   "_" +
			   getDictionary().getSequence(this.qInterval2.referenceIndex).getSequenceName()
			   ;
		
		HicReaderImpl.this.seekableStream.seek(HicReaderImpl.this.masterIndexPosition);
		final LittleEndianInputStream fin = HicReaderImpl.this.streamToEndian();
		final int nBytes = fin.readInt();

		// loop over the entries to find the chr-chr data
		  int nEntries = fin.readInt();
		  for (int i=0; i<nEntries; i++) {
		    final String str = fin.readString();
		    final long fpos = fin.readLong();
		    fin.skip(Integer.SIZE);//sizeinbytes
		    if (!str.equals(key)) continue;
		    this.chr_chri_fpos = fpos;
		    break;
		    }
		 // not found 
		 if ( this.chr_chri_fpos <0L ) {
			this.errorMessage = "File "+getSource()+" doesn't have the given "+key+" map" ;
		    return;
		  	}
		  
		  if (Normalization.NONE.equals(this.normalization)) return; // no need to read norm vector index
		  
		  // read in and ignore expected value maps; don't store; reading these to 
		  // get to norm vector index
		  int nExpectedValues = fin.readInt();
		  
		  for (int i=0; i<nExpectedValues; i++) {
		    fin.readString();//ignore unit
		    fin.skip(Integer.BYTES );//binSize
		    
		    final int nValues = fin.readInt();
		    fin.skip(Double.BYTES * nValues);

		    final int nNormalizationFactors = fin.readInt();
		    fin.skip((Integer.BYTES /* chrIdx */ + Double.BYTES /* value */) * nNormalizationFactors);
		    }
		 
		  
		  nExpectedValues = fin.readInt();
		  
		  for (int i=0; i<nExpectedValues; i++) {
		    fin.readString(); //typeString
		    fin.readString(); //unit
		    fin.skip(Integer.BYTES );//binSize
		    

		    final int nValues = fin.readInt();
		    fin.skip(Double.BYTES * nValues);
		    
		    final int nNormalizationFactors = fin.readInt();
		    fin.skip((Integer.BYTES /* chrIdx */ + Double.BYTES /* value */) * nNormalizationFactors);
		    }
		  
		  // Index of normalization vectors
		  nEntries = fin.readInt();
		  
		  
		  for (int i = 0; i < nEntries && (this.normEntry1==null || this.normEntry2==null); i++) {
		    final Normalization normtype = Normalization.valueOf(fin.readString());
		    final int chrIdx = fin.readInt();
		    final Unit unit1 = Unit.valueOf(fin.readString());
		    final int resolution1 = fin.readInt();
		    
		    
		    final long filePosition = fin.readLong();
		    final int sizeInBytes= fin.readInt();
		    if (chrIdx == this.qInterval1.referenceIndex && normtype.equals(this.normalization) && unit1.equals(this.unit) && resolution1 == this.binsize) {
		      this.normEntry1 = new IndexEntry(sizeInBytes, filePosition);
		    }
		    if (chrIdx == this.qInterval2.referenceIndex && normtype.equals(this.normalization) && unit1.equals(this.unit) && resolution1 == this.binsize) {
		      this.normEntry2 = new IndexEntry(sizeInBytes, filePosition);
		    }
		  }
		  if (this.normEntry1==null || this.normEntry2==null) {
		   this.errorMessage =  "Normalization vectors not found for one or both chromosomes at " + this.binsize + " " + this.unit+" available ";
		  	}
		}
		
		
		
		 private List<ContactRecord> readBlock(final int blockNumber) throws IOException {
			 IndexEntry idx = this.blockMap.get(blockNumber);
			 if (idx==null) {
				 return Collections.emptyList();
			 }

			 if (idx.size == 0) {
				 return Collections.emptyList();
			 }
			 final byte compressedBytes[] = new byte[idx.size];

			 seekableStream.seek(idx.position);
			 seekableStream.readFully(compressedBytes);
			 final InflaterInputStream zipIn = new InflaterInputStream(new ByteArrayInputStream(compressedBytes));


			 // create stream from buffer for ease of use

			 @SuppressWarnings("resource")
			 final LittleEndianInputStream bufferin=new LittleEndianInputStream(zipIn);
			 final int nRecords = bufferin.readInt();
			 final List<ContactRecord> contactRecords = new ArrayList<>(nRecords);

			 int binXOffset =bufferin.readInt();
			 int binYOffset =bufferin.readInt();
			 byte useShort = bufferin.readByte();
			 int type = (int)bufferin.readByte();

			 switch(type) {
			 case 1: {
				 // List-of-rows representation
				 short rowCount = bufferin.readShort();
				 for (int i = 0; i < rowCount; i++) {
					 short y = bufferin.readShort();
					 int binY = y + binYOffset;
					 short colCount = bufferin.readShort();
					 for (int j = 0; j < colCount; j++) {
						 short x = bufferin.readShort();
						 int binX = binXOffset + x;
						 float counts;
						 if (useShort == 0) { // yes this is opposite of usual
							 short c = bufferin.readShort();
							 counts = c;
						 } 
						 else {
							 counts = bufferin.readFloat();
						 }
						 ContactRecord record = new ContactRecord();
						 record.binX = binX;
						 record.binY = binY;
						 record.counts = counts;
						 contactRecords.add(record);
					 }
				 }
				 break;
			 }
			 case 2 : { // have yet to find test file where this is true, possibly entirely deprecated
				 int nPts = bufferin.readInt();
				 short w = bufferin.readShort();

				 for (int i = 0; i < nPts; i++) {
					 //int idx = (p.y - binOffset2) * w + (p.x - binOffset1);
					 int row = i / w;
					 int col = i - row * w;
					 int bin1 = binXOffset + col;
					 int bin2 = binYOffset + row;

					 float counts;
					 if (useShort == 0) { // yes this is opposite of the usual
						 short c = bufferin.readShort();
						 if (c != -32768) {
							 ContactRecord record = new ContactRecord();
							 record.binX = bin1;
							 record.binY = bin2;
							 record.counts = c;
							 contactRecords.add(record);
						 }
					 } 
					 else {

						 counts  = bufferin.readFloat();
						 if (counts != 0x7fc00000) { // not sure this works
							 //	  if (!Float.isNaN(counts)) {
							 ContactRecord record = new ContactRecord();
							 record.binX = bin1;
							 record.binY = bin2;
							 record.counts = counts;
							 contactRecords.add(record);
						 }
					 }
				 }
				 break;
			     }

			 }

		return contactRecords;
		 }
	}
	
	}
