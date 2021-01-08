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


A large part of this code is based on https://github.com/aidenlab/Juicebox (Aiden Lab)

The MIT License (MIT)

Copyright (c) 2019 Aiden Lab

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
package com.github.lindenb.jvarkit.hic;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.zip.InflaterInputStream;

import com.github.lindenb.jvarkit.lang.Paranoid;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalParserFactory;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.util.LittleEndianInputStream;

/**
HIC format: 
	https://github.com/aidenlab/Juicebox/blob/5a56089c63957cb15401ea7906ab77e242dfd755/HiCFormatV8.md
	https://github.com/aidenlab/straw/blob/24a50c4777e8992270402fa4465cd5f7d1dad8ca/C%2B%2B/straw.cpp
    https://www.encodeproject.org/files/ENCFF784GFP/@@download/ENCFF784GFP.hic <- doesn't work
    https://github.com/aidenlab/Juicebox/blob/5a56089c63957cb15401ea7906ab77e242dfd755/HiC_format_v8.docx
	https://github.com/igvteam/hic-straw/blob/d428ee7e6df5488dd1295b33a81eeb06adbf7a51/src/hicFile.js#L81
*/
public class HicReaderImpl implements HicReader {
	private final Paranoid paranoid = Paranoid.createThrowingInstance();
	private static final Logger LOG = Logger.build(HicReaderImpl.class).make();
	private static final String MAGIC="HIC";
	private static final int DEFAULT_VERSION = 8;
	
	private static final boolean DEBUG = true || System.getProperty("jvarkit.hic.debug","").equals("true");
	private static final boolean SKIP = false;
	
	
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
	
	private static class ContactRecord
		{
		final int binX;
		final int binY;
		final float counts;
		ContactRecord(int binX,int binY, float counts) {
			this.binX = binX;
			this.binY =  binY;
			this.counts =  counts;
			}
		}
	
	private LittleEndianInputStream streamToEndian() throws IOException  {
		return new LittleEndianInputStream(new BufferedInputStream(this.seekableStream));
		}
	
	private void debug(Object o) {
		if(!DEBUG) return;
		LOG.debug(o);
		}
	
	private void fullySkip(LittleEndianInputStream in,long n) throws IOException {
		long n2 = in.skip(n);
		if(n2!=n) throw new IOException("Cannot skip "+n+" bytes (got "+n2+")");
		}
	
	/** called by HicReaderFactory */
	HicReaderImpl(final Object source,final SeekableStream seekableStream) throws IOException {
		this.source = source;
		this.seekableStream = seekableStream;
		
		@SuppressWarnings("resource")
		LittleEndianInputStream lis = this.streamToEndian();
		final String magic = lis.readString();
		
		/* parse magic */
		if(!magic.equals(MAGIC)) {
			this.seekableStream.close();
			throw new IOException("Bad magic for " + source);
			}
		
		/* parse version */
		this.version = lis.readInt();
		
		if(this.version != DEFAULT_VERSION) {
			this.seekableStream.close();
			throw new IOException("Version "+this.version+" not supported. for " + source+" expected "+DEFAULT_VERSION);
			}
		
		/* master index */
		this.masterIndexPosition = lis.readLong();
		paranoid.assertGt(this.masterIndexPosition, 0);

		/* genome name */
		this.genomeId = lis.readString();

		
		/* attributes */
		final int n_attributes = lis.readInt();

		
		final Map<String,String> atts = new LinkedHashMap<>(n_attributes);
		for(int i=0;i< n_attributes;i++) {
			final String key = lis.readString();
			
			final String value = lis.readString();
			
			atts.put(key, value);
			}
		this.attributes = Collections.unmodifiableMap(atts);
		
		/* dictionary */
		final int n_chromosomes = lis.readInt();
		paranoid.assertGt(n_chromosomes, 0);

		
		final List<SAMSequenceRecord> ssrList = new ArrayList<>(n_chromosomes);
		for(int i=0;i< n_chromosomes;i++) {
			final String chromName = lis.readString();

			
			final int chromLen = lis.readInt();
			paranoid.assertGt(chromLen, 0);
			
			final SAMSequenceRecord ssr = new SAMSequenceRecord(chromName, chromLen);
			ssr.setAssembly(this.genomeId);
			ssrList.add(ssr);
			}
		this.dictionary = new SAMSequenceDictionary(ssrList);

		
		/* number of base pair resolutions */
		final int nBpResolution = lis.readInt();
		paranoid.assertGe(nBpResolution, 0);
		
		final Set<Integer> resBPSet = new TreeSet<>();
		for(int i=0;i< nBpResolution;i++) {
			final int resBP = lis.readInt();
			resBPSet.add(resBP);
			}
		this.basePairResolutions = Collections.unmodifiableSet(resBPSet);
		
		/* number of fragment resolutions */
		final int nFragResolution = lis.readInt();
		paranoid.assertGe(nBpResolution, 0);
		
		final Set<Integer> resFragSet = new TreeSet<>();
		for(int i=0;i< nFragResolution;i++) {
			final int resF = lis.readInt();
			resFragSet.add(resF);
			}
		this.fragmentResolutions = Collections.unmodifiableSet(resFragSet);

		// 
		
		
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
	public Optional<Locatable> parseInterval(final String s) {
		if(StringUtils.isBlank(s)) return Optional.empty();
		final Function<String, Optional<SimpleInterval>> intervalParser =
				IntervalParserFactory.newInstance().
				dictionary(this.getDictionary()).
				enableWholeContig().
				make();
		final Optional<SimpleInterval> r= intervalParser.apply(s);
		return r.isPresent()?
				Optional.of(r.get()):
				Optional.empty()
				;
		}

	private QueryInterval convertLocatableToQueryInterval(final Locatable loc) {
		final SAMSequenceRecord ssr = getDictionary().getSequence(loc.getContig());
		return ssr==null?null: new QueryInterval(ssr.getSequenceIndex(),loc.getStart(),loc.getEnd());
		}
	
	
	@Override
	public boolean query(
			final Locatable interval1,
			final Locatable interval2,
			final Normalization norm,
			final int binsize, 
			final Unit unit,
			final HicReader.QueryCallBack callback
			)
		{
		try {
			if(callback==null) throw new IllegalArgumentException("callback is null");
			final Query q = new Query();
			q.interval1 = interval1;
			q.interval2 = interval2;
			q.normalization = norm;
			q.unit = unit;
			q.binsize = binsize;
			q.callback = callback;
			
			final Function<Locatable, QueryInterval> interval2query = (R)->{
				final QueryInterval q1= convertLocatableToQueryInterval(R);
				if(q1==null) {
					q.callback.error("unknown contig in \""+R+"\". Available are: "+
						getDictionary().getSequences().stream().map(SSR->SSR.getSequenceName()).collect(Collectors.joining(";")));
					}
				return q1;
				};
			
			q.qInterval1 = interval2query.apply(interval1);
			if(q.qInterval1==null) return false;
			
			q.qInterval2 = interval2query.apply(interval2);
			if(q.qInterval2==null) return false;
			
			/* swap if needed */
			if(q.qInterval1.referenceIndex < q.qInterval2.referenceIndex ||
				(q.qInterval1.referenceIndex == q.qInterval2.referenceIndex && q.qInterval1.start > q.qInterval2.start)
				) {
				debug("swap "+q.qInterval1 +" "+q.qInterval2);
				final Locatable tmp1 = q.interval2;
				q.interval2 = q.interval1;
				q.interval1 = tmp1;
				
				final QueryInterval tmp2 = q.qInterval2;
				q.qInterval2 = q.qInterval1;
				q.qInterval1 = tmp2;
				}
			
			
			q.scanFooter();
			if(q.chr_chri_fpos<0L) {
				q.callback.warning("cannot find chri_chrj_fpos");
				return false;
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
			
		   q.readMatrix(q.chr_chri_fpos); 
			
		  final Set<Integer> blockNumbers = q.getBlockNumbersForRegionFromBinPosition(); 
			
			
			
		 
		  
		  for (final Integer it:blockNumbers) {
		    // get contacts in this block
		    for(final ContactRecord rec:q.readBlockId(it)) {     
		      final int x = rec.binX * binsize;
		      final int y = rec.binY * binsize;
		      
		      if(!CoordMath.encloses(q.qInterval1.start, q.qInterval1.end, x, x)) continue;
		      if(!CoordMath.encloses(q.qInterval2.start, q.qInterval2.end, y, y)) continue;
		      
		      if(q.qInterval1.referenceIndex==q.qInterval2.referenceIndex) {
			      if(!CoordMath.encloses(q.qInterval1.start, q.qInterval1.end, y, y)) continue;
			      if(!CoordMath.encloses(q.qInterval2.start, q.qInterval2.end, x, x)) continue;
		      	}
		      
		      float c = rec.counts;
		      if (!norm .equals(Normalization.NONE)) {
				  c = (float)(c / (c1Norm[rec.binX] * c2Norm[rec.binY]));
			      }
		     
		     q.callback.reportContact(
		    		 q.interval1.getContig(), x, x+binsize,
		    		 q.interval2.getContig(), y, y+binsize, norm, unit, binsize, c);
		    }
		  }
			return false;
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
    /*
	public void scan(
			final Locatable interval,
			final Normalization norm,
			final int binsize, 
			final Unit unit
			)
		{
		try {
			final ScanRegion q = new ScanRegion();
			q.interval1 = interval;
			q.normalization = norm;
			q.unit = unit;
			q.binsize = binsize;
			
			q.qInterval1 = convertLocatableToQueryInterval(interval);
			
			
			q.scanFooter();
			if(q.chri_chrj_fpos.isEmpty()) {
				return  ;
				}
			final Map<Integer,double[]> tid2norm;
			
			if (!q.normalization.equals(Normalization.NONE)) {
				tid2norm = new HashMap<>(q.tid2indexentry.size());
				for(final Integer tid: q.tid2indexentry.keySet()) {
					tid2norm.put(tid, readNormalizationVector(q.tid2indexentry.get(tid)));
					}
			    }
			else
				{
				tid2norm = null;
				}
		for(final Long offset: q.chri_chrj_fpos) {
			   q.readMatrix(offset); 
			   final Set<Integer> blockNumbers = q.getBlockNumbersForRegionFromBinPosition(); 
	
			   
			 
				
				
				
				
			  // getBlockIndices
			  final HicMatrixImpl matrixImpl = new HicMatrixImpl();
			  matrixImpl.interval1 = q.interval1;
			  matrixImpl.interval2 = q.interval2;
			  matrixImpl.binSize = q.binsize;
			  matrixImpl.normalization = q.normalization;
			  matrixImpl.unit = q.unit;
			  
			  for (final Integer it:blockNumbers) {
			    // get contacts in this block
			    for(final ContactRecord rec:q.readBlock(it)) {     
			      final int x = rec.binX * binsize;
			      final int y = rec.binY * binsize;
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
				
			    	  final HicContact xyv = new HicContactImpl(x,y,c);
			    	 
			    	  matrixImpl.contacts.add(xyv);
			      }
			    }
			  }
			  }
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}*/

	private static class IndexEntry {
		  final int size;
		  final long position;
		  IndexEntry(final int size,final long position)  {
			  this.size = size;
			  this.position = position;
		  }
		  @Override
		public String toString() {
			return "IndexEntry( size:"+size+" position:"+position+")";
		  }
		}
	
	// reads the normalization vector from the file at the specified location
	private double[] readNormalizationVector(final IndexEntry entry) throws IOException {
		  debug("read normalisation " + entry); 
		  final byte buf[] = new byte[entry.size];
		  seekableStream.seek(entry.position);
		  seekableStream.readFully(buf);

	 	  @SuppressWarnings("resource")
		  final LittleEndianInputStream in = new LittleEndianInputStream(new ByteArrayInputStream(buf));
		  final int nValues = in.readInt();
		  paranoid.assertGe(nValues,0);
		  double values[] = new double[nValues];
		  for (int i = 0; i < nValues; i++) {
		     values[i] = in.readDouble();
		  	 }
		  return values;
		  }
	
	
	
	/**
	https://github.com/igvteam/juicebox.js/blob/55bd6c7815f9abee74368c14a9d9403d2998313f/js/hicDataset.js#L95 	 
	https://github.com/igvteam/hic-straw/blob/d428ee7e6df5488dd1295b33a81eeb06adbf7a51/src/hicFile.js#L307 */
	private List<ContactRecord> readBlock(final IndexEntry indexEntry) throws IOException {
		debug("read block " + indexEntry); 
		if (indexEntry==null) {
			 return Collections.emptyList();
		 	}
		 debug(indexEntry);
		 if (indexEntry.size == 0) {
			 return Collections.emptyList();
		 }
		 final byte compressedBytes[] = new byte[indexEntry.size];

		 seekableStream.seek(indexEntry.position);
		 seekableStream.readFully(compressedBytes);
		 final InflaterInputStream zipIn = new InflaterInputStream(new ByteArrayInputStream(compressedBytes));


		 // create stream from buffer for ease of use

		 @SuppressWarnings("resource")
		 final LittleEndianInputStream bufferin=new LittleEndianInputStream(zipIn);
		 final int nRecords = bufferin.readInt();
		 paranoid.assertGe(nRecords, 0);
		 final List<ContactRecord> contactRecords = new ArrayList<>(nRecords);


		 if(getVersion()<7) {
			 final int binX = bufferin.readInt();
			 final int binY = bufferin.readInt();
			 final float counts = bufferin.readFloat();
             contactRecords.add(new ContactRecord(binX, binY, counts));
		 }

		 
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
				 final ContactRecord record = new ContactRecord(binX,binY,counts);
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
					 if (c != Short.MIN_VALUE) {
						 final ContactRecord record = new ContactRecord(bin1,bin2,c);
						 contactRecords.add(record);
					 }
				 } 
				 else {

					 counts  = bufferin.readFloat();
					 if (!Float.isNaN(counts)) { // not sure this works
						 //	  if (!Float.isNaN(counts)) {
						 final ContactRecord record = new ContactRecord(bin1,bin2,counts);
						 contactRecords.add(record);
					 }
				 }
			 }
			 break;
		     }
		 default: throw new IOException("unknown block type");
		 }
	 return contactRecords;
	 }
	
	/* data to ignore in the header */
	private void skipExpectedValuesMaps(final LittleEndianInputStream fin) throws IOException  {
	  debug("skipExpectedVaues");
	  // read in and ignore expected value maps; don't store; reading these to 
	  // get to norm vector index
	  final int nExpectedValues1 = fin.readInt();
	  paranoid.assertGe(nExpectedValues1, 0);
	  
	  for (int i=0; i< nExpectedValues1; i++) {
	    Unit.valueOf(fin.readString());
	    fin.readInt();//binSize
	    
	    final int nValues = fin.readInt();
	    
	    paranoid.assertGe(nValues, 0);
	    
		if(SKIP) {
			fullySkip(fin,nValues*Double.BYTES);
		    }
	    else
		    {
		    for(int j=0;j< nValues;j++)
		    	{
		    	fin.readDouble();
		    	}
		    }

	    final int nNormalizationFactors = fin.readInt();
	    paranoid.assertGe(nNormalizationFactors, 0);
	    
	    if(SKIP) {
			fullySkip(fin,nNormalizationFactors*(Double.BYTES+Integer.BYTES));
		    }
	    else
		    {
		    for(int j=0;j< nNormalizationFactors;j++)
		    	{
		    	fin.readInt();
		    	fin.readDouble();
		    	}		    
		    }
	    }
	  
	  final int nExpectedValues2 = fin.readInt();
	  paranoid.assertGe(nExpectedValues2, 0);
	  for (int i=0; i<nExpectedValues2; i++) {
	    fin.readString(); //typeString
	    Unit.valueOf(fin.readString()); //unit
	    fin.readInt();//binSize
	    

	    final int nValues = fin.readInt();
	    paranoid.assertGe(nValues, 0);
	    if(SKIP) {
	    	fullySkip(fin,nValues*(Double.BYTES));
	    	}
	    else
		    {
		    for(int j=0;j< nValues;j++)
		    	{
		    	fin.readDouble();
		    	}
		    }
	    
	    final int nNormalizationFactors = fin.readInt();
	    paranoid.assertGe(nNormalizationFactors, 0);
	    if(SKIP) {
			fullySkip(fin,nNormalizationFactors*(Double.BYTES+Integer.BYTES));
	    	} 
	    else
		    {
		    for(int j=0;j< nNormalizationFactors;j++)
		    	{
		    	fin.readInt();//tid
		    	fin.readDouble();//value
		    	}
		    }
	    }
	  debug("skipExpectedVaues: done");
	}

	
	private abstract class AbstractQuery
		{
		/** user query */
		Normalization normalization;
		/** user unit */
		Unit unit;
		/** user interval1 */
		Locatable interval1;
		/** user resolution */
		int binsize;//resolution
		
		
		/** normalized intervals from query */
		QueryInterval qInterval1;
		

		
		

	}
	
	private class Query extends AbstractQuery {
		/** user interval1 */
		Locatable interval2;

		/** normalized intervals from query */
		QueryInterval qInterval2;
				
		HicReader.QueryCallBack callback;
		
		/** chr-chr filepos found in scanHeader */
		long chr_chri_fpos = -1L;

		
		/** set in readFooter */
		IndexEntry normEntry1 = null;
		/** set in readFooter */
		IndexEntry normEntry2 = null;
		
		

		/** zoom data */
		int blockBinCount;
		int blockColumnCount;
		final Map<Integer, IndexEntry> blockMap = new TreeMap<>();

		
		
	
		/* reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count 
		https://github.com/igvteam/hic-straw/blob/d428ee7e6df5488dd1295b33a81eeb06adbf7a51/src/hicFile.js#L727
		*/
		
		
		
		
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
		
		void scanFooter() throws IOException {
		
		
		debug("seek to "+ masterIndexPosition);
		HicReaderImpl.this.seekableStream.seek(HicReaderImpl.this.masterIndexPosition);
		final LittleEndianInputStream fin = HicReaderImpl.this.streamToEndian();
		final int nBytes = fin.readInt();//nBytes Total size, master index + expected values
		paranoid.assertGe(nBytes,0);
		
		// loop over the entries to find the chr-chr data
		final int nEntries1 = fin.readInt();
		paranoid.assertGe(nEntries1, 0);
		
		// loop over master index
		for (int i=0; i<nEntries1; i++) {
		    final String str = fin.readString();
		    final long fpos = fin.readLong();
		    paranoid.assertGe(fpos, 0L);
		    if(SKIP) {
		    	fullySkip(fin,Integer.BYTES);//sizeinbytes
		    	}
		    else
		    	{
		    	fin.readInt();
		    	}
		    
		    final int u  = str.indexOf('_');
		    if(u==-1) throw new IllegalStateException("Cannot find underscore in "+u);
		    final int tid1 = Integer.parseInt(str.substring(0,u));
		    final int tid2 = Integer.parseInt(str.substring(u+1));
		    paranoid.assertLe(tid1, tid2);
		    
		    if (tid1==this.qInterval1.referenceIndex && tid2==this.qInterval2.referenceIndex )
		    	{
		    	paranoid.assertLt(this.chr_chri_fpos,0);
		    	this.chr_chri_fpos = fpos;
		    	// no 'continue', must read whole header
		    	}
		    //no 'break', must read whole header
		    }
		
		 // not found 
		 if ( this.chr_chri_fpos < 0L ) {
			this.callback.warning( "File "+getSource()+
					" doesn't have the given key "+
					this.qInterval1.referenceIndex + "_" +
					this.qInterval2.referenceIndex + " map"
					);
		    return;
		  	}
		  
		  if (Normalization.NONE.equals(this.normalization)) return; // no need to read norm vector index
		  
		  skipExpectedValuesMaps(fin);
		  
		  // Index of normalization vectors
		  final int nEntries2 = fin.readInt();
		  paranoid.assertGe(nEntries2, 0);
		  
		  for (int i = 0; i < nEntries2 && (this.normEntry1==null || this.normEntry2==null); i++) {
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
		   this.callback.warning( "Normalization vectors not found for one or both chromosomes at " + this.binsize + " " + this.unit+" available ");
		  	}
		}
		

		/** https://github.com/igvteam/hic-straw/blob/d428ee7e6df5488dd1295b33a81eeb06adbf7a51/src/hicFile.js#L226 */
		protected void readMatrix(final long offset) throws IOException {
			debug("seek matrix at "+offset);
			seekableStream.seek(offset);  
			final LittleEndianInputStream in = streamToEndian();
			if(SKIP)
				{
				fullySkip(in,Integer.BYTES * 2);//ignore c1 + c2
				}
			else
				{
				in.readInt();//c1
				in.readInt();//c2
				}
			
			//  # of resolution levels (bp and frags)
			final int nResolutions = in.readInt();
			paranoid.assertGe(nResolutions, 0);
			
			for(int i=0; i<nResolutions;i++) {
				if(this.readMatrixZoomData(in)) return;
				}
			throw new IOException("Error finding block data");
			}
		
		private boolean readMatrixZoomData(final LittleEndianInputStream fin) throws IOException 
		  {
		  debug("read zoom data");
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
		  paranoid.assertGe(nBlocks, 0);
		  
		  for (int i = 0; i < nBlocks; i++) {
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
		
		/**
		https://github.com/igvteam/juicebox.js/blob/55bd6c7815f9abee74368c14a9d9403d2998313f/js/hicDataset.js#L95 	 
		https://github.com/igvteam/hic-straw/blob/d428ee7e6df5488dd1295b33a81eeb06adbf7a51/src/hicFile.js#L307 */
		private List<ContactRecord> readBlockId(final int blockNumber) throws IOException {
			return readBlock(this.blockMap.get(blockNumber));
			
		 	}
		}
	
	}
