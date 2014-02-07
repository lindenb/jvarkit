package com.github.lindenb.jvarkit.tools.redon;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.broad.tribble.readers.LineIterator;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.Function;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;


class ListOfDouble implements Iterable<Double>,Cloneable
	{
	private double _array[];
	private int _len=0;
	private int _extend=-1;
	ListOfDouble(ListOfDouble copy)
		{
		this._len=copy._len;
		this._array=Arrays.copyOf(copy._array, copy._len);
		this._extend=copy._extend;
		}
	ListOfDouble(int capacity)
		{
		this(capacity,-1);
		}
	ListOfDouble(int capacity,int extend)
		{
		this._array=new double[capacity];
		this._extend=extend;
		Arrays.fill(this._array, 0);
		}
	
	public void push_back(double v)
		{
		if(_len>=this._array.length)
			{
			int x=this._extend;
			if(x==-1) x=1+(_array.length)/2;
			this._array=Arrays.copyOf(_array, _len+x);
			}	
		this._array[_len]=v;
		++_len;
		}
	public int size()
		{
		return _len;
		}
	public double get(int index)
		{
		if(index<0 || index >=_len) throw new IndexOutOfBoundsException("0<="+index+"<"+_len);
		return _array[index];
		}
	public double set(int index,double v)
		{
		if(index<0 || index >=_len) throw new IndexOutOfBoundsException("0<="+index+"<"+_len);
		double old= _array[index];
		_array[index]=v;
		return old;
		}
	@Override
	public Iterator<Double> iterator() {
		return new MyIterator();
		}
	@Override
	public ListOfDouble clone()
		{
		return new ListOfDouble(this);
		}
	public double calc(ProbabilityDistribution algo)
		{
		return algo.calc(this._array,0,this._len);
		}		
	
	public double getMin()
		{
		double minDepth=Double.MAX_VALUE;
		for(int i=0;i< _len;++i)
			{
			double v=_array[i];
			if(Double.isNaN(v)) continue;
			minDepth=Math.min(minDepth,v);
			}
		return minDepth;
		}
	
	public void substract(double v)
		{
		for(int i=0;i< this._len;++i)
			{
			this._array[i]+=v;
			}
		}
	
	public void divide(double v)
		{
		if(v==0) throw new IllegalArgumentException();
		for(int i=0;i< this._len;++i)
			{
			this._array[i]/=v;
			}
		}

	
	private class MyIterator implements Iterator<Double>
		{
		private int index=0;
		@Override
		public boolean hasNext() {
			return index < ListOfDouble.this.size();
			}
		@Override
		public Double next()
			{
			return ListOfDouble.this.get(index++);
			}
		
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}
	
	}
abstract class  ProbabilityDistribution 
	{
	public <IN >double calc(Collection<IN> L,Function<IN, Double> converter)
		{
		double v[]=new double[L.size()];
		int i=0;
		for(IN o:L)
			{
			v[i++]=converter.apply(o);
			}
		return calc(v);
		}
	public double calc(Collection<Double> L)
		{
		double v[]=new double[L.size()];
		int i=0;
		for(Double d:L)
			{
			v[i++]=d;
			}
		return calc(v);
		}
	public double calc(double v[])
		{
		return calc(v,0,v.length);
		}
	public abstract double calc(double v[],int start,int len);
	}

/* use a MEAN */
class  Mean extends ProbabilityDistribution 
	{
	public double calc(double v[],int beg, int len)
		{
		if(len==0) throw new IllegalArgumentException();
		double sum=0;
		for(int i=0;i< len;++i) sum+=v[beg+i];
		return sum/len;
		}
	}
/* use the median  */
class  Median extends ProbabilityDistribution 
	{
	public double calc(double v[],int beg, int len)
		{
		if(len==0) throw new IllegalArgumentException();
		if(len==1) return v[beg];
		double copy[]=new double[len];
		System.arraycopy(v, beg, copy, 0, len);
		Arrays.sort(copy);
		int mid=copy.length/2;
		if(copy.length%2==1)
			{
			return copy[mid];
			}
		else
			{
			return (copy[mid-1]+copy[mid])/2.0;
			}
		}
	}


/**
 * 
 * CopyNumber01
 *
 */
public class CopyNumber01 extends AbstractCommandLineProgram
	{
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private int windowSize=50;
	private int windowStep=25;
	private List<Replicate> replicates=new ArrayList<Replicate>();
	private Set<String> sexualChromosomes=new HashSet<String>();
	private List<List<RegionCaptured>> chrom2capture=null;
	private int num_windows=0;
	private Map<String,String> resolveChromName=new HashMap<String, String>();

	

	private static class GCAndDepth
		implements Comparable<GCAndDepth>
		{
		double depth=0;
		double gc=0;
		@Override
		public int compareTo(final GCAndDepth o)
			{
			if(depth < o.depth) return -1;
			if(depth > o.depth) return  1;

			if(gc < o.gc) return -1;
			if(gc > o.gc) return  1;
			return 0;
			}
		
		}

	
	
	/** a Sample, can contain multiple BAM */
	private class Replicate
		{
		int id=0;
		String filename;
		File tmpFile;
		File createTmpFile(String prefix) throws IOException
			{
			return CopyNumber01.this.createTempFile("pred"+(id+1)+"_"+prefix+"_");
			}
		}
	
	/** A bed segment from the catpure */
	private class RegionCaptured
		implements Comparable<RegionCaptured>,Iterable<RegionCaptured.SlidingWindow>
		{
		int tid;
		int start0;
		int end0;
		
		public SAMSequenceRecord getSAMSequenceRecord()
			{
			return indexedFastaSequenceFile.getSequenceDictionary().getSequence(this.tid);
			}
		
		public String getChromosome()
			{
			return getSAMSequenceRecord().getSequenceName();
			}
		
		public int getStart()
			{
			return this.start0;
			}
		public int getEnd()
			{
			return this.end0;
			}
		@SuppressWarnings("unused")
		public int length()
			{
			return this.getEnd()-this.getStart();
			}
		
		@Override
		public int compareTo(RegionCaptured o)
			{
			int i=tid-o.tid;
			if(i!=0) return i;
			i=start0-o.start0;
			if(i!=0) return i;
			i=end0-o.end0;
			if(i!=0) return i;
			return 0;
			}
		
		Iterator<SlidingWindow> windows()
			{
			return new MyIter();
			}
		
		@Override
		public Iterator<SlidingWindow> iterator()
			{
			return windows();
			}
		
		private class MyIter
			implements Iterator<SlidingWindow>
			{
			int index_in_roi=0;
			private SlidingWindow make()
				{
				return new SlidingWindow(index_in_roi);
				}
			@Override
			public boolean hasNext()
				{
				return make().isValid();
				}
			@Override
			public SlidingWindow next()
				{
				SlidingWindow w=make();
				index_in_roi++;
				return w;
				}
			
			@Override
			public void remove()
				{
				throw new UnsupportedOperationException();
				}
			}
		
		/** a Sliding window from the RegionCaptured */
		public class SlidingWindow
			{
			int index_in_roi;
			private SlidingWindow(int index_in_roi)
				{
				this.index_in_roi=index_in_roi;
				}
			public String getChromosome()
				{
				return RegionCaptured.this.getChromosome();
				}
			public int getStart()
				{
				return index_in_roi*windowStep + RegionCaptured.this.start0;
				}
			public int getEnd()
				{
				return Math.min(
					getStart()+windowSize ,
					RegionCaptured.this.getEnd()
					);
				}
			public int length()
				{
				return getEnd()-getStart();
				}
			boolean isValid()
				{
				if(getStart()+windowSize<   RegionCaptured.this.getEnd())
					{
					return true;
					}
				/** avoid side effect */
				if(getStart()+windowSize >= RegionCaptured.this.getSAMSequenceRecord().getSequenceLength())
					{
					return false;
					}
				int overhang=getStart()+windowSize-RegionCaptured.this.getEnd();
				return overhang>windowSize/2;
				}
			}
		}
	
	
	/** constructor */
	private CopyNumber01()
		{
		}
	
	private File createTempFile(String prefix)
			throws IOException
		{
		File tmpFile=File.createTempFile(prefix, ".data.gz", getTmpDirectories().get(0));
		info("TmpFile"+tmpFile);
		tmpFile.deleteOnExit();
		return tmpFile;
		}
	
	private DataOutputStream openDataStreamForWriting(File f)
		throws IOException
		{
		info("Opening "+f);
		DataOutputStream daos=new DataOutputStream(
				new BufferedOutputStream(
						new GZIPOutputStream(new FileOutputStream(f))));
		return daos;
		}
	
	private DataInputStream openDataStreamForReading(File f)
			throws IOException
		{
		info("Opening "+f);
		IoUtil.assertFileIsReadable(f);
		return new DataInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(f))));
		}
	
	private SAMSequenceDictionary getDicionary()
		{
		return this.indexedFastaSequenceFile.getSequenceDictionary();
		}
	
	/** get a GC% , ignore the sexual chromosomes */
	private ListOfDouble createGCPercent() throws Exception
		{
		ListOfDouble gcPercent=new ListOfDouble(this.num_windows);
		for(List<RegionCaptured> rois:this.chrom2capture)
			{
			if(rois.isEmpty()) continue;
			final String chrom=rois.get(0).getChromosome();
			if(this.sexualChromosomes.contains(chrom))
				{
				info("Ignoring chromosome "+chrom);
				continue;
				}
			GenomicSequence genomic=new GenomicSequence(
					this.indexedFastaSequenceFile,
					chrom
					);
			info("gc% for "+chrom);
			
			for(RegionCaptured roi:rois)
				{
				for(RegionCaptured.SlidingWindow win: roi)
					{
					double total=0f;
					int countN=0;
					for(int pos=win.getStart();pos<win.getEnd() && countN==0;++pos)
						{
						switch(genomic.charAt(pos))
							{
							case 'c':case 'C':
							case 'g':case 'G':		
							case 's':case 'S':
								{
								total++;
								break;
								}
							case 'n':case 'N':countN++;break;
							default:break;
							}
						}
					gcPercent.push_back(countN>0?
							Double.NaN:
							total/(double)win.length()
							);
					}
				}
			}
		info("GC% done size:"+gcPercent.size());
		return gcPercent;
		}
	
	
	private void scan(Replicate replicate,final ListOfDouble gcPercentArray) throws Exception
		{
		List<SAMFileReader> samFileReaders=new ArrayList<SAMFileReader>();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(getDicionary());
		/** open all BAM for this replicate */
		for(String samfile:replicate.filename.split("[\\:]"))
			{
			if(samfile.isEmpty()) continue;
			File bamFile=new File(samfile);
			info("opening "+bamFile);
			SAMFileReader samFileReader=new SAMFileReader(bamFile);
			samFileReaders.add(samFileReader);
			samFileReader.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=samFileReader.getFileHeader();
			
			/* check same dictionaries */
			if(!SequenceUtil.areSequenceDictionariesEqual(
					header.getSequenceDictionary(),
					getDicionary()))
				{
				warning("not the same sequence dictionaries : "+bamFile+"/REFERENCE");
				}
			}

		if(samFileReaders.isEmpty())
			{
			throw new PicardException("empty input for replicate \""+replicate.filename+"\"");
			}
		
		ListOfDouble depth0List=new ListOfDouble(this.num_windows);
		Iterator<Double> gcIter=gcPercentArray.iterator();
		List<GCAndDepth> gcAndDepth=new ArrayList<GCAndDepth>(this.num_windows);
		for(List<RegionCaptured> rois:this.chrom2capture)
			{
			if(rois.isEmpty()) continue;
			final String chrom=rois.get(0).getChromosome();
			boolean is_sexual_chromosome = this.sexualChromosomes.contains(chrom);
			
			info("Alloc for "+chrom);
			short wholeDepth[]=new short[getDicionary().getSequence(chrom).getSequenceLength()];
			Arrays.fill(wholeDepth, (short)0);
			info("Scanning whole "+chrom);
			for(SAMFileReader sfr:samFileReaders)
				{
				String samChromName=chrom;
				
				/* try to resolve chromosome name */
				if(sfr.getFileHeader().getSequenceDictionary().getSequence(samChromName)==null)
					{
					info("chromosome "+samChromName+" is not present in dictionary");
					String altName=resolveChromName.get(samChromName);
					if(altName==null)
						{
						throw new PicardException("in SamFile: unknown chromosome \""+samChromName+"\"");
						}
					if(sfr.getFileHeader().getSequenceDictionary().getSequence(altName)==null)
						{
						throw new PicardException("in SamFile: unknown chromosome \""+samChromName+" or "+altName+ "\"");
						}
					info(samChromName+ " Resolved as "+altName);
					samChromName=altName;
					}
				
				SAMRecordIterator iter=sfr.query(
						samChromName,
						0,//whole chrom
						0,//whole chrom
						true
						);
				while(iter.hasNext())
					{
					SAMRecord rec=iter.next();
					progress.watch(rec);
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getNotPrimaryAlignmentFlag()) continue;
					if(rec.getSupplementaryAlignmentFlag()) continue;
					if(rec.getDuplicateReadFlag()) continue;//check this ?
					
					for(int pos1= rec.getAlignmentStart();
							pos1<=rec.getAlignmentEnd() && pos1<=wholeDepth.length ;
							++pos1
							)
						{
						
						int pos0=pos1-1;
						if(wholeDepth[pos0]<Short.MAX_VALUE)
							{
							wholeDepth[pos0]++;
							}
						else
							{
							warning("High depth "+chrom+":"+pos0);
							}
						}
					}
				iter.close();
				}
			info("End scanning whole chromosome. "+chrom);
			
			for(RegionCaptured roi:rois)
				{
				
				for(RegionCaptured.SlidingWindow win: roi)
					{
					double depth=0;
					
					for(int pos0=win.getStart();
							pos0<win.getEnd() ;
							++pos0
							)
						{
						depth+=wholeDepth[pos0];
						}
					
					double mean_depth=(float)depth/(float)win.length();
					depth0List.push_back(mean_depth);
					
					
					if(!is_sexual_chromosome)
						{
						double G=gcIter.next();
						
						if(!Double.isNaN(G))
							{
							GCAndDepth pair=new GCAndDepth();
							pair.depth=(float)mean_depth;
							pair.gc=G;
							gcAndDepth.add(pair);
							}
						}

					}
				}
			}
		info("Closing Sam File for replicate "+replicate.filename);
		progress.finish();
		CloserUtil.close(samFileReaders);
		if(gcIter.hasNext()) throw new IllegalStateException();

		
		info("sorting gc/depth N="+ gcAndDepth.size());
		Collections.sort(gcAndDepth);
		info("sort done");
			
			{
			ProbabilityDistribution pDist=new Median();
			int j=0;
			int  i1=0;
			while(i1< gcAndDepth.size())
				{
				
				double i1_depth=gcAndDepth.get(i1).depth;
				ListOfDouble all_gc=new ListOfDouble(100);
				all_gc.push_back(gcAndDepth.get(i1).gc);
				int i2=i1+1;
				while(	i2 < gcAndDepth.size() &&
						gcAndDepth.get(i2).depth == i1_depth)
					{
					
					all_gc.push_back(gcAndDepth.get(i2).gc);
					i2++;
					}
				i1=i2;
				
				gcAndDepth.get(j).depth=(float)i1_depth;
				gcAndDepth.get(j).gc=all_gc.calc(pDist);
				j++;
				}
			while(gcAndDepth.size()>j) gcAndDepth.remove(gcAndDepth.size()-1);
			}
			
		info("compacted same depth N="+ gcAndDepth.size() );
		double xval[]=new double[gcAndDepth.size()];
		double yval[]=new double[gcAndDepth.size()];
		
		for(int i=0;i<gcAndDepth.size();++i)
			{
			xval[i]=gcAndDepth.get(i).depth;
			yval[i]=gcAndDepth.get(i).gc;
			}
		
		PolynomialSplineFunction spline=null;
		try
			{
			LoessInterpolator loessInterpolator=new LoessInterpolator();
			info("interpolate");
			spline=loessInterpolator.interpolate(xval,yval);
			}
		catch(org.apache.commons.math3.exception.NumberIsTooSmallException err)
			{
			spline=null;
			warning(err);
			}
		
		info("apply loess");
		for(int i=0;i< depth0List.size();++i)
            {
            if(spline!=null)
            	{
            	try
	            	{
	            	double smoothed=spline.value(depth0List.get(i));
	            	depth0List.set(i,smoothed);
	            	}
            	catch(Exception err)
            		{
            		warning(err);
            		}
            	}
            }
		spline=null;
		xval=null;
		yval=null;
		
		
		
		info("get min depth");
		double minDepth=depth0List.getMin();
		info("substract depth :"+minDepth);
		depth0List.substract(minDepth);
		
		
		
		
		for(int i=0;i< depth0List.size();++i)
			{
			double v=depth0List.get(i);
			depth0List.set(i, v-minDepth);
			}
		
		
		info("get median");
		double median=depth0List.calc(new Median());

		if(median!=0)
			{
			info("div median:="+median);
			depth0List.divide(median);
			}
		replicate.tmpFile=replicate.createTmpFile("depth");
		DataOutputStream daos=openDataStreamForWriting(replicate.tmpFile);
		for(int i=0;i< depth0List.size();++i)
			{
			daos.writeDouble(depth0List.get(i));
			}
		daos.flush();
		daos.close();
		if(depth0List.size()!=this.num_windows) throw new IllegalStateException();
		info("filesize: "+replicate.tmpFile.length());
		}
	
	
	
	private class OuputRow
		{
		RegionCaptured.SlidingWindow win;
		double all_samples[];
		OuputRow(RegionCaptured.SlidingWindow win)
			{
			this.win=win;
			this.all_samples=new double[replicates.size()];
			Arrays.fill(this.all_samples, 0);
			}
		
		}
	private static final int SMOOTH_WINDOW=5;
	/** print smoothing values with neighbours */
	
	
	private void digestAll() throws Exception
		{
		System.out.print("#chrom");
		System.out.print('\t');
		System.out.print("start");
		System.out.print('\t');
		System.out.print("en");
		for(int i=0;i< this.replicates.size();++i)
			{
			System.out.print('\t');
			System.out.print(this.replicates.get(i).filename);
			}
		System.out.println();

		
		List<DataInputStream> disL=new ArrayList<>(this.replicates.size()); 
		for(int i=0;i< this.replicates.size();++i)
			{
			info("Open "+this.replicates.get(i).tmpFile);
			DataInputStream dis=openDataStreamForReading(this.replicates.get(i).tmpFile);
			disL.add(dis);
			}
		
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(getDicionary());
		double all_samples[]=new double[this.replicates.size()];
		ProbabilityDistribution calcMedian=new Median();
		ProbabilityDistribution smoothSamples=new Median();
		for(List<RegionCaptured> rois:this.chrom2capture)
			{
			if(rois.isEmpty()) continue;
			String chrom=rois.get(0).getChromosome();
			info("Saving "+chrom+" rois:"+rois.size());
			for(RegionCaptured roi:rois)
				{
				List<OuputRow> rowBuffer=new ArrayList<OuputRow>();
				for(RegionCaptured.SlidingWindow win:roi)
					{
					progress.watch(chrom, win.getStart());
					for(int i=0;i< this.replicates.size();++i)
						{
						all_samples[i]=disL.get(i).readDouble();
						}
					double median=calcMedian.calc(all_samples);
					if(median>0.2)
						{
						OuputRow newrow=new OuputRow(win);
						rowBuffer.add(newrow);
						for(int i=0;i< this.replicates.size();++i)
							{
							newrow.all_samples[i]=(all_samples[i]/median);
							}
						}
					}
				info("N rows: "+rowBuffer.size());
				for(int rowIndex=0;rowIndex< rowBuffer.size();++rowIndex)
					{
					OuputRow row=rowBuffer.get(rowIndex);
					System.out.print(row.win.getChromosome());
					System.out.print('\t');
					System.out.print(row.win.getStart());
					System.out.print('\t');
					System.out.print(row.win.getEnd());
					
					for(int i=0;i<replicates.size();++i)
						{
						System.out.print('\t');
						ListOfDouble array=new ListOfDouble(SMOOTH_WINDOW*3);
						for(int x= Math.max(0,rowIndex-SMOOTH_WINDOW);
								x<=(rowIndex+SMOOTH_WINDOW) && x< rowBuffer.size();
								++x)
							{
							array.push_back(rowBuffer.get(x).all_samples[i]);
							}
						System.out.print(array.calc(smoothSamples));
						}
					System.out.println();
					if(System.out.checkError()) break;
					}
				}
			info("Done "+chrom);
			}
		CloserUtil.close(disL);
		}
	
	@Override
	public String getProgramDescription() {
		return "CopyNumber01";
		}
	
	 @Override
	protected String getOnlineDocUrl() {
		return "";
	 	}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -R (fasta) reference file indexed with samtools. Required");
		out.println(" -T (dir) add tmp directory.");
		out.println(" -B (file) BED capture file (optional)");
		out.println(" -w (window size) default:"+this.windowSize);
		out.println(" -s (window shft) default:"+this.windowStep);
		out.println(" -X (chrom) add this sexual chromosome.");
		out.println(" -N (file) chrom name helper (name1)(tab2)(name2).");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File bedFile=null;
		File refFile=null;
		String chromNameFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"T:R:w:s:B:X:N:"))!=-1)
			{
			switch(c)
				{
				case 'w': this.windowSize=Integer.parseInt(opt.getOptArg());break;
				case 'B': bedFile=new File(opt.getOptArg());break;
				case 's': this.windowStep=Integer.parseInt(opt.getOptArg());break;
				case 'R': refFile=new File(opt.getOptArg());break;
				case 'T': this.addTmpDirectory(new File(opt.getOptArg()));break;
				case 'X': this.sexualChromosomes.add(opt.getOptArg());break;
				case 'N':
					{
					chromNameFile=opt.getOptArg();
					break;
					}
				default:
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(this.windowSize<=0)
			{
			error("Bad window size.");
			return -1;
			}
		if(this.windowStep<=0)
			{
			error("Bad window step.");
			return -1;
			}
		
		if(refFile==null)
			{
			error("Undefined REF file");
			return -1;
			}
		
		
		if(opt.getOptInd()==args.length)
			{
			error("Illegal Number of arguments.");
			return -1;
			}
		try
			{
			info("Loading "+refFile);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(refFile);
			SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
			this.chrom2capture=new ArrayList<List<RegionCaptured>>(dict.size());
			while(this.chrom2capture.size() < dict.size())
				{
				this.chrom2capture.add(new ArrayList<CopyNumber01.RegionCaptured>());
				}
			
			if(chromNameFile!=null)
				{
				info("Reading "+chromNameFile);
				LineIterator r=IOUtils.openURIForLineIterator(chromNameFile);
				while(r.hasNext())
					{
					String tokens[]=r.next().split("[\t]");
					if(tokens.length<2) continue;
					this.resolveChromName.put(tokens[0], tokens[1]);
					this.resolveChromName.put(tokens[1], tokens[1]);
					}
				CloserUtil.close(r);
				}
			
			
			
			if(bedFile!=null)
				{
				int N=0;
				info("Reading BED:" +bedFile);
				Pattern tab=Pattern.compile("[\t]");
				LineIterator r=IOUtils.openFileForLineIterator(bedFile);
				while(r.hasNext())
					{
					String line=r.next();
					if(line.startsWith("#") || line.isEmpty()) continue;
					String tokens[]=tab.split(line,4);
					if(tokens.length<3)
						{
						warning("No enough column in "+line+" "+bedFile);
						continue;
						}
					RegionCaptured roi=new RegionCaptured();
					roi.tid=dict.getSequenceIndex(tokens[0]);
					if(roi.tid==-1)
						{
						String altName=this.resolveChromName.get(tokens[0]);
						if(altName!=null) roi.tid=dict.getSequenceIndex(altName);
						}
					
					
					if(roi.tid==-1)
						{
						warning("not in reference: chromosome "+tokens[0]+" in "+line+" "+bedFile);
						continue;
						}
					roi.start0=Integer.parseInt(tokens[1]);
					roi.end0=Integer.parseInt(tokens[2]);
					this.chrom2capture.get(roi.tid).add(roi);
					++N;
					}
				CloserUtil.close(r);
				info("end Reading BED:" +bedFile+"  N="+N);
				for(List<RegionCaptured> roi:this.chrom2capture)
					{
					Collections.sort(roi);
					}
				}
			else
				{
				info("No capture, peeking everything");
				for(SAMSequenceRecord ssr:dict.getSequences())
					{
					RegionCaptured roi=new RegionCaptured();
					roi.tid=ssr.getSequenceIndex();
					roi.start0=0;
					roi.end0=ssr.getSequenceLength();
					this.chrom2capture.get(roi.tid).add(roi);
					}
				}
			
			this.num_windows=0;
			for(List<RegionCaptured> rois:this.chrom2capture)
				{
				for(RegionCaptured roi:rois)
					{
					for(@SuppressWarnings("unused") RegionCaptured.SlidingWindow win:roi)
						{
						this.num_windows++;
						}
					}
				}
			info("Number of windows : "+this.num_windows);
			
			ListOfDouble gcPercentArray=createGCPercent();
			for(int optind=opt.getOptInd();optind< args.length;++optind)
				{
				Replicate replicate=new Replicate();
				replicate.id=this.replicates.size();
				replicate.filename=args[optind];
				info("Read replicate: "+args[optind]+" : "+(1+optind-opt.getOptInd())+"/"+(args.length-opt.getOptInd()));
				this.replicates.add(replicate);
				scan(replicate,gcPercentArray);
				}
			gcPercentArray=null;
			
			digestAll();
		
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}	
		}
	
	public static void main(String[] args) {
		new CopyNumber01().instanceMainWithExit(args);
		}
	}