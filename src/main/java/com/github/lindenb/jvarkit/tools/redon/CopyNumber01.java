package com.github.lindenb.jvarkit.tools.redon;

import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.broad.tribble.readers.LineIterator;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;

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
	private SamSequenceRecordTreeMap<Boolean> capture=null;
	
	private static class GCAndDepth
		implements Comparable<GCAndDepth>
		{
		float depth=0;
		float gc=0;
		@Override
		public int compareTo(final GCAndDepth o)
			{
			if((double)depth < (double)o.depth) return -1;
			if((double)depth > (double)o.depth) return  1;

			if((double)gc < (double)o.gc) return -1;
			if((double)gc > (double)o.gc) return  1;
			return 0;
			}
		}
	
	
	private static class Replicate
		{
		int id=0;
		String filename;
		List<SAMFileReader> samFileReaders=new ArrayList<>();
		File tmpFile;
		}
	private List<Replicate> replicates=new ArrayList<Replicate>();
	
	private CopyNumber01() {
		
		}
	
	private boolean accept(SAMSequenceRecord ssr)
		{
		return (this.capture==null?true:capture.containsChromosome(ssr.getSequenceIndex()));
		}
	
	private boolean accept(SAMSequenceRecord ssr,int pos0)
		{
		return(this.capture==null?true:capture.containsOverlapping(ssr.getSequenceIndex(), pos0+1, pos0+1));
		}
	private void scan(Replicate replicate) throws Exception
		{
		File tmpDir=getTmpDirectories().get(0);
		replicate.tmpFile=File.createTempFile("replicate"+replicate.id+"_", ".data.gz", tmpDir);
		replicate.tmpFile.deleteOnExit();
		info("Writing to "+replicate.tmpFile);
		DataOutputStream daos=new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(replicate.tmpFile))));
		daos.writeUTF(replicate.filename);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.indexedFastaSequenceFile.getSequenceDictionary());
		//step 1, get coverage
		for(SAMSequenceRecord ssr:this.indexedFastaSequenceFile.getSequenceDictionary().getSequences())
			{
			if(!accept(ssr)) continue;
			progress.watch(ssr.getSequenceName(), 1);
			info("step 1: "+ssr.getSequenceName()+" length:"+ssr.getSequenceLength());
			short depth[]=new short[ssr.getSequenceLength()];
			Arrays.fill(depth, (short)0);
			for(SAMFileReader sfr:replicate.samFileReaders)
				{
				SAMRecordIterator iter=null;
				iter=sfr.query(ssr.getSequenceName(), 0/*start*/, 0/*end*/, true);
				while(iter.hasNext())
					{
					SAMRecord rec=iter.next();
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getNotPrimaryAlignmentFlag()) continue;
					if(rec.getSupplementaryAlignmentFlag()) continue;
					if(rec.getDuplicateReadFlag()) continue;//check this ?
					for(int pos1=rec.getAlignmentStart();
							pos1<=rec.getAlignmentEnd() && (pos1-1)< depth.length;
							++pos1)
						{
						if(!accept(ssr,pos1-1)) continue;
						if(depth[pos1-1]<Short.MAX_VALUE)
							{
							depth[pos1-1]++;
							}
						}
					}
				iter.close();
				}
			
			info("end depth for "+ssr.getSequenceName());
			final int num_windows=depth.length/windowStep;
			info("n-windows :"+num_windows+" win:"+windowStep);
			List<GCAndDepth> gcAndDepthList=new ArrayList<>(num_windows);
			double smooth[]=new double[num_windows];
			for(int i=0;i< num_windows;i++)
				{
				double total=0f;
				int pos0=i*windowStep;
				for(int j=0;j< windowSize && pos0+j< depth.length;++j)
					{
					total+=depth[pos0+j];
					}
				GCAndDepth gcdepth=new GCAndDepth();
				gcdepth.depth=(float)(total/windowSize);
				gcAndDepthList.add(gcdepth);
				smooth[i]=gcdepth.depth;
				}
			depth=null;
			
			info("getting gc%");
			GenomicSequence genomic=new GenomicSequence(this.indexedFastaSequenceFile, ssr.getSequenceName());
			for(int i=0;i< gcAndDepthList.size(); i++)
				{
				int pos0=i*windowStep;
				if(!accept(ssr,pos0)) continue;
				GCAndDepth gcAndDepth=gcAndDepthList.get(i);
				/* if(gcAndDepthList.get(i).depth<=0) continue; */
				double total=0f;
				
				for(int j=0;j< windowSize && pos0+j< genomic.length();++j)
					{
					switch(genomic.charAt(pos0+j))
						{
						case 'c':case 'C':
						case 'g':case 'G':		
						case 's':case 'S':
							{
							total++;
							break;
							}
						default:break;
						}
					}
				gcAndDepth.gc=(float)total/windowSize;
				}
			
			int j=0;
			for(int i=0;i<  gcAndDepthList.size();++i)
				{
				if(!accept(ssr,i*windowStep)) continue;
				gcAndDepthList.set(j,gcAndDepthList.get(i));
				j++;
				}
			while(gcAndDepthList.size()>j) gcAndDepthList.remove(gcAndDepthList.size()-1);

			info("sorting gc/depth N="+ gcAndDepthList.size());
			Collections.sort(gcAndDepthList);
			info("sort done");
				
				{
				j=0;
				int  i1=0;
				while(i1< gcAndDepthList.size())
					{
					double i1_depth=gcAndDepthList.get(i1).depth;
					double total=gcAndDepthList.get(i1).gc;
					int n=1;
					int i2=i1+1;
					while(	i2 < gcAndDepthList.size() &&
							(double)gcAndDepthList.get(i2).depth == i1_depth)
						{
						total+= gcAndDepthList.get(i2).gc;
						++n;
						i2++;
						}
					i1=i2;
					gcAndDepthList.get(j).depth=(float)i1_depth;
					gcAndDepthList.get(j).gc=(float)(total/n);
					j++;
					}
				while(gcAndDepthList.size()>j) gcAndDepthList.remove(gcAndDepthList.size()-1);
				}
				
			info("compacted same depth N="+ gcAndDepthList.size() );
			double xval[]=new double[gcAndDepthList.size()];
			double yval[]=new double[gcAndDepthList.size()];
			
			for(int i=0;i<gcAndDepthList.size();++i)
				{
				xval[i]=gcAndDepthList.get(i).depth;
				yval[i]=gcAndDepthList.get(i).gc;
				}
			try
				{
				LoessInterpolator loessInterpolator=new LoessInterpolator();
				info("calcal loess");
				PolynomialSplineFunction spline=loessInterpolator.interpolate(xval,yval);
				info("end loess");
				for(int i=0;i< smooth.length;++i)
	                {
					if(!accept(ssr,i*windowStep)) continue;
	                smooth[i]=spline.value(smooth[i]);
	                }
				spline=null;
				}
			catch(org.apache.commons.math3.exception.NumberIsTooSmallException err)
				{
				warning(err);
				}
			xval=null;
			yval=null;
			

			info("get min depth");
			double minDepth=Double.MAX_VALUE;
			for(int i=0;i< smooth.length;++i)
				{
				if(!accept(ssr,i*windowStep)) continue;
				minDepth=Math.min(minDepth, smooth[i]);
				}
			
			info("substract depth :"+minDepth);
			List<Double> median_array=new ArrayList<Double>();
			for(int i=0;i< smooth.length;++i)
				{
				if(!accept(ssr,i*windowStep)) continue;
				smooth[i]-=minDepth;
				
				median_array.add(smooth[i]);
				}
			
			info("get median");
			Collections.sort(median_array);
			double median=median_array.get(median_array.size()/2);
			median_array=null;
			
			if(median!=0)
				{
				info("div median:="+median);
				for(int i=0;i< smooth.length;++i)
					{
					if(!accept(ssr,i*windowStep)) continue;
					smooth[i]/=median;
					}
				}
			
			info("Writing "+ssr.getSequenceName());
			//write tid
			daos.writeInt(ssr.getSequenceIndex());
			daos.writeInt(smooth.length);
			for(int i=0;i< smooth.length;++i)
				{
				daos.writeDouble(smooth[i]);
				}
			daos.flush();
			}
		daos.flush();
		daos.close();
		progress.finish();
		info("filesize: "+replicate.tmpFile.length());
		}
	
	private class OuputRow
		{
		int index;
		double all_samples[];
		OuputRow(int index)
			{
			this.index=index;
			this.all_samples=new double[replicates.size()];
			Arrays.fill(this.all_samples, 0);
			}
		
		}
	private static final int SMOOTH_WINDOW=5;
	/** print smoothing values with neighbours */
	void print(final SAMSequenceRecord ssr,final List<OuputRow> rows,int rowIndex)
		{
		OuputRow row=rows.get(rowIndex);
		System.out.print(ssr.getSequenceName());
		System.out.print('\t');
		System.out.print(row.index*windowStep);
		System.out.print('\t');
		System.out.print((row.index*windowStep)+windowSize);
		for(int i=0;i<replicates.size();++i)
			{
			System.out.print('\t');
			List<Double> array=new ArrayList<Double>(SMOOTH_WINDOW*2+1);
			for(int x=Math.max(0,rowIndex-SMOOTH_WINDOW);x<=(rowIndex+SMOOTH_WINDOW) && x< rows.size();++x)
				{
				array.add(rows.get(x).all_samples[i]);
				}
			Collections.sort(array);
			System.out.print(array.get(array.size()/2));
			}
		System.out.println();
		}
	
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
			DataInputStream dis=new DataInputStream(new GZIPInputStream(new FileInputStream(this.replicates.get(i).tmpFile)));
			String name=dis.readUTF();
			info("name is "+name);
			disL.add(dis);
			}
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.indexedFastaSequenceFile.getSequenceDictionary());
		for(SAMSequenceRecord ssr:this.indexedFastaSequenceFile.getSequenceDictionary().getSequences())
			{
			if(!accept(ssr)) continue;  
			List<OuputRow> rowBuffer=new ArrayList<OuputRow>(SMOOTH_WINDOW*3);
			final int num_windows=ssr.getSequenceLength()/windowStep;
			double normalize[]=new double[num_windows];
			Arrays.fill(normalize, 0.0);
			
			for(int i=0;i< this.replicates.size();++i)
				{
				int tid=disL.get(i).readInt();
				if(tid!=ssr.getSequenceIndex()) throw new IllegalStateException();
				if(disL.get(i).readInt()!=num_windows)  throw new IllegalStateException();
				}
			double all_samples[]=new double[this.replicates.size()];
			double copy[]=new double[all_samples.length];
			for(int j=0;j< num_windows;++j)
				{
				//int pos0=j*windowStep;
				progress.watch(ssr.getSequenceName(), j*windowStep);
				
				for(int i=0;i< this.replicates.size();++i)
					{
					all_samples[i]=disL.get(i).readDouble();
					}
				
				if(!accept(ssr,j*windowStep))
					{
					//dump buffer
					for(int x=0;x<rowBuffer.size();++x) print(ssr, rowBuffer, x);
					rowBuffer.clear();
					if(System.out.checkError()) break;
					continue;
					}
				
				System.arraycopy(all_samples, 0, copy, 0, all_samples.length);
				Arrays.sort(copy);
				double median=copy[copy.length/2];
				if(median>0.2)
					{
					OuputRow newrow=new OuputRow(j);
					rowBuffer.add(newrow);
					for(int i=0;i< this.replicates.size();++i)
						{
						newrow.all_samples[i]=(all_samples[i]/median);
						}
					}
				if(System.out.checkError()) break;
				}
			//dump buffer
			for(int x=0;x<rowBuffer.size();++x) print(ssr, rowBuffer, x);
			rowBuffer.clear();
			if(System.out.checkError()) break;
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
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File bedFile=null;
		File refFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"T:R:w:s:B:"))!=-1)
			{
			switch(c)
				{
				case 'w': this.windowSize=Integer.parseInt(opt.getOptArg());break;
				case 'B': bedFile=new File(opt.getOptArg());break;
				case 's': this.windowStep=Integer.parseInt(opt.getOptArg());break;
				case 'R':refFile=new File(opt.getOptArg());break;
				case 'T':this.addTmpDirectory(new File(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
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
			
			if(bedFile!=null)
				{
				int N=0;
				info("Reading BED:" +bedFile);
				this.capture=new SamSequenceRecordTreeMap<>(this.indexedFastaSequenceFile.getSequenceDictionary());
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
					if(this.indexedFastaSequenceFile.getSequenceDictionary().getSequence(tokens[0])==null)
						{
						warning("not in reference: chromosome "+tokens[0]+" in "+line+" "+bedFile);
						continue;
						}
					this.capture.put(tokens[0],Integer.parseInt(tokens[1])+1, Integer.parseInt(tokens[2]),Boolean.TRUE);
					++N;
					}
				CloserUtil.close(r);
				info("end Reading BED:" +bedFile+"  N="+N);
				}
			
			for(int optind=opt.getOptInd();optind< args.length;++optind)
				{
				Replicate replicate=new Replicate();
				replicate.id=this.replicates.size();
				replicate.filename=args[optind];
				info("Read replicate: "+args[optind]+" : "+(1+optind-opt.getOptInd())+"/"+(args.length-opt.getOptInd()));
				this.replicates.add(replicate);
				
				for(String samfile:replicate.filename.split("[\\:]"))
					{
					if(samfile.isEmpty()) continue;
					 File bamFile=new File(samfile);
					info("opening "+bamFile);
					SAMFileReader samFileReader=new SAMFileReader(bamFile);
					replicate.samFileReaders.add(samFileReader);
					samFileReader.setValidationStringency(ValidationStringency.SILENT);
					SAMFileHeader header=samFileReader.getFileHeader();
					
					if(!SequenceUtil.areSequenceDictionariesEqual(
							header.getSequenceDictionary(),
							indexedFastaSequenceFile.getSequenceDictionary()))
						{
						error("not the same sequence dictionaries : "+bamFile+"/"+refFile);
						return -1;
						}
					}
				
				if(replicate.samFileReaders.isEmpty())
					{
					error("empty input");
					return -1;
					}
				scan(replicate);
				info("Closing replicate");
				CloserUtil.close(replicate.samFileReaders);
				}
			
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
			for(Replicate rep:this.replicates) CloserUtil.close(rep.samFileReaders);
			}	
		}
	
	public static void main(String[] args) {
		new CopyNumber01().instanceMainWithExit(args);
		}
	}
