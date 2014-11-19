/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.redon;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;


/**
 * 
 * CopyNumber01
 *
 */
public class CopyNumber01 extends AbstractCommandLineProgram
	{
	/** sample Name */
	private String sampleName="SAMPLE";
	/** reference */
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;	
	/** chrom name helper */
	private Map<String,String> resolveChromName=new HashMap<String, String>();
	/** global sam dict */
	private SAMSequenceDictionary samDictionary=null;
	/** map interval to depths and GC */
	private List<GCAndDepth> interval2row=new ArrayList<GCAndDepth>(1000);
	/** size of a window */
	private int windowSize=100;
	
	
	

	/* fact: Y=depth X=GC% */
	private class GCAndDepth
		{
		int tid;
		int start;
		int end;
		double depth=0;
		double gc=0;

		double getX()
			{
			return gc;
			}
		double getY()
			{
			return depth;
			}
		
		public long getGenomicIndex()
			{
			long n=this.start;
			for(int i=0;i< this.tid;++i) n+=  samDictionary.getSequence(i).getSequenceLength();
			return n;
			}
		
		public String getChrom()
			{
			return samDictionary.getSequence(this.tid).getSequenceName();
			}
		}

	private static final Comparator<GCAndDepth> sortOnXY=new Comparator<CopyNumber01.GCAndDepth>()
			{
			@Override
			public int compare(GCAndDepth a, GCAndDepth b)
				{
				if(a.getX() < b.getX()) return -1;
				if(a.getX() > b.getX()) return  1;

				if(a.getY() < b.getY()) return -1;
				if(a.getY() > b.getY()) return  1;
				return 0;
				}
			};
	private static final Comparator<GCAndDepth> sortOnPosition=new Comparator<CopyNumber01.GCAndDepth>()
			{
			@Override
			public int compare(GCAndDepth a, GCAndDepth b)
				{
				int i= a.tid - b.tid;
				if(i!=0) return i;
				return a.start - b.start;
				}
			};

			
			
	
	
	/** constructor */
	private CopyNumber01()
		{
		}
		
	
	private boolean ignoreChromosomeName(String chrom)
		{
		return !chrom.matches("(chr)?([0-9]+|X|Y)");
		}
	
	private void prefillGCPercent(
			GenomicSequence genomic,
			int chromStart,
			int chromEnd) throws Exception
			{
			int pos=chromStart;
			while(pos+this.windowSize<chromEnd && pos< genomic.length())
				{
				char c=genomic.charAt(pos);
				if(c=='n' || c=='N')
					{
					++pos;
					continue;
					}
				int total=0;
				int n=0;
				boolean foundN=false;
				for(n=0;n<this.windowSize && pos +n < genomic.length() && !foundN;++n)
					{
					switch(genomic.charAt(pos+n))
						{
						case 'c':case 'C':
						case 'g':case 'G':		
						case 's':case 'S':
							{
							total++;
							break;
							}
						case 'n':case 'N':foundN=true;break;
						default:break;
						}
					}
				if(n != this.windowSize)
					{
					pos++;
					continue;
					}
				GCAndDepth dataRow=new GCAndDepth();
				dataRow.tid=genomic.getSAMSequenceRecord().getSequenceIndex();
				dataRow.start=pos+1;
				dataRow.end=pos+this.windowSize;
				
				dataRow.gc=total/(double)this.windowSize;
				
				this.interval2row.add(dataRow);
				
				pos+=this.windowSize;
				}
			}
	
	/** get a GC% */
	private void prefillGCPercentWithCapture(File bedFile) throws Exception
		{
		Pattern tab=Pattern.compile("[\t]");
		BufferedReader in= IOUtils.openFileForBufferedReading(bedFile);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.trim().isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line,4);
			String chrom=tokens[0];
			if(this.samDictionary.getSequence(chrom)==null)
				{
				chrom = this.resolveChromName.get(chrom);
				if(chrom==null)
					{
					info("Cannot resolve "+chrom);
					continue;
					}
				}
			
			if(ignoreChromosomeName(chrom))
				{
				info("Ignoring "+chrom);
				continue;
				}
			
			
			GenomicSequence genomic=new GenomicSequence(
				this.indexedFastaSequenceFile,
				chrom
				);
			int bedStart=Integer.parseInt(tokens[1]);
			int bedEnd=Integer.parseInt(tokens[2]);
			prefillGCPercent(genomic, bedStart, bedEnd);
			}
		}
	
	
	/** get a GC% */
	private void prefillGCPercentWithoutCapture() throws Exception
		{
		for(SAMSequenceRecord ssr:this.indexedFastaSequenceFile.getSequenceDictionary().getSequences())
			{
			String chrom=ssr.getSequenceName();
			if(this.samDictionary.getSequence(chrom)==null)
				{
				chrom = this.resolveChromName.get(chrom);
				if(chrom==null)
					{
					info("Cannot resolve "+chrom);
					continue;
					}
				}
			
			if(ignoreChromosomeName(chrom))
				{
				info("Ignoring "+ssr.getSequenceName());
				continue;
				}
			
			
			GenomicSequence genomic=new GenomicSequence(
				this.indexedFastaSequenceFile,
				chrom
				);
		
			prefillGCPercent(genomic,0, ssr.getSequenceLength());
			}
		}
	
	private void scanCoverage(SamReader sr)
		throws IOException
		{
	
		for(GCAndDepth row:this.interval2row)
			{
			double sum = 0;
			SAMRecordIterator sri=sr.query(
					this.samDictionary.getSequence(row.tid).getSequenceName(),
					row.start,
					row.end,
					false);
			while(sri.hasNext())
				{
				SAMRecord rec = sri.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag())continue;
				Cigar c= rec.getCigar();
				int refStart= rec.getAlignmentStart();
				for(CigarElement ce:c.getCigarElements())
					{
					if(!ce.getOperator().consumesReferenceBases()) continue;
					if(ce.getOperator().consumesReadBases())
						{
						for(int x=0;x< ce.getLength();++x)
							{
							if(refStart+x>= row.start && refStart+x<=row.end)
								{
								sum++;
								}
							}
						}
					refStart+=ce.getLength();
					}		
				}
			sri.close();
			row.depth += sum /(double)this.windowSize;
			}
		}
	
	private boolean isSexualChrom(String chrom)
		{
		return chrom.matches("(chr?)(X|Y)");
		}
	
	
	private void normalizeCoverage()
		{
		int autosome_count=0;
		Collections.sort(this.interval2row,CopyNumber01.sortOnXY);
		
		for(int j=0;j< this.interval2row.size();++j)
			{
			GCAndDepth r=this.interval2row.get(j);
			if(isSexualChrom(r.getChrom())) continue;
			autosome_count++;
			}

		
		
		double x[]=new double[autosome_count];
		double y[]=new double[autosome_count];

		int i=0;
		for(int j=0;j< this.interval2row.size();++j)
			{
			GCAndDepth r=this.interval2row.get(j);
			if(isSexualChrom(r.getChrom())) continue;
			x[i] = r.getX();
			y[i] = r.getY();
			++i;
			}
		LoessInterpolator interpolator=new LoessInterpolator();
		PolynomialSplineFunction  spline = interpolator.interpolate(x, y);
		
		for(GCAndDepth gc:this.interval2row)
			{
			gc.depth = spline.value(gc.getX());
			}
			
		}
	
	
	private void smoothCoverage()
		{		
		Collections.sort(this.interval2row,CopyNumber01.sortOnPosition);

		double x[]=new double[this.interval2row.size()];
	
		for(int j=0;j< x.length;++j)
			{
			x[j] = this.interval2row.get(j).getX();
			}

		
		
		for(int j=0;j< x.length;++j)
			{
			int curr_tid= this.interval2row.get(j).tid;
			for(int y=j-SMOOTH_WINDOW;y<=j+SMOOTH_WINDOW && y< x.length;++y)
				{
				if(this.interval2row.get(j).tid!=curr_tid) continue;
				}
			x[j] = this.interval2row.get(j).getX();
			}
			
		}

	
	private void saveCoverage(GZIPOutputStream zout)
		{
		info("Dumping coverage ");
		PrintWriter pw=new PrintWriter(zout);
		
		/* header */
		pw.println("ID\tCHROM\tSTART\tEND\tGC\t"+this.sampleName);
		
		/* get data */
		for(GCAndDepth r:this.interval2row)
			{
			pw.print(r.getGenomicIndex());
			pw.print('\t');
			pw.print( this.samDictionary.getSequence(r.tid).getSequenceName());
			pw.print('\t');
			pw.print( r.start);
			pw.print('\t');
			pw.print( r.end);
			pw.print('\t');
			pw.print(r.gc);
			pw.print('\t');
			pw.print(r.depth);
			pw.println();
			}
		pw.flush();
		}
	
	
	private static final int SMOOTH_WINDOW=5;
	/** print smoothing values with neighbours */
	
	
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
		out.println(" -R (fasta) "+getMessageBundle("reference.faidx")+". Required");
		out.println(" -b (file) BED capture file (optional)");
		out.println(" -w (window size) default:"+this.windowSize);
		out.println(" -N (file) chrom name helper (name1)(tab2)(name2).");
		out.println(" -o  output base name.");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		String outfile="output";
		File bedFile=null;
		File refFile=null;
		String chromNameFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:w:b:N:o:"))!=-1)
			{
			switch(c)
				{
				case 'w': this.windowSize=Integer.parseInt(opt.getOptArg());break;
				case 'b': bedFile=new File(opt.getOptArg());break;
				case 'R': refFile=new File(opt.getOptArg());break;
				case 'o': outfile = opt.getOptArg();break;
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
		
		
		if(refFile==null)
			{
			error("Undefined REF file");
			return -1;
			}
		
		
		if(outfile==null)
			{
			error("Undefined output file.");
			return -1;
			}
	
		

		SamReader samReader = null;
		
		try
			{
			File bamFile=null;
			if(opt.getOptInd()+1==args.length)
				{
				bamFile=new File(args[opt.getOptInd()]);
				}
			else
				{
				error("Illegal Number of arguments.");
				return -1;
				}
			
			final SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			info("get Dict for "+bamFile);
			samReader = srf.open(bamFile);
			this.samDictionary=samReader.getFileHeader().getSequenceDictionary();
			for(SAMReadGroupRecord rg:samReader.getFileHeader().getReadGroups())
				{
				String s = rg.getSample();
				if(s==null || s.trim().isEmpty()) continue;
				this.sampleName=s;
				break;
				}
			samReader.close();
			
			
			
			/* loading REF Reference */
			info("Loading "+refFile);
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(refFile);
			SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null)
				{
				error("Cannot get sequence dictionary for "+refFile);
				return -1;
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
					this.resolveChromName.put(tokens[1], tokens[0]);
					}
				CloserUtil.close(r);
				}
			
			if(bedFile!=null)
				{
				prefillGCPercentWithCapture(bedFile);
				}
			else
				{
				prefillGCPercentWithoutCapture();
				}
			
			samReader = srf.open(bamFile);
			scanCoverage(samReader);
			samReader.close();
				
			/* save raw coverage */
			GZIPOutputStream zout = new GZIPOutputStream(new FileOutputStream(outfile+"_raw.tsv.gz"));
			saveCoverage(zout);
			zout.finish();zout.close();
			
			normalizeCoverage();
			
			/* save normalized coverage */
			zout = new GZIPOutputStream(new FileOutputStream(outfile+"_normalized.tsv.gz"));
			saveCoverage(zout);
			zout.finish();zout.close();

		
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
