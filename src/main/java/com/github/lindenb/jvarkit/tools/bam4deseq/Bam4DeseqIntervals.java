package com.github.lindenb.jvarkit.tools.bam4deseq;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;


import com.github.lindenb.jvarkit.util.picard.cmdline.CommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;

public class Bam4DeseqIntervals extends CommandLineProgram
	{
	private static final Log LOG=Log.getInstance(Bam4DeseqIntervals.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" create a table for DESEQ with the number of reads within a sliding window for multiple BAMS.";
    
	@Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process",optional=false,minElements=1)
	public List<File> IN=new ArrayList<File>();
	@Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="output filename. Default stdout. ",optional=true)
	public File OUT=null;
	@Option(shortName="WIN",doc="size of the observed window. ",optional=true)
	public int WINDOW_SIZE=500;
	@Option(shortName="SHIFT",doc="shift window by SHIFT pb ",optional=true)
	public int WINDOW_SHIFT=250;
	@Option(shortName="COV",doc="ignore regions with NO coverage",optional=true)
	public boolean ONLY_COVERED=false;
	@Option(shortName="HEAD",doc="print header",optional=true)
	public boolean HEADER=true;
	@Option(shortName="CHR", doc="limit to this chromosome",optional=true,minElements=0)
	public Set<String> CHROMOSOME=new HashSet<String>();
	@Option(shortName="BED", doc="limit to this bed file",optional=true)
	private File BEDFILE=null;
	
	
	private static class Cov
		{
		int sampleId=-1;
		int tid=-1;
		int pos=-1;
		int cov=0;
		}
	
	private static class CovCmp implements Comparator<Cov>
		{
		@Override
		public int compare(Cov o1, Cov o2)
			{
			int i=o1.tid-o2.tid;
			if(i!=0) return i;
			i=o1.pos-o2.pos;
			if(i!=0) return i;
			i=o1.sampleId-o2.sampleId;
			return i;
			}
		}

	
	private static class CovCodec extends AbstractDataCodec<Cov>
		{
		@Override
		public Cov decode(DataInputStream dis) throws IOException
			{
			Cov c=new Cov();
			try {
				c.sampleId=dis.readInt();
			} catch (Exception e) {
				return null;
				}
			c.tid=dis.readInt();
			c.pos=dis.readInt();
			c.cov=dis.readInt();
			
			return c;
			}
		
		@Override
		public void encode(DataOutputStream dos, Cov o) throws IOException
			{
			dos.writeInt(o.sampleId);
			dos.writeInt(o.tid);
			dos.writeInt(o.pos);
			dos.writeInt(o.cov);
			
			}
		@Override
		public CovCodec clone()
			{
			return new CovCodec();
			}
		}
	
	
	@Override
	public String getVersion() {
		return "1.0";
		}
	
	@Override
	protected int doWork()
		{
		List<SAMFileReader> samFileReaders=new ArrayList<SAMFileReader>(IN.size());
		PrintWriter out=new PrintWriter(System.out);
		SortingCollection<Cov> covCollections=null;
		CloseableIterator<Cov> covIter=null;
		try {

			SAMSequenceDictionary ssDict=null;
			for(File in:IN)
				{
				LOG.info("opening "+in);
				SAMFileReader sfr=new SAMFileReader(in);
				sfr.setValidationStringency(super.VALIDATION_STRINGENCY);
				samFileReaders.add(sfr);
				
				SAMSequenceDictionary dict=sfr.getFileHeader().getSequenceDictionary();
				
				if(dict==null || dict.isEmpty())
					{
					LOG.error("No dictionary in "+in);
					return -1;
					}
				if(ssDict==null)
					{
					ssDict=dict;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, ssDict))
					{
					LOG.error("Not the same sequence dictionaries "+IN.get(0)+" vs "+in);
					return -1;
					}
				}
			
			List<SamRecordFilter> samRecordFilters=new ArrayList<SamRecordFilter>();
			samRecordFilters.add(new SamRecordFilter()
				{
				@Override
				public boolean filterOut(SAMRecord first, SAMRecord second) {
						return filterOut(first) || filterOut(second);
					}
				@Override
				public boolean filterOut(SAMRecord record) {
					if(record.getReadFailsVendorQualityCheckFlag()) return true;
					if(record.getReadUnmappedFlag()) return true;
					if(record.getDuplicateReadFlag()) return true;
					if(record.getReadPairedFlag())
						{
						if(record.getMateUnmappedFlag()) return true;
						if(!record.getProperPairFlag()) return true;
						}
					return false;
					}
				});
			
			
			covCollections=SortingCollection.newInstance(
					Cov.class,
					new CovCodec(), 
					new CovCmp(),
					super.MAX_RECORDS_IN_RAM
					);
			covCollections.setDestructiveIteration(true);

			
			for(SAMSequenceRecord ssr:ssDict.getSequences())
				{
				if(!CHROMOSOME.isEmpty() && !CHROMOSOME.contains(ssr.getSequenceName()))
					{
					LOG.info("ignoring "+ssr.getSequenceName());
					continue;
					}
				
				int coverage[]=new int[1+ssr.getSequenceLength()];
				
				for(int t=0;t< samFileReaders.size();++t)
					{
					Arrays.fill(coverage, 0);
					SAMFileReader sfr=samFileReaders.get(t);
					
					SamLocusIterator sli=null;
					
					if(BEDFILE==null)
		                {
						sli=new SamLocusIterator(sfr);
		                }
	                else
	                	{
		    	     	IntervalList intervalList=new IntervalList(sfr.getFileHeader());
		    	     	BufferedReader in=new BufferedReader(new FileReader(BEDFILE));
	    	     	    String line=null;
	    	     	    while((line=in.readLine())!=null)
	    	     	    	{
	    	     	    	if(line.isEmpty() || line.startsWith("#")) continue;
	    	     	    	String tokens[]=line.split("[\t]");
	    	     	    	if(!CHROMOSOME.isEmpty() && !CHROMOSOME.contains(tokens[0]))
		    					{
		    					continue;
		    					}
	    	     	    	Interval interval=new Interval(tokens[0], 1+Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
	    	     	    	intervalList.add(interval);
	    	     	    	}
	    	     	    in.close();
		    	        intervalList.sort();	                    
	                    sli=new SamLocusIterator(sfr,intervalList);
	                	}
					
					
					sli.setSamFilters(samRecordFilters);
					sli.setEmitUncoveredLoci(true);
					
					
					while(sli.hasNext())
						{
						LocusInfo rec=sli.next();
						int pos=rec.getPosition();
						if(pos<1 || pos>=coverage.length) continue;
						coverage[pos]=rec.getRecordAndPositions().size();
						}
					sli.close();
					
					for(int i=1;i+this.WINDOW_SIZE <=coverage.length;i+=this.WINDOW_SHIFT)
						{
						Cov c=new Cov();
						c.sampleId=t;
						c.tid=ssr.getSequenceIndex();
						c.pos=i;
						int n=0;
						double depth=0;
						for(int j=0;j<this.WINDOW_SIZE && i+j< coverage.length;++j,++n)
							{
							depth+=coverage[i+j];
							}
						c.cov=(int)((depth)/n);
						covCollections.add(c);
						}
					}
				}
			covCollections.doneAdding();
			
			if(OUT!=null)
				{
				out=new PrintWriter(IOUtils.openFileForBufferedWriting(OUT));
				}
			if(HEADER)
				{
				out.println("NAME");
				for(File in:IN)
					{
					out.print('\t');
					out.print(in);
					}
				out.println();
				}
			
			covIter=covCollections.iterator();
			Map<Integer,Cov> curr=new HashMap<Integer,Cov>(this.IN.size());
			for(;;)
				{
				Cov c=null;
				if(covIter.hasNext())
					{
					c=covIter.next();
					}
				
				Cov first=(curr.isEmpty()?null:curr.values().iterator().next());
				
				if(c!=null && first!=null && c.pos==first.pos && c.tid==first.tid)
					{
					curr.put(c.sampleId,c);
					continue;
					}
				
				if(first!=null)
					{
					out.print(ssDict.getSequence(first.tid).getSequenceName()+"_"+first.pos+"_"+(first.pos+this.WINDOW_SIZE));
					for(int i=0;i< IN.size();++i)
						{
						out.print('\t');
						Cov c1=curr.get(i);
						out.print(c1==null?0:c1.cov);
						}
					curr.clear();
					out.println();
					}
					
				if(c==null) break;
				curr.put(c.sampleId,c);
				}
			
			

			
			return 0;
			} 
		catch (Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			if(covIter!=null) covIter.close();
			for(SAMFileReader sf:samFileReaders) sf.close();
			if(out!=null) out.flush();
			if(out!=null) out.close();
			}
		
		}
	
	public static void main(String[] args) {
		new Bam4DeseqIntervals().instanceMainWithExit(args);
		}
	}
