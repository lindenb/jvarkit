package com.github.lindenb.jvarkit.tools.impactdup;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.RuntimeEOFException;
import htsjdk.samtools.util.SamRecordIntervalIteratorFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;

@Program(name="impactofduplicates",
description="Impact of Duplicates per BAM.",
keywords={"bam"}
)
public class ImpactOfDuplicates extends Launcher
    {
	private static final Logger LOG = Logger.build(ImpactOfDuplicates.class).make();

    
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-B", description="BED File")
	private File BEDFILE = null;

	@ParametersDelegate
	private WritingSortingCollection sortingCollectionArgs=new WritingSortingCollection();
	
	
    /** sam file dict, to retrieve the sequences names  */
    private List< SAMSequenceDictionary> samFileDicts=new ArrayList<SAMSequenceDictionary>();
    /** buffer for Duplicates */
    private List<Duplicate> duplicatesBuffer=new ArrayList<Duplicate>();
    /** output */
    private PrintStream out=System.out;
    
    
    /* current index in BAM list */
    private int bamIndex;
    /* all duplicates, sorted */
    private SortingCollection<Duplicate> duplicates;
    
    private class Duplicate implements Comparable<Duplicate>
        {
        int tid;
        int pos;
        int size;
        int bamIndex;
        
        public String getReferenceName()
        	{
        	return samFileDicts.get(this.bamIndex).getSequence(this.tid).getSequenceName();
        	}
        public int compareChromPosSize(final Duplicate o)
        	{
        	int i=getReferenceName().compareTo(o.getReferenceName());
        	if(i!=0) return i;
        	i=pos-o.pos;
        	if(i!=0) return i;
        	i=size-o.size;
        	return i;
        	}
        
        
        
        @Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + bamIndex;
			result = prime * result + pos;
			result = prime * result + size;
			result = prime * result + tid;
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Duplicate other = (Duplicate) obj;
			return compareTo(other)==0;
			}
		@Override
        public int compareTo(final Duplicate o)
        	{
        	int i=compareChromPosSize(o);
        	if(i!=0) return i;
        	return bamIndex-o.bamIndex;
        	}
		
		@Override
		public String toString() {
			return "(bamIndex:"+bamIndex+" pos:"+pos+" size:"+size+" tid:"+tid+")";
			}
		
        }
    
    

    private  class DuplicateCodec
    	extends BinaryCodec
        implements SortingCollection.Codec<Duplicate >
        {
    	
        @Override
        public void encode(final Duplicate d)
            {
        	
	        	this.writeInt(d.tid);
	        	this.writeInt(d.pos);
	        	this.writeInt(d.size);
	        	this.writeInt(d.bamIndex);
	        	
            }
        
    
        
        @Override
        public Duplicate decode()
            {
        	Duplicate d=new Duplicate();
        	try
	        	{
	            d.tid=this.readInt();
	        	}
        	catch(RuntimeEOFException err)
        		{
        		return null;
        		}
            d.pos=this.readInt();
            d.size=this.readInt();
            d.bamIndex=this.readInt();
	        return d;
            }
        
        @Override
        public DuplicateCodec clone()
            {
            return new DuplicateCodec();
            }
        }
    
   
    
    private void dumpDuplicatesBuffer(final List<File> INPUT)
		{
    	if(this.duplicatesBuffer.isEmpty()) return;
    	int counts[]=new int[INPUT.size()];
    	Arrays.fill(counts, 0);
    	int maxDup=0;
    	for(int i=0;i< this.duplicatesBuffer.size();++i)
    		{
    		Duplicate di=this.duplicatesBuffer.get(i);
    		counts[di.bamIndex]++;
    		maxDup=Math.max(maxDup,counts[di.bamIndex]);
    		}
    	
    	
    	
    	if(maxDup<10)
    		{
    		this.duplicatesBuffer.clear();
    		return;
    		}
    	
    	int total=0;
    	for(int i:counts) total+=i; 
    	
    	Duplicate front=this.duplicatesBuffer.get(0);
    	out.print(
    			front.getReferenceName()+":"+
    			front.pos+"-"+
    			(front.pos+front.size)
    			);
    	out.print("\t"+maxDup+"\t"+(int)(total/(1.0*counts.length)));
    	for(int i=0;i< counts.length;++i)
    		{
    		out.print('\t');
    		out.print(counts[i]);
    		}

    	out.println();
    	this.duplicatesBuffer.clear();
		}
    
    
    @Override
    public int doWork(final List<String> args) {
    
       
       CloseableIterator<Duplicate> dupIter=null;
       final List<File> INPUT = args.stream().map(S->new File(S)).collect(Collectors.toList());
       
        try
            {
        this.duplicates=SortingCollection.newInstance(
                    Duplicate.class,
                    new DuplicateCodec(),
                    new Comparator<Duplicate>()
    	            	{
    	        		@Override
    	        		public int compare(Duplicate o1, Duplicate o2)
    	        			{
    	        			return o1.compareTo(o2);
    	        			}
    	            	},
                    this.sortingCollectionArgs.getMaxRecordsInRam(),
                    this.sortingCollectionArgs.getTmpPaths()
                    );
        	
        	
            for(this.bamIndex=0;
        		this.bamIndex<  INPUT.size();
        		this.bamIndex++)
                {
            	int prev_tid=-1;
            	int prev_pos=-1;
            	long nLines=0L;
            	File inFile= INPUT.get(this.bamIndex);
            	LOG.info("Processing "+inFile);
                IOUtil.assertFileIsReadable(inFile);
                SamReader samReader=null;
                CloseableIterator<SAMRecord> iter=null;
                try
	                {
	                samReader= SamReaderFactory.make().validationStringency(ValidationStringency.LENIENT).open(inFile);
	                final SAMFileHeader header=samReader.getFileHeader();
	                this.samFileDicts.add(header.getSequenceDictionary());
	                if(BEDFILE==null)
		                {
		                iter=samReader.iterator();
		                }
	                else
	                	{
		    	     	   IntervalList intervalList=new IntervalList(header);

		    	     	    BufferedReader in=new BufferedReader(new FileReader(BEDFILE));
		    	     	    String line=null;
		    	     	    while((line=in.readLine())!=null)
		    	     	    	{
		    	     	    	if(line.isEmpty() || line.startsWith("#")) continue;
		    	     	    	String tokens[]=line.split("[\t]");
		    	     	    	Interval interval=new Interval(tokens[0], 1+Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
		    	     	    	intervalList.add(interval);
		    	     	    	}
		    	     	    in.close();
		    	        intervalList=intervalList.sorted();
	                    List<Interval> uniqueIntervals=IntervalList.getUniqueIntervals(intervalList,false);

	             	   SamRecordIntervalIteratorFactory sriif=new  SamRecordIntervalIteratorFactory();
	            	    iter=sriif.makeSamRecordIntervalIterator(samReader, uniqueIntervals, false);
	                	}
	                while(iter.hasNext())
	                    {
	                    SAMRecord rec=iter.next();
	                    if(rec.getReadUnmappedFlag()) continue;
	                    if(!rec.getReadPairedFlag()) continue;
	                    if(rec.getReferenceIndex()!=rec.getMateReferenceIndex()) continue;
	                    if(!rec.getProperPairFlag()) continue;
	                    if(!rec.getFirstOfPairFlag()) continue;
	                    
	                    if(prev_tid!=-1 )
	                    	{
	                    	if(prev_tid> rec.getReferenceIndex())
	                    		{
	                    		throw new IOException("Bad sort order from "+rec);
	                    		}
	                    	else if(prev_tid==rec.getReferenceIndex() && prev_pos>rec.getAlignmentStart())
	                    		{
	                    		throw new IOException("Bad sort order from "+rec);
	                    		}
	                    	else
	                    		{
	                    		prev_pos=rec.getAlignmentStart();
	                    		}
	                    	}
	                    else
	                    	{
	                    	prev_tid=rec.getReferenceIndex();
	                    	prev_pos=-1;
	                    	}
	                    
	                    
	                    if((++nLines)%1000000==0)
	                    	{
	                    	LOG.info("In "+inFile+" N="+nLines);
	                    	}
	                    Duplicate dup=new Duplicate();
	                	dup.bamIndex=this.bamIndex;
	                	dup.pos=Math.min(rec.getAlignmentStart(),rec.getMateAlignmentStart());
	                	dup.tid=rec.getReferenceIndex();
	                	dup.size=Math.abs(rec.getInferredInsertSize());
	                	this.duplicates.add(dup);
	                    }
	                }
                finally
	                {
	                if(iter!=null) iter.close();
	                if(samReader!=null) samReader.close();
	                }
                LOG.info("done "+inFile);
                }
            /** loop done, now scan the duplicates */
            
            LOG.info("doneAdding");
            this.duplicates.doneAdding();
            
           
            this.out= super.openFileOrStdoutAsPrintStream(outputFile);
            	
            
        	out.print("#INTERVAL\tMAX\tMEAN");
        	for(int i=0;i< INPUT.size();++i)
        		{
        		out.print('\t');
        		out.print(INPUT.get(i));
        		}
        	out.println();

           dupIter=this.duplicates.iterator();
           while(dupIter.hasNext())
            	{
        	    Duplicate dup=dupIter.next();
	        	if( this.duplicatesBuffer.isEmpty() ||
	        		dup.compareChromPosSize(this.duplicatesBuffer.get(0))==0)
	            	{
	        		this.duplicatesBuffer.add(dup);
	            	}
	            else
	            	{
	            	dumpDuplicatesBuffer(INPUT);
	            	this.duplicatesBuffer.add(dup);
	            	}
            	}
            dumpDuplicatesBuffer(INPUT);
            LOG.info("end iterator");
            out.flush();
            out.close();
            }
        catch (Exception e) {
        	LOG.error(e);
            return -1;
            }
        finally
        	{
        	if(dupIter!=null) dupIter.close();
        	LOG.info("cleaning duplicates");
        	this.duplicates.cleanup();
        	}
        return 0;
        }
    public static void main(final String[] argv)
		{
	    new ImpactOfDuplicates().instanceMainWithExit(argv);
		}	

    }
