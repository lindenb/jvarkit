/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Log.LogLevel;
import htsjdk.samtools.util.PeekableIterator;


public class CompareBams4  extends AbstractCompareBams4
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(CompareBams4.class);
	private LiftOver liftOver=null;
	private enum OnlyIn { ONLY_IN_FIRST,ONLY_IN_SECOND,BOTH};
	private enum LiftOverStatus {NoLiftOver,SourceUnmapped,DestNotInDict,LiftOverFailed,DiscordantChrom,SameChrom};
	private enum CompareContig {BothUnmapped,GainMapping,LostMapping,SameContig,DiscordantContig};
	private enum Shift {Zero,Le10,Le20,Le100,Gt100};
	private enum Flag { Discordant,Same};
	
	/* see https://github.com/samtools/hts-specs/issues/5 */
	private static int strnum_cmp(final String _a,final String _b)
	{

		char ca='\0',cb='\0';
	    int ia = 0, ib = 0;
	    while(ia< _a.length() && ib<_b.length()) {
	    	ca=(ia< _a.length()?_a.charAt(ia):'\0');
	    	cb=(ib< _b.length()?_b.charAt(ib):'\0');
	    	
	    	
	        if (Character.isDigit(ca) && Character.isDigit(cb))
	        	{
	            while (ca=='0'){
	            	++ia;
	            	ca=(ia< _a.length()?_a.charAt(ia):'\0');
	            	}
	            while (cb=='0'){
	            	++ib;
	            	cb=(ib< _b.length()?_b.charAt(ib):'\0');
	            	}
	            
	            
	            while (Character.isDigit(ca) && Character.isDigit(cb) && ca==cb)
	            	{
	            	++ia;
	            	++ib;
			    	ca=(ia< _a.length()?_a.charAt(ia):'\0');
			    	cb=(ib< _b.length()?_b.charAt(ib):'\0');
	            	}
	            
	            if (Character.isDigit(ca) && Character.isDigit(cb)) {
	                int i = 0;
	                while((ia+i)< _a.length() &&
	                	  (ib+i)< _b.length() &&
	                	  Character.isDigit(_a.charAt(ia+i)) &&
	                	  Character.isDigit(_b.charAt(ib+i))
	                	  ) {
	                	++i;
	                	}
			    	final char c1=((ia+i)< _a.length()?_a.charAt(ia+i):'\0');
			    	final char c2=((ib+i)< _b.length()?_b.charAt(ib+i):'\0');
	                return Character.isDigit(c1)? 1 : Character.isDigit(c2)? -1 : (int)ca - (int)cb;
	            }
	          else if (Character.isDigit(ca))
	        	  {
	        	  return 1;
	        	  }
	         else if (Character.isDigit(cb))
	        	 {
	        	 return -1;
	        	 }
	         else if (ia != ib)
	        	 {
	        	 return ia - ib;
	        	 }
	        	}/* end of is digit */
	        else
	        	{
	            if (ca != cb) {
	            	return (int)ca - (int)cb;
	            }
	            ++ia; ++ib;
	        	}
	    	}
	    if(ca==cb) return 0;
	    return ca=='\0'?1:cb=='\0'?-1: (int)ca - (int)cb;
	}

	
	private static int hash(final int prev,final Object o) {
		return prev*31 +(o==null?0:o.hashCode());
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private static int cmp(final Comparable a,final Comparable b) {
		if(a==null && b!=null) return -1;
		if(a!=null && b==null) return 1;
		if(a==null && b==null) return 0;
		return a.compareTo(b);
	}
	
	private static void str(final StringBuilder sb,final Object o) {
		sb.append(o==null?".":o.toString()).append("\t");
	}
	
	private static class Couple<T extends Comparable<T>> implements Comparable<Couple<T>>
	{
		final T v1;
		final T v2;
		public Couple(final T v1,final T v2) {
			this.v1 = v1;
			this.v2 = v2;
		}
		
		@Override
		public int compareTo(final Couple<T> o) {
			int i=v1.compareTo(o.v1);if(i!=0) return i;
			i=v2.compareTo(o.v2);if(i!=0) return i;
			return 0;
			}
		
		public int hashCode() { return v1.hashCode()*31+v2.hashCode();}
		@SuppressWarnings("unchecked")
		@Override
		public boolean equals(final Object obj) {
			if( obj==this) return true;
			if(obj == null || !(obj instanceof Couple)) return false;
			return this.compareTo(Couple.class.cast(obj)) == 0;
		}

		@Override
		public String toString() {
			return String.valueOf(v1)+"/"+v2;
			}
		
		
	}

	
	private static final class Diff
		implements Comparable<Diff>
		{
		OnlyIn onlyIn = null;
		LiftOverStatus liftover = null;
		CompareContig compareContig = null;
		Shift shift = null;
		Integer diffCigarOperations = null;
		Integer diffNM = null;
		Couple<Integer> diffFlags = null;
		Couple<String> diffChrom = null;
		Flag diffFlag;
		
		@Override
		public int hashCode() {
			int h = 0;
			h = hash(h,onlyIn);
			h = hash(h,liftover);
			h = hash(h,compareContig);
			h = hash(h,shift);
			h = hash(h,diffCigarOperations);
			h = hash(h,diffNM);
			h = hash(h,diffFlags);
			h = hash(h,diffChrom);
			h = hash(h,diffFlag);
			return h;
		}
		
		@Override
		public int compareTo(final Diff o) {
			if( o==this) return 0;
			int i;
			i = cmp(this.onlyIn,o.onlyIn); if(i!=0) return i;
			i = cmp(this.liftover,o.liftover); if(i!=0) return i;
			i = cmp(this.compareContig,o.compareContig); if(i!=0) return i;
			i = cmp(this.shift,o.shift); if(i!=0) return i;
			i = cmp(this.diffCigarOperations,o.diffCigarOperations); if(i!=0) return i;
			i = cmp(this.diffNM,o.diffNM); if(i!=0) return i;
			i = cmp(this.diffFlags,o.diffFlags); if(i!=0) return i;
			i = cmp(this.diffChrom,o.diffChrom); if(i!=0) return i;
			i = cmp(this.diffFlag,o.diffFlag); if(i!=0) return i;
			return 0;
			}
		
		@Override
		public boolean equals(final Object obj) {
			if( obj==this) return true;
			if(obj == null || !(obj instanceof Diff)) return false;
			return this.compareTo(Diff.class.cast(obj)) == 0;
			}
		
		@Override
		public String toString() {
			final StringBuilder sb = new StringBuilder();
			str(sb,this.onlyIn);
			str(sb,this.liftover);
			str(sb,this.compareContig);
			str(sb,this.shift);
			str(sb,this.diffCigarOperations);
			str(sb,this.diffNM);
			str(sb,this.diffFlags);
			str(sb,this.diffChrom);
			str(sb,this.diffFlag);
			return sb.toString();
			}
		}
	
	private static Interval interval(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return null;
		return new Interval(
			rec.getReferenceName(),
			rec.getAlignmentStart(),
			rec.getAlignmentEnd(),
			rec.getReadNegativeStrandFlag(),
			null
			);
	}
	
	private static int side(final SAMRecord rec) {
		if(rec.getReadPairedFlag()) {
			if(rec.getFirstOfPairFlag()) return 1;
			if(rec.getSecondOfPairFlag()) return 2;
			throw new RuntimeException("Side for flag "+rec.getReadName()+":"+rec.getFlags()+"?");
		}
		else {
			return 0;
		}
	}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		PrintWriter out=null;
		final List<String> args = super.getInputFiles();
		final SamReader samFileReaders[]=new SamReader[]{null,null};
		final SAMSequenceDictionary dicts[]=new SAMSequenceDictionary[]{null,null};
		@SuppressWarnings("unchecked")
		final PeekableIterator<SAMRecord> iters[]=new PeekableIterator[]{null,null};
		final List<List<SAMRecord>> recordLists = Arrays.asList(new ArrayList<SAMRecord>(),new ArrayList<SAMRecord>());
		final Counter<Diff> diffs = new Counter<>();
		final short NM_TAG = SAMTagUtil.getSingleton().NM;
		if(args.size() !=2)
			{
			return wrapException("Expected two and only two bams please, but got "+args.size());
			}
		
		try
			{
			
			
			final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			for(int i=0;i< args.size() && i< samFileReaders.length;++i)
				{
				final File samFile = new File(args.get(i));
				LOG.info("opening "+samFile);
				samFileReaders[i]=srf.open(samFile);
				final SAMFileHeader header = samFileReaders[i].getFileHeader();
				if(header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
					return wrapException("Expected "+samFile+" to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder());
					}
				dicts[i] = header.getSequenceDictionary();
				if( dicts[i]==null || dicts[i].isEmpty()) {
					return wrapException("In "+samFile+": No SAMSequenceDictionary in header.");
					}
				iters[i] = new PeekableIterator<>(samFileReaders[i].iterator());
				}
			
			if(super.chainFile!=null) {
				LOG.info("opening chain file "+super.chainFile+".");
				this.liftOver =new LiftOver(super.chainFile);
				this.liftOver.setLiftOverMinMatch(super.liftOverMismatch<=0?LiftOver.DEFAULT_LIFTOVER_MINMATCH:super.liftOverMismatch);

				if(!super.disableChainValidation) {
				this.liftOver.validateToSequences(dicts[1]);
				}
			}
			
			final Comparator<SAMRecord> comparator;
			
			if( super.samtoolsquerysort) {
				LOG.info("using the samtools sort -n comparator");
				comparator = new Comparator<SAMRecord>() {
					@Override
					public int compare(final SAMRecord o1, final SAMRecord o2) {
						int i= strnum_cmp(o1.getReadName(), o2.getReadName());
						if(i!=0) return i;
						return side(o1)-side(o2);
						}
					};			
				}
			else
				{
				comparator = new Comparator<SAMRecord>() {
					@Override
					public int compare(final SAMRecord o1, final SAMRecord o2) {
						int i=o1.getReadName().compareTo(o2.getReadName());
						if(i!=0) return i;
						return side(o1)-side(o2);
						}
					};
				}
			
			for(;;) {
				for(int i=0;i< 2;++i) {
					if(recordLists.get(i).isEmpty())
						{
						while(iters[i].hasNext()) {
							final SAMRecord rec = iters[i].peek();
							if(rec.isSecondaryOrSupplementary()) {
								iters[i].next();
								continue;
								}
							if(!recordLists.get(i).isEmpty() && comparator.compare(recordLists.get(i).get(0),rec)>0)
								{
								throw new SAMException("Something is wrong in sort order of "+args.get(i)+" : got\n\t"
										+rec+" "+side(rec)+"\nafter\n\t"+ recordLists.get(i).get(0)+" "+side(recordLists.get(i).get(0))+"\nSee also option -"+OPTION_SAMTOOLSQUERYSORT
										);
								}
							else if( recordLists.get(i).isEmpty() ||
								comparator.compare(recordLists.get(i).get(0),rec)==0)
								{
								recordLists.get(i).add(iters[i].next());
								}
							else
								{	
								break;
								}
							}
						}
				}
				if(recordLists.get(0).isEmpty() && recordLists.get(1).isEmpty() ) break;
				
				if(recordLists.get(0).size()>1) LOG.warn("size>2 for 1:"+recordLists.get(0));
				if(recordLists.get(1).size()>1) LOG.warn("size>2 for 2:"+recordLists.get(1));
				
				final Diff diff = new Diff();
				
								
				if((recordLists.get(0).isEmpty() && !recordLists.get(1).isEmpty()) ||
				   (!recordLists.get(0).isEmpty() && !recordLists.get(1).isEmpty() && comparator.compare(recordLists.get(0).get(0), recordLists.get(1).get(0))>0)
					) {
					diff.onlyIn = OnlyIn.ONLY_IN_SECOND;
					recordLists.get(1).clear();
					}
				else if((!recordLists.get(0).isEmpty() && recordLists.get(1).isEmpty()) ||
						(!recordLists.get(0).isEmpty() && !recordLists.get(1).isEmpty() && comparator.compare(recordLists.get(0).get(0), recordLists.get(1).get(0))<0))
					{
					diff.onlyIn = OnlyIn.ONLY_IN_FIRST;
					recordLists.get(0).clear();
					}
				else
				{
					final SAMRecord rec0 = recordLists.get(0).get(0);
					final SAMRecord rec1 = recordLists.get(1).get(0);
					final Interval rgn0a = interval(rec0);
					final Interval rgn0b;
					
					diff.onlyIn = OnlyIn.BOTH;
					diff.diffFlags = new Couple<Integer>(rec0.getFlags(),rec1.getFlags());
					diff.diffFlag = (rec0.getFlags()==rec1.getFlags()?Flag.Same:Flag.Discordant);

					
					if(this.liftOver==null) {
						rgn0b = rgn0a;
						diff.liftover = LiftOverStatus.NoLiftOver;
						}
					else if(rgn0a==null) {
						diff.liftover = LiftOverStatus.SourceUnmapped;
						rgn0b = rgn0a;
						}
					else {
						rgn0b = this.liftOver.liftOver(rgn0a);
						if(rgn0b==null) {
							diff.liftover = LiftOverStatus.LiftOverFailed;
							}
						else if(dicts[1].getSequence(rgn0b.getContig())==null)
							{
							diff.liftover = LiftOverStatus.DestNotInDict;
							}
						else if(rgn0a.getContig().equals(rgn0b.getContig())) {
							diff.liftover = LiftOverStatus.SameChrom;
							}
						else
							{
							diff.liftover = LiftOverStatus.DiscordantChrom;
							}
						}
					
					final Interval rgn1 = interval(rec1);
					if(rgn0b==null && rgn1==null) {
						diff.compareContig = CompareContig.BothUnmapped;
					} else if(rgn0b==null && rgn1!=null) {
						diff.compareContig = CompareContig.GainMapping;
						diff.diffChrom = new Couple<String>("*",rgn1.getContig());

					} else if(rgn0b!=null && rgn1==null) {
						diff.compareContig = CompareContig.LostMapping;
						diff.diffChrom = new Couple<String>(rgn0b.getContig(),"*");
					} else if(rgn0b.getContig().equals(rgn1.getContig()))
					{
						diff.compareContig = CompareContig.SameContig;
						diff.diffChrom = new Couple<String>(rgn0b.getContig(),rgn1.getContig());
						final int shift = Math.abs(rgn0b.getStart() - rgn1.getStart());
						if(shift == 0 ) {
							diff.shift = Shift.Zero;
						}
						else if(shift <=10) {
							diff.shift = Shift.Le10;
							
						}else if(shift <=20) {
							diff.shift = Shift.Le20;
						}
						else if( shift<= 100)
						{
							diff.shift = Shift.Le100;
						} else {
							diff.shift = Shift.Gt100;
						}
						
					} else
					{
						
						if(dicts[1].getSequence(rgn0b.getContig())==null)
							{
							diff.diffChrom = new Couple<String>("?",rgn1.getContig());
							}
						else
							{
							diff.diffChrom = new Couple<String>(rgn0b.getContig(),rgn1.getContig());
							}
						diff.compareContig = CompareContig.DiscordantContig;
					}
					
					if(rgn0b!=null && rgn1!=null) {
						diff.diffCigarOperations = (
								rec0.getCigar().numCigarElements() - rec1.getCigar().numCigarElements()
								);
						final Object nm0 = rec0.getAttribute(NM_TAG);
						final Object nm1 = rec1.getAttribute(NM_TAG);
						if(nm0!=null && nm1!=null) {
							diff.diffNM = (Number.class.cast(nm0).intValue()-Number.class.cast(nm1).intValue());
						}
					}
					recordLists.get(1).clear();
					recordLists.get(0).clear();
					}
				
				diffs.incr(diff);
			}
			
			
			LOG.info("done");
			
			final StringBuilder sb = new StringBuilder();
			str(sb,"onlyIn");
			str(sb,"liftover");
			str(sb,"compareContig");
			str(sb,"shift");
			str(sb,"diffCigarOperations");
			str(sb,"diffNM");
			str(sb,"diffFlags");
			str(sb,"diffChroms");
			str(sb,"diffFlag");
			sb.append("Count");
			out = super.openFileOrStdoutAsPrintWriter();
			out.println(sb);
			for(final Diff key:diffs.keySet()) {
				out.print(key.toString());
				out.println(diffs.count(key));
			}
			out.flush();
			out.close();
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			this.liftOver = null;
			for(int i=0;i< samFileReaders.length;++i){
				CloserUtil.close(iters[i]);
				CloserUtil.close(samFileReaders[i]);
				}
			CloserUtil.close(out);
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		 Log.setGlobalLogLevel(LogLevel.ERROR);
		new CompareBams4().instanceMainWithExit(args);
		}
}
