/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.log;

import java.text.DecimalFormat;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public abstract class ProgressFactory {
private static final Logger LOG=Logger.build( ProgressFactory.class).make();
private SAMSequenceDictionary _dictionary = null;
private Logger _logger = LOG;
private int everySeconds = 10;
private String _prefix;
private boolean _validateContigInDict = true;
private boolean _validateOrderInDict = false;
private boolean _runInBackground  = true;

protected ProgressFactory() {
}

/** run logger in background */
public ProgressFactory threaded(final boolean b) {
	this._runInBackground = b;
	return this;
	}

public boolean threaded() {
	return this._runInBackground;
	}


/** validate contig: check the contig is in the dictionary, if the dictionary is defined */
public ProgressFactory validateContig(final boolean b) {
	this._validateContigInDict = b;
	return this;
	}

public boolean validateContig() {
	return this._validateContigInDict;
	}

/** validate contig: check items are sorted in the dictionary, it the dictionary is defined */
public ProgressFactory validateSortOrder(final boolean b) {
	this._validateOrderInDict = b;
	return this;
	}

public boolean validateSortOrder() {
	return this._validateOrderInDict;
	}

/** define log prefix */
public ProgressFactory prefix(final String prefix) {
	this._prefix = prefix;
	return this;
	}

public String prefix() {
	return this._prefix;
	}



public ProgressFactory dictionary(final VCFHeader header) {
	return dictionary(header==null?null:header.getSequenceDictionary());
	}

public ProgressFactory dictionary(final SAMFileHeader header) {
	return dictionary(header==null?null:header.getSequenceDictionary());
	}

public ProgressFactory dictionary(final SamReader r) {
	return dictionary(r==null?null:r.getFileHeader());
	}

public ProgressFactory dictionary(final VCFFileReader r) {
	return dictionary(r==null?null:r.getFileHeader());
	}

public ProgressFactory dictionary(final SAMSequenceDictionary dictionary) {
	this._dictionary = dictionary;
	return this;
	}

public SAMSequenceDictionary dictionary() {
	return this._dictionary;
	}

public ProgressFactory logger(final Logger logger) {
	this._logger = logger;
	return this;
	}

public Logger logger() {
	return (this._logger==null?LOG:this._logger);
	}

public ProgressFactory setEverySeconds(int everySeconds) {
	this.everySeconds = everySeconds;
	return this;
}

public int everySeconds() {
	return everySeconds;
	}

public static ProgressFactory newInstance()
	{
	return new ProgressFactory()
		{
		};
	}
public <T extends Locatable> Stream<T> stream(final CloseableIterator<T> delegate) {
	final CloseableIterator<T> iter= build(delegate);
	final Stream<T> st = iter.stream();
	st.onClose(()->iter.close());
	return st;
	}

public <T extends Locatable> CloseableIterator<T> build(final CloseableIterator<T> delegate) {
	final AbstractLogIter<T> iter=new AbstractLogIter<>();
	iter._dictionary = this.dictionary();
	iter._logger = this.logger();
	iter._delegate = delegate;
	iter._everySeconds= this.everySeconds();
	iter._logPrefix= this.prefix();
	iter._checkDictContig = this.dictionary()!=null && this.validateContig();
	iter._checkSorted = this.validateSortOrder();
	iter._threaded = this.threaded();
	//
	iter._logPrefix=(StringUtil.isBlank( this.prefix())?"":"["+ this.prefix()+"]");
	return iter;
	}

private class AbstractLogIter<T extends Locatable>
	extends AbstractIterator<T>
	implements CloseableIterator<T>,Runnable
	{

	private SAMSequenceDictionary _dictionary = null;
	private Logger _logger = LOG;
	private CloseableIterator<T> _delegate = null;
	private boolean firstCall = true;
	private ScheduledExecutorService scheduledExecutorService=null;
	private int _everySeconds = 10;
	private boolean _checkDictContig = false;
	private boolean _checkSorted = false;
	private boolean _threaded = true;
	private T previous=null;
	private long count_items = 0L;
	private long startMillisec = System.currentTimeMillis();
	private transient long lastCallMillisec = startMillisec;
	private boolean dataAreSorted=true;
	private long cumulLengthDone[]=null;
	private long referenceLength=0L;
	private String _logPrefix = null;
	private final DecimalFormat niceInt = new DecimalFormat("###,###");

	
	@Override
	public void run() {
		final long now = System.currentTimeMillis();
		final long diff_millisec = now - this.lastCallMillisec;
		if(this._threaded ) {
			if(this.scheduledExecutorService==null) return;
			}
		else
			{
			if(diff_millisec* 1000 <this._everySeconds) return;
			}
		final T last = this.previous;
		this.lastCallMillisec = now;
		final long count = this.count_items;
		final String pfx= this._logPrefix;
		if(last==null) {
			_logger.info(pfx+"No data received. "+duration(diff_millisec));
			return ;
			}
		if(this._dictionary==null || !this.dataAreSorted) {
			_logger.info(pfx+"Last "+loc2str(last)+" "+duration(diff_millisec));
			return ;
			}
		
		final int tid = _dictionary.getSequenceIndex(last.getContig());
		final int pos = last.getStart();
		if(tid==-1 || pos<1)
			{
			_logger.info(last.toString()+"."+duration(diff_millisec));
			return;
			}
			
		final long numBasesDone=(tid==0?0:this.cumulLengthDone[tid-1])+pos;
		final long numBasesRemains=Math.max(0,referenceLength-numBasesDone);
		
		final double percentDone=numBasesDone/(double)this.referenceLength;
		final double millisecPerBase=(double)(diff_millisec)/numBasesDone;
		final long timeRemain=(long)(numBasesRemains*millisecPerBase);
		
		
		this._logger.info(
				String.format("%sCount: %d Elapsed: %s(%.2f%%) Remains: %s(%.2f%%) Last: %s:%s",
				
				pfx,
				format(count),
				
				duration(diff_millisec),
				(percentDone*100.0),
				
				duration(timeRemain),
				(100-percentDone*100.0),
				
				last.getContig(),
				this.format(pos)
				));
		}
	
	@Override
	protected T advance() {
		if(this.firstCall)
			{
			this.firstCall=false;
			this.startMillisec = System.currentTimeMillis();
			this.lastCallMillisec = this.startMillisec;
			if(this._dictionary==null) {
				this._checkDictContig = false;
				this.cumulLengthDone=new long[0];
				this.referenceLength=0L;
				}
			else
				{	
				this.cumulLengthDone=new long[this._dictionary.size()];
				long prev_cumul=0L;
				this.referenceLength=this._dictionary.getReferenceLength();
				for(int i=0;i< this._dictionary.size();++i)
					{
					cumulLengthDone[i]=prev_cumul;
					final SAMSequenceRecord ssr=this._dictionary.getSequence(i);
					prev_cumul += ssr.getSequenceLength();
					}
				}
			if(this._everySeconds>0 && this._threaded) {
				this.scheduledExecutorService = Executors.newSingleThreadScheduledExecutor();
				this.scheduledExecutorService.scheduleAtFixedRate(this,
					this._everySeconds,
					this._everySeconds,
					TimeUnit.SECONDS
					);
				}
			}
		if(this._delegate==null || !this._delegate.hasNext()) {
			close();
			return null;
		}
		final T item = this._delegate.next();
		if(item==null) {
			close();
			return null;
			}
		
		this.count_items++;
		
		if((this._checkDictContig || this._checkSorted) && this._dictionary.getSequence(item.getContig())==null) {
			close();
			throw new JvarkitException.ContigNotFoundInDictionary(item.getContig(), this._dictionary);
			}
		
		if(this.previous!=null) {
			
				if(this.dataAreSorted) {
				final int tid1 =  this._dictionary.getSequenceIndex(this.previous.getContig());
				final int tid2 =  this._dictionary.getSequenceIndex(item.getContig());
				if((tid1!=-1 && tid2!=-1) && (tid1 > tid2 || (tid1==tid2 && this.previous.getStart()>item.getStart()))) {
					if(this._checkSorted) {
						close();
						throw new JvarkitException.BadLocatableSortOrder(this.previous, item, this._dictionary);
						}
					this.dataAreSorted=false;
					}
				}
			
			this.previous=item;
			if(!this._threaded) run();
			}
		
		return item;
		}
	@Override
	public void close() {
		if(this._delegate!=null) {
			_delegate.close();
			_delegate=null;
			}
		if(this.scheduledExecutorService!=null)
			{
			this.scheduledExecutorService.shutdown();
			this.scheduledExecutorService=null;
			}
		this.previous=null;
		this._logger.info(this._logPrefix +". Completed. N="+format(count_items)+". That took:"+duration(System.currentTimeMillis()-this.startMillisec));
		}
	
	private String format(final long loc) {
		return this.niceInt.format(loc);
	}
	
	private  String loc2str(final Locatable loc) {
		if(loc==null) return "(null)";
		if(loc instanceof VariantContext) {
			final VariantContext ctx = VariantContext.class.cast(loc);
			return ctx.getContig()+":"+format(ctx.getStart())+":"+ctx.getReference().getDisplayString()+":"+ctx.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"));
			}
		if(loc instanceof SAMRecord) {
			final SAMRecord rec = SAMRecord.class.cast(loc);
			return rec.getReadName()+";"+rec.getContig()+";"+format(rec.getStart());
			}
		return loc.getContig()+":"+format(loc.getStart())+"-"+format(loc.getEnd());
		}

	}

private static String duration(long millisecs)
	{
	long n =millisecs/1000;
	if(n<60)
		{
		return n+" second"+(n<2?"":"s");
		}
	n/=60;//minutes
	if(n< 60)
		{
		return n+" minute"+(n<2?"":"s");
		}
	n/=60;//hours
	if(n< 24)
		{
		return n+" hour"+(n<2?"":"s");
		}
	n/=24;
	
	if(n< 365)
		{
		return n+" day"+(n<2?"":"s");
		}
	n/=365;
	
	return n+" year"+(n<2?"":"s");
	}


}
