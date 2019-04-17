/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import java.io.Closeable;
import java.util.Iterator;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.readers.DelegateVcfIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class ProgressFactory {
private static final Logger LOG=Logger.build( ProgressFactory.class).make();
private SAMSequenceDictionary _dictionary = null;
private Logger _logger = LOG;
private int everySeconds = ProgressFactory.getDefaultEverySeconds();
private String _prefix;
private boolean _validateContigInDict = true;
private boolean _validateOrderInDict = false;
private boolean _runInBackground  = ProgressFactory.isDefaultRunningInBackground();
private boolean _silent  = ProgressFactory.isDefaultSilent();

private ProgressFactory() {
}


/** default is false, unless ` -Djvarkit.progress.silent=true` is defined */
public static boolean isDefaultSilent() {
	try {
		if("true".equals(System.getProperty("jvarkit.progress.silent", ""))) return true;
	} catch (Throwable e) {
		//ignore
	}
	return false;
	}

/* private stuff at cea */
private static boolean isRunningAtTgcc() {
	try {
		final String s = System.getenv("CCC_ORG");
		if(s!=null && s.startsWith("fg")) return true;
		}
	catch (Throwable e) {
		}
	return false;
	}

/** default is false, unless ` -Djvarkit.progress.background=true` is defined */
public static boolean isDefaultRunningInBackground() {
	try {
		if("true".equals(System.getProperty("jvarkit.progress.background", ""))) return true;
		if(isRunningAtTgcc()) return false;
		}
	catch (final Throwable e) {
		}
	return false;
	}

/** default is 10, unless ` -Djvarkit.progress.everysecs=1234` is defined */
public static int getDefaultEverySeconds() {
	try {
		final String s=System.getProperty("jvarkit.progress.everysecs", null);
		if(s!=null) {
			int n = Integer.parseInt(s);
			if(n>=1) return n;
			}
		}
	catch(final NumberFormatException err) {
		//ignore
		}
	if(isRunningAtTgcc()) return 60;
	return 10;
	}


/** run logger in background */
public ProgressFactory threaded(final boolean b) {
	this._runInBackground = b;
	return this;
	}

public boolean isThreaded() {
	return this._runInBackground;
	}

/** do not log anything */
public ProgressFactory setSilent(boolean _silent) {
	this._silent = _silent;
	return this;
	}

public boolean isSilent() {
	return _silent;
	}

/** validate contig: check the contig is in the dictionary, if the dictionary is defined */
public ProgressFactory validatingContig(final boolean b) {
	this._validateContigInDict = b;
	return this;
	}

public boolean isValidatingContig() {
	return this._validateContigInDict;
	}

/** validate contig: check items are sorted in the dictionary, it the dictionary is defined */
public ProgressFactory validatingSortOrder(final boolean b) {
	this._validateOrderInDict = b;
	return this;
	}

public boolean isValidatingSortOrder() {
	return this._validateOrderInDict;
	}

/** define log prefix */
public ProgressFactory prefix(final String prefix) {
	this._prefix = prefix;
	return this;
	}

public String getPrefix() {
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

public SAMSequenceDictionary getDictionary() {
	return this._dictionary;
	}

public ProgressFactory logger(final Logger logger) {
	this._logger = logger;
	return this;
	}

public Logger getLogger() {
	return (this._logger==null?LOG:this._logger);
	}

public ProgressFactory setEverySeconds(int everySeconds) {
	this.everySeconds = everySeconds;
	return this;
}

public int getEverySeconds() {
	return everySeconds;
	}

public static ProgressFactory newInstance()
	{
	return new ProgressFactory();
	}

public <T extends Locatable> Stream<T> stream(final CloseableIterator<T> delegate) {
	final CloseableIterator<T> iter= build(delegate);
	final Stream<T> st = iter.stream();
	st.onClose(()->iter.close());
	return st;
	}

public <T extends Locatable> CloseableIterator<T> build(final Iterator<T> delegate) {
	return new LogIter<>(delegate,build());
	}

public VCFIterator build(final VCFIterator delegate) {
	return new VCFIteratorWatcher(delegate,build());
	}

public <T extends Locatable> Watcher<T> build() {
	if(this._silent) return new SilentWatcher<>();
	
	
	final WatcherImpl<T> w = new WatcherImpl<>();
	w._dictionary = this.getDictionary();
	w._logger = this.getLogger();
	w._everySeconds= this.getEverySeconds();
	
	w._checkDictContig = this.getDictionary()!=null && this.isValidatingContig();
	w._checkSorted = this.isValidatingSortOrder();
	w._threaded = this.isThreaded();

	w._logPrefix=(StringUtil.isBlank( this.getPrefix())?"":"["+ this.getPrefix()+"]");
	return w;
	}

public static interface Watcher<T extends Locatable>
	extends Closeable,Function<T, T>
	{
	@Override
	public void close();
	}

private static class SilentWatcher<T extends Locatable>
	implements Watcher<T>
	{
	@Override
	public T apply(final T t) {
		return t;
		}
	@Override
	public void close() {
		
		}
	}

private static class WatcherImpl<T extends Locatable>
	implements Watcher<T>,Runnable
	{
	private boolean EOF_flag = false;
	private SAMSequenceDictionary _dictionary = null;
	private Logger _logger = LOG;
	private boolean firstCall = true;
	private ScheduledExecutorService scheduledExecutorService=null;
	private int _everySeconds = 10;
	private boolean _checkDictContig = false;
	private boolean _checkSorted = false;
	private boolean _threaded = true;
	private T previousLocatable = null;
	private long count_items = 0L;
	private long startMillisec = System.currentTimeMillis();
	private transient long lastCallMillisec = startMillisec;
	private boolean dataAreSorted=true;
	private long cumulLengthDone[]=null;
	private long referenceLength=0L;
	private String _logPrefix = null;

	
	@Override
	public void run() {
		if(this.EOF_flag) return;
		final long now = System.currentTimeMillis();
		final long diff_millisec = now - this.lastCallMillisec;
		if(this._threaded ) {
			if(this.scheduledExecutorService==null) return;
			}
		else
			{
			if(diff_millisec / 1_000 <= this._everySeconds) return;
			}
		final T last = this.previousLocatable;
		this.lastCallMillisec = now;
		final long count = this.count_items;
		final String pfx= this._logPrefix;
		if(last==null) {
			_logger.info(pfx+"No data received. Elapsed time: "+duration(now-this.startMillisec));
			return ;
			}
		final int tid = getTid(last);
		final int pos = tid<0?-1:last.getStart();
		if(this._dictionary==null || !this.dataAreSorted || tid==-1 || pos<1) {
			_logger.info(pfx+"Last "+loc2str(last)+" "+duration(diff_millisec));
			return ;
			}
			
		final long numBasesDone=(tid==0?0:this.cumulLengthDone[tid-1])+pos;
		final long numBasesRemains=Math.max(0,referenceLength-numBasesDone);
		
		final double percentDone=numBasesDone/(double)this.referenceLength;
		final double millisecPerBase=(double)(now-this.startMillisec)/numBasesDone;
		final long timeRemain=(long)(numBasesRemains*millisecPerBase);
		
		
		final String msg = String.format(
				"%sCount: %s Elapsed: %s(%.2f%%) Remains: %s(%.2f%%) Last: %s:%s",
				pfx,
				this.format(count),
				
				duration(now-this.startMillisec),
				(percentDone*100.0),
				
				duration(timeRemain),
				(100-percentDone*100.0),
				
				last.getContig(),
				this.format(pos)
				);
		this._logger.info(msg);
		}
	
		@Override
		public T apply(final T item) {
		if(this.EOF_flag) throw new IllegalStateException("Walker was closed");
		
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
		

		if(item==null) {
			close();
			throw new IllegalStateException("Walker invoked with 'null'");
			}
		
		this.count_items++;
		
		//e.g: SAMRecord not mapped
		if(StringUtil.isBlank(item.getContig()))
			{
			this.previousLocatable = item;
			return item;
			}
		final int tid2 = getTid(item);
		
		
		if(this.previousLocatable!=null) {
				if(this.dataAreSorted) {
					final int tid1 =  getTid(this.previousLocatable);
					if((tid1!=-1 && tid2!=-1) && (tid1 > tid2 || (tid1==tid2 && this.previousLocatable.getStart()>item.getStart()))) {
						if(this._checkSorted) {
							close();
							throw new JvarkitException.BadLocatableSortOrder(this.previousLocatable, item, this._dictionary);
							}
						this.dataAreSorted=false;
						}
				}
			}
		this.previousLocatable = item;
		if(!this._threaded) run();
		return item;
		}
	@Override
	public void close() {
		if(EOF_flag) return;
		
		if(this.scheduledExecutorService!=null)
			{
			this.scheduledExecutorService.shutdown();
			this.scheduledExecutorService=null;
			}
		this.previousLocatable=null;
		this._logger.info(this._logPrefix +". Completed. N="+format(count_items)+". That took:"+duration(System.currentTimeMillis()-this.startMillisec));
		this.EOF_flag=true;
		}
	
	private String format(final long loc) {
		return StringUtils.niceInt(loc);
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
	
	private int getTid(final Locatable loc) {
		if(this._dictionary==null || loc==null) return -1;
		if(StringUtil.isBlank(loc.getContig())) return -1;//can happen for SamRecord
		final int tid = this._dictionary.getSequenceIndex(loc.getContig());
		if(tid<0 && this._checkDictContig)
			{
			throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), this._dictionary);
			}
		return tid;
		}
	}

private static class LogIter<T extends Locatable>
	extends AbstractIterator<T>
	implements CloseableIterator<T>
	{
	final Iterator<T> delegate;
	final Watcher<T> watcher;
	LogIter(final Iterator<T> delegate,final Watcher<T> watcher) {
		this.delegate = delegate;
		this.watcher = watcher;
		}
	@Override
	protected T advance() {
		if(!delegate.hasNext()) {
			close();
			return null;
			}
		return this.watcher.apply(delegate.next());
		}
	@Override
	public void close() {
		CloserUtil.close(delegate);
		CloserUtil.close(watcher);
		}
	}

private static class VCFIteratorWatcher extends DelegateVcfIterator
	{
	final Watcher<VariantContext> watcher;
	VCFIteratorWatcher(final VCFIterator delegate,final Watcher<VariantContext> watcher) {
		super(delegate);
		this.watcher=watcher;
		}
	@Override
	public boolean hasNext() {
		boolean b = super.hasNext();
		if(!b) close();
		return b;
		}
	@Override
	public VariantContext next() {
		return watcher.apply(super.next());
		}
	@Override
	public void close() {
		super.close();
		watcher.close();
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
