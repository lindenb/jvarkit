package com.github.lindenb.jvarkit.tools.bigwigmerge;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.wig.BigWigReader;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
/**
BEGIN_DOC

# motivation

merge several Bigwig files using different descriptive statistics (mean, median, etc..)

Output is a BedGraph file.

Input is a set of bigwig file or a file with the '.list' suffix containing the path to the bigwig


## Example

```
find DIR -type f -name "*.bigWig" > tmp.list
java -jar jvarkit.jar bigwigmerge -R genome.fa tmp.list --interval "chr1:234-567" --header --method median  > bedGraph.bed
```


END_DOC

**/
@Program(
		name="bigwigmerge",
		description="merge several Bigwig files using different descriptive statistics (mean, median, etc..)",
		keywords={"wig","bigwig"},
		creationDate="20240417",
		modificationDate="20240417",
		jvarkit_amalgamion =  true
		)
public class BigwigMerge extends Launcher {
	private static final Logger LOG = Logger.of(BigwigMerge.class);
	private enum Method {median,average,count,sum,min,max,random};
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE,required = true)
	private Path referenceFile = null;
	@Parameter(names={"--interval","-r"},description=" process only this interval")
	private String intervalStr= null;
	@Parameter(names={"--method","-m"},description=" how to compute merge of wiggle values ?")
	private Method method = Method.median;
	@Parameter(names={"--min-value"},description=" skip output data with values lower than 'x'")
	private Double treshold_value = null;
	@Parameter(names={"--min-item-count"},description=" skip output data if thre is less than 'x' bigwig file at a location")
	private int min_item_count=0;
	@Parameter(names={"--header"},description="write track header")
	private boolean with_track_header=false;

	private final Random random = new Random(System.currentTimeMillis());
	
	
	private static class OneWigBase {
		final String contig;
		final int pos1;
		final float value;
		OneWigBase(final String contig,int pos1,float value) {
			this.contig = contig;
			this.pos1=pos1;
			this.value=value;
			}
		@Override
		public String toString() {
			return contig+":"+pos1+":"+value;
			}
		}
	
	private static class WigItemImpl implements ExtendedLocatable {
		final String contig;
		final int start;
		int end;
		final float value;
		WigItemImpl(List<OneWigBase> L,float value) {
			if(L.isEmpty()) throw new IllegalArgumentException("empty list");
			this.contig = L.get(0).contig;
			this.start = L.get(0).pos1;
			this.end= this.start;
			this.value=value;
			if(L.stream().anyMatch(it->!it.contig.equals(this.contig) || it.pos1!=start)) throw new IllegalArgumentException("boum contig");
			}
		@Override public String getContig() {
			return this.contig;
			};
		public int getStart() {
			return start;
			}
		@Override
		 public int getEnd() {
			return end;
		 	}
		}


	
	private static class OneWigBaseIterator extends AbstractIterator<OneWigBase> {
		final BigWigReader.WigItem item;
		int pos;
		OneWigBaseIterator(final BigWigReader.WigItem item) {
			this.item=item;
			this.pos= item.getStart();
			}
		protected OneWigBase advance() {
			if(pos> this.item.getEnd()) return null;
			final OneWigBase b = new OneWigBase(item.getContig(),this.pos,this.item.getValue());
			pos++;
			return b;
			}
		}
	private static class OneWigBaseIterator2 extends AbstractIterator<OneWigBase> {
		private final PeekableIterator<BigWigReader.WigItem> delegate;
		private OneWigBaseIterator currentIter =null;
		OneWigBaseIterator2(final CloseableIterator<BigWigReader.WigItem> delegate) {
			this.delegate=new PeekableIterator<>(delegate);
			}
		protected OneWigBase advance() {
			for(;;)
				{
				if(currentIter==null) {
					if(!delegate.hasNext()) {
						return null;
						}
					final BigWigReader.WigItem wi = delegate.next();
					//paranoid
					if(delegate.hasNext()) {
						final BigWigReader.WigItem wj = delegate.peek();
						if(wi.overlaps(wj)) throw new IllegalStateException("overlapping items in the same are not supported "+wi+" "+wj);
						}
					currentIter = new OneWigBaseIterator(wi);
					}
				if(!currentIter.hasNext()) {
					currentIter=null;
					continue;
					}
				return  currentIter.next();
				}
			}
		}
	/** merging iterator, choose items with lowest position */
	private static class OneWigBaseIterator3 extends AbstractIterator<OneWigBase> {
		final  List<OneWigBaseIterator2> delegates;
		OneWigBaseIterator3(final List<OneWigBaseIterator2> delegates) {
			this.delegates= new ArrayList<>(delegates);
			}
		protected OneWigBase advance() {
			for(;;)
				{
				if(delegates.isEmpty()) return null;
				int best_idx =-1;
				OneWigBase best=null;
				int i=0;
				while(i< delegates.size()) {
					OneWigBaseIterator2 iter = delegates.get(i);
					if(!iter.hasNext()) {
						delegates.remove(i);
						continue;
						}
					else
						{
						final OneWigBase b = iter.peek();
						if(best==null || b.pos1< best.pos1) {
							best = b;
							best_idx=i;
							}
						i++;
						}
					}
				if(best_idx<0) return null;
				return delegates.get(best_idx).next();
				}
			}
		}

	private float computeValue(final List<OneWigBase> row) {
		switch(this.method) {
			case random: {
				final double[] array = row.stream().mapToDouble(it->(double)it.value).toArray();
				return (float)array[random.nextInt(array.length)];
				}
			case count: return row.size();
			case median:{
				return (float)new Median().evaluate(
						row.stream().mapToDouble(it->(double)it.value).toArray()
						);
				}
			case average:
				return (float)row.stream().mapToDouble(it->(double)it.value).average().orElse(0.0);
			case sum:
				return (float)row.stream().mapToDouble(it->(double)it.value).sum();
			case min:
				return (float)row.stream().mapToDouble(it->(double)it.value).min().getAsDouble();
			case max:
				return (float)row.stream().mapToDouble(it->(double)it.value).max().getAsDouble();
			default: throw new IllegalStateException();
			}
		}
	
	private void mergeWiggle(final PrintWriter pw ,final Locatable loc,	final List<BigWigReader> readers) {
		final List<CloseableIterator<BigWigReader.WigItem>> iterators1 = new ArrayList<>(readers.size());
		for(BigWigReader bwr: readers) {
			iterators1.add(bwr.query(loc));
			}
		final List<OneWigBaseIterator2> iterators2 = iterators1.stream().map(it->new OneWigBaseIterator2(it)).collect(Collectors.toList());
		final OneWigBaseIterator3 iterator3 = new OneWigBaseIterator3(iterators2); 
		try(final EqualRangeIterator<OneWigBase> iter = new EqualRangeIterator<>(iterator3, (A,B)->Integer.compare(A.pos1, B.pos1))) {
			WigItemImpl prev=null;
			for(;;) {
				WigItemImpl curr=null;
				if(iter.hasNext()) {
					final List<OneWigBase> row = iter.next();
					if(row.size() < this.min_item_count) {
						continue;
						}
					curr = new WigItemImpl(row, computeValue(row));
					}
				if(prev!=null && curr!=null && prev.getEnd()+1==curr.getStart() && prev.value==curr.value) {
					prev.end=curr.getStart();
					continue;
					}
				
				if(curr==null || prev!=null)
					{
					if(prev!=null) {
						if(treshold_value==null || prev.value >= treshold_value.doubleValue()) {
							pw.print(prev.toBed3());
							pw.print("\t");
							pw.print(prev.value);
							pw.println();
							}
						prev=null;
						}
					if(curr==null) break;
					}
				prev=curr;
				}
			}
		iterators1.stream().forEach(IT->IT.close());
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<BigWigReader> readers = new ArrayList<>();
		try {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(referenceFile);
			final Locatable loc = intervalStr==null?null:
				new IntervalParser(dict).apply(this.intervalStr).
				orElseThrow(IntervalParser.exception("Cannot parse "+intervalStr));

			
			for(String wiggleURI: IOUtils.unrollStrings(args)) {
				readers.add(new BigWigReader(wiggleURI));
				}
			if(readers.isEmpty()) {
				LOG.error("input missing");
				return -1;
				}
			if(readers.size()< this.min_item_count) {
				LOG.error("min-item-count < coun(bigwig)");
				return -1;
				}
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				if(with_track_header) {
					pw.println("track type=bedGraph name=\"__TITLE__\" description=\""+String.join(" ", args)+" merged using jvarkit.bigwigmerge method:"+this.method+"\"");
					}
				if(loc==null) {
					for(final SAMSequenceRecord ssr: dict.getSequences()) {
						mergeWiggle(pw,ssr,readers);
						}
					}
				else
					{
					mergeWiggle(pw,loc,readers);
					}
				pw.flush();
				}
			
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			readers.stream().forEach(R->R.close());
		}
	}	
	
public static void main(String[] args) {
	new BigwigMerge().instanceMainWithExit(args);
	}
}
