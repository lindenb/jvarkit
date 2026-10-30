/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.log;


import java.util.Map;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;


/** progress log reloaded... */
public class ProgressDict<T extends Locatable> extends AbstractProgress implements UnaryOperator<T> {
	private static Logger LOG=Logger.of(ProgressDict.class);
	
	protected long genome_length = 0L;
	private final Map<String,Long> contig2pos;

	public ProgressDict(final Logger log,final SAMSequenceDictionary dict) {
		super(LOG);
		this.contig2pos = SequenceDictionaryUtils.makeContigToIndexMap(dict);
		this.genome_length = dict.getReferenceLength();
		}
	public ProgressDict(final SAMSequenceDictionary dict) {
		this(LOG,dict);
		}
	
	/** convert current item to label */
	public String toString(final T loc) {
		if(loc==null) return "null";
		if(loc instanceof VariantContext) {
			final VariantContext vc=VariantContext.class.cast(loc);
			return new StringBuilder()
					.append(vc.getContig())
					.append(":")
					.append(StringUtils.niceInt(vc.getStart()))
					.append(":")
					.append(vc.getReference().getBaseString())
					.toString()
					;
				}
		else if(loc instanceof SAMRecord) {
			final SAMRecord rr=SAMRecord.class.cast(loc);
			return new StringBuilder()
					.append(rr.getContig())
					.append(":")
					.append(StringUtils.niceInt(rr.getStart()))
					.append("-")
					.append(StringUtils.niceInt(rr.getEnd()))
					.append(":")
					.append(rr.getReadName())
					.toString()
					;
				}
		
		
		return new StringBuilder()
			.append(loc.getContig())
			.append(":")
			.append(StringUtils.niceInt(loc.getStart()))
			.append("-")
			.append(StringUtils.niceInt(loc.getEnd()))
			.toString()
			;
		}
	
	/** end item, with no message */
	public T apply(final T loc) {
		return apply(loc,null);
		}
	
	/** end item, with message */
	public T apply(final T loc,final Object optional_msg) {
		int every_x_seconds = 60;
		final long now = System.currentTimeMillis();
		final long diff_millisec = now-super.last_messsage_millisec;
		boolean should_display=	isEnabled() && 
				loc!=null && 
				loc.getStart()>0 &&
				this.contig2pos.containsKey(loc.getContig()) && 
				( super.last_messsage_millisec<0L || (diff_millisec / 1_000 <= every_x_seconds));
		if(should_display) {
			final long index_in_genome= this.contig2pos.get(loc.getContig()).longValue()+loc.getStart();
			final long millisec_per_item  = (now - super.startMilliSec)/(index_in_genome+1L);
			final StringBuilder sb=new StringBuilder();
			sb.append(toString(loc))
				.append(". That took: ")
				.append(StringUtils.niceDuration(System.currentTimeMillis() - super.startMilliSec))
				.append(". Remain: ")
				.append(StringUtils.niceDuration((this.genome_length-(index_in_genome))*millisec_per_item))
				;
			if(optional_msg!=null) {
				sb.append(". ").append(optional_msg);
				}
			final String msg= sb.toString();
			getLogger().info(msg);
			super.last_messsage_millisec = now;
			}
		return loc;
		}

}
