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


import java.util.function.Consumer;

import com.github.lindenb.jvarkit.lang.StringUtils;


/** progress log reloaded... */
public class Progress extends AbstractProgress {
	private static Logger LOG=Logger.of(Progress.class);
	
	private long n_expected = 0L;
	private long current_index = 0L;

	public Progress(final Logger log,final long n_expected) {
		this(MSG->log.info(MSG),n_expected);
		}
	
	public Progress(final Consumer<String> log,final long n_expected) {
		super(log);
		this.n_expected=n_expected;
		}
	
	public Progress(final long n_expected) {
		this(LOG,n_expected);
		}
	/** end item, with no message */
	public void mark() {
		mark(this.current_index);
		this.current_index++;
		}
	
	/** end item, with message */
	public void mark(final long index) {
		if(!isEnabled()) return;
		final long now = System.currentTimeMillis();
		final int every_x_seconds = super.getEveryXSeconds();
		boolean should_display=	
				this.n_expected >0 && 
				index >= 0 &&
				index< this.n_expected &&
				(now-super.last_messsage_millisec) / 1_000 >= every_x_seconds;
				
		if(should_display) {
			final StringBuilder sb=new StringBuilder();
			sb.append(StringUtils.niceInt(this.current_index+1))
				.append(" / ")
				.append(StringUtils.niceInt(this.n_expected))
				.append(".");
			if(index>0) {
				final double millisec_per_item  = (now - super.start_millisec)/(double)index;
				sb.append(" That took: ")
				.append(StringUtils.niceDuration(now - super.start_millisec))
				.append(". Remain: ")
				.append(StringUtils.niceDuration((long)((this.n_expected-(index))*millisec_per_item)))
				;
				}
			final String msg= sb.toString();
			super.emit(msg);
			super.last_messsage_millisec =now;
			}
		}



}
