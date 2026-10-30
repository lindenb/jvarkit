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


import com.github.lindenb.jvarkit.lang.StringUtils;


/** progress log reloaded... */
public class Progress extends AbstractProgress {
	private static Logger LOG=Logger.of(Progress.class);
	
	protected long n_expected = 0L;
	private long current_index = 0L;

	protected Progress(final Logger log,final long n_expected) {
		super(log);
		this.n_expected=n_expected;
	}
	
public static Progress of(long n_expected) {
	return of(LOG,n_expected);
	}
public static Progress of(Logger log,long n_expected) {
	return new Progress(log,n_expected);
	}


/** end item, with no message */
public void endItem() {
	endItem(null);
	}

/** end item, with message */
public void endItem(final Object optional_msg) {
	int every_x_seconds = 60;
	final long now = System.currentTimeMillis();
	final long diff_millisec = now- super.last_messsage_millisec;
	boolean should_display=	isEnabled() && 
			this.n_expected >0 && 
			this.current_index>=0 &&
			( super.last_messsage_millisec<0L || (diff_millisec / 1_000 <= every_x_seconds));
	if(should_display) {
		final long millisec_per_item  = (now - super.startMilliSec)/(current_index+1L);
		final StringBuilder sb=new StringBuilder();
		sb.append("N=").append(StringUtils.niceInt(this.current_index+1)).append(" / ")
			.append(StringUtils.niceInt(this.n_expected))
			.append(". That took: ")
			.append(StringUtils.niceDuration(System.currentTimeMillis() - super.startMilliSec))
			.append(". Remain: ")
			.append(StringUtils.niceDuration((this.n_expected-(current_index+1))*millisec_per_item))
			;
		if(optional_msg!=null) {
			sb.append(". ").append(optional_msg);
			}
		final String msg= sb.toString();
		getLogger().info(msg);
		super.last_messsage_millisec = now;
		}
	// do not increment for genomic data
	this.current_index++;	
	}



}
