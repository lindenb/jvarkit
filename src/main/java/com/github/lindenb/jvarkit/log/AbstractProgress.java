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

import htsjdk.samtools.util.StopWatch;

/* base class for Progress and ProgressDict */
public class AbstractProgress {
	private static Logger LOG=Logger.of(AbstractProgress.class);
	
	private final Consumer<String> log_emiter;
	protected final long start_millisec;
	protected long last_messsage_millisec;
	private boolean enabled=true;

	protected AbstractProgress(final Consumer<String> log_emiter) {
		this.log_emiter=log_emiter;
		this.start_millisec = System.currentTimeMillis();
		this.last_messsage_millisec = this.start_millisec;
		}
	
	public void setEnabled(boolean enabled) {
		this.enabled = enabled;
		}
	public boolean isEnabled() {
		return enabled;
		}
	
	protected int getEveryXSeconds() {
		return 60;
		}
	
	/** return the logger, never null */
	protected void emit(final String msg) {
		if(!isEnabled())return;
		if(msg==null) return;
		if(this.log_emiter==null) {
			LOG.info(msg);
			}
		else
			{
			this.log_emiter.accept(msg);
			}
		}
	@Override
	public String toString() {
		return StringUtils.niceDuration(System.currentTimeMillis()-this.start_millisec);
		}
	}
