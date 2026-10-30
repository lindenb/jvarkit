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


/* base class for Progress and ProgressDict */
public class AbstractProgress {
	private static Logger LOG=Logger.of(AbstractProgress.class);
	
	protected long startMilliSec = System.currentTimeMillis();
	protected long last_messsage_millisec = -1L;
	private Logger logger = LOG;
	private boolean enabled=true;

	protected AbstractProgress(final Logger log) {
		this.logger=LOG;
	}
	

public final AbstractProgress setLogger(Logger logger) {
	this.logger = logger;
	return this;
	}

public AbstractProgress setEnabled(boolean enabled) {
	this.enabled = enabled;
	return this;
	}
public boolean isEnabled() {
	return enabled;
	}

/** return the logger, never null */
public final Logger getLogger() {
	return this.logger==null?LOG:this.logger;
	}



}
