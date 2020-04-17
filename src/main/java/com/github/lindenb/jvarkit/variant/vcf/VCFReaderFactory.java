/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.File;
import java.nio.file.Path;

import htsjdk.variant.vcf.VCFFileReader;

/** Factory creating VCFReader */
public abstract class VCFReaderFactory {
	private boolean requireIndex = true; /* because default is true in VCFFileReader */
	
	/** initialize with requireIndex = true */
	protected VCFReaderFactory() {	
		}
	
	public static VCFReaderFactory makeDefault() {
		return new VCFReaderFactory() {
			};
		}
	
	public VCFReaderFactory setRequireIndex(boolean b) {
		this.requireIndex = b;
		return this;
		}
	
	public boolean isRequireIndex() {
		return requireIndex;
		}
	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public VCFFileReader open(final Path p) {
		return open(p,isRequireIndex());
		}
	
	/** open new VCFReader */
	public VCFFileReader open(final Path path,boolean requireIndex) {
		return new VCFFileReader(path, requireIndex);
		}
	
	/** open new VCFReader */
	public final VCFFileReader open(final File path,boolean requireIndex) {
		return open(path.toPath(),requireIndex);
		}
	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public final VCFFileReader open(final File p) {
		return open(p,isRequireIndex());
		}

}
