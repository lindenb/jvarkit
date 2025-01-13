/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.hts;

import java.net.URL;
import java.nio.file.Path;
import java.util.Optional;


import htsjdk.samtools.util.FileExtensions;

/**
 * common HTS format
 * @author lindenb
 *
 */
public enum HtsFileType {
	FASTA,
	FASTA_INDEX,
	DICT,
	BAM,
	CRAM,
	SAM,
	VCF,
	COMPRESSED_VCF,
	BCF,
	COMPRESSED_INTERVAL_LIST,
	INTERVAL_LIST,
	WIG,
	BIGWIG,
	BED,
	COMPRESSED_BED
	;
	
	public boolean hasDictionary() {
		switch(this) {
		case FASTA:
		case FASTA_INDEX:
		case DICT:
		case BAM: case CRAM: case SAM:
		case VCF: case COMPRESSED_VCF: case BCF:
		case INTERVAL_LIST: case COMPRESSED_INTERVAL_LIST:
			return true;
		default: return false;
		}
	}
	
	public boolean isTypeFor(URL url) {
		 final Optional<HtsFileType> opt=of(url);
		 if(!opt.isPresent()) return false;
		 return opt.get().equals(this);
		}
	
	public boolean isTypeFor(Path path) {
		 final Optional<HtsFileType> opt=of(path);
		 if(!opt.isPresent()) return false;
		 return opt.get().equals(this);
		}
	
	public boolean isTypeFor(String fname) {
		 final Optional<HtsFileType> opt=of(fname);
		 if(!opt.isPresent()) return false;
		 return opt.get().equals(this);
		}
	
	public static Optional<HtsFileType> of(final URL url) {
		return of(url.getPath());
		}
	
	public static Optional<HtsFileType> of(final Path path) {
		return of(path.getFileName().toString());
		}
	
	public static Optional<HtsFileType> of(final String fname) {
    	if(FileExtensions.FASTA.stream().anyMatch(EXT->fname.endsWith(EXT))) {
    		return Optional.of(FASTA);
    		}
    	else if(fname.endsWith(FileExtensions.BCF)) {
    		return Optional.of(BCF);
    		}
    	else if(fname.endsWith(".vcf.bgz") || fname.endsWith(FileExtensions.COMPRESSED_VCF)) {
    		return Optional.of(COMPRESSED_VCF);
    		}
    	else if(fname.endsWith(FileExtensions.VCF)) {
    		return Optional.of(VCF);
    		}
    	else if(fname.endsWith(FileExtensions.DICT)) {
    		return Optional.of(DICT);
    		}
    	else if(fname.endsWith(FileExtensions.FASTA_INDEX)) {
    		return Optional.of(FASTA_INDEX);
    		}
    	else if(fname.endsWith(FileExtensions.CRAM))  {
    		return Optional.of(CRAM);
    		}
    	else if(fname.endsWith(FileExtensions.BAM))  {
    		return Optional.of(BAM);
    		}
    	if(fname.endsWith(FileExtensions.SAM))  {
    		return Optional.of(SAM);
    		}
    	if(fname.endsWith(FileExtensions.INTERVAL_LIST))  {
    		return Optional.of(INTERVAL_LIST);
    		}
    	if(fname.endsWith(FileExtensions.BED))  {
    		return Optional.of(BED);
    		}
    	if(fname.endsWith(".bed.gz") || fname.endsWith(".bed.bgz") || fname.endsWith(".bed.bgzf"))  {
    		return Optional.of(COMPRESSED_BED);
    		}
    	if(fname.endsWith(FileExtensions.COMPRESSED_INTERVAL_LIST))  {
    		return Optional.of(COMPRESSED_INTERVAL_LIST);
    		}
    	if(fname.toLowerCase().endsWith(".bigwig") || fname.toLowerCase().endsWith(".bw")) {
    		return Optional.of(BIGWIG);
    		}
    	if(fname.toLowerCase().endsWith(".wig") || fname.toLowerCase().endsWith(".wiggle")) {
    		return Optional.of(WIG);
    		}
    	return Optional.empty();
    	}
}
