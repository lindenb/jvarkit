/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.dict;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Optional;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.hts.HtsFileType;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.variant.vcf.VcfHeaderExtractor;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFHeader;
/**
 * My implementation of htsjdk.variant.utils.SAMSequenceDictionaryExtractor
 * works with BCF, FAI, remote files
 *
 */
public class SequenceDictionaryExtractor  {
	public SequenceDictionaryExtractor() {
		}
	

	
    
    protected InputStream openStream(final URL url) throws IOException {
    	return ParsingUtils.getURLHelper(url).openInputStream();
    	}
   
    
    private Optional<SAMSequenceDictionary> extractDictionary(final HtsFileType type,InputStream in,String source) {
    	try { 
	    	switch(type) {
		    	case VCF:
		    	case COMPRESSED_VCF:
	    		case BCF:
	    			{
	    			final VCFHeader h=VcfHeaderExtractor.decode(in);
	    			return h==null?Optional.empty():Optional.ofNullable(h.getSequenceDictionary());
	    			}
	    		case DICT: {
					try (BufferedLineReader bufferedLineReader = new BufferedLineReader(in)) {
						final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
						final SAMFileHeader header = codec.decode(bufferedLineReader,source);
						return Optional.of(header.getSequenceDictionary());
						} 
					catch (final Throwable e) {
						return Optional.empty();
						}
	    			}
	    		case INTERVAL_LIST:
	    		case COMPRESSED_INTERVAL_LIST:
	    			{
	    			try(BufferedReader br= new BufferedReader(new InputStreamReader(in))) {
	    				final SAMFileHeader h= IntervalList.fromReader(br).getHeader();
	    				return h==null?Optional.empty():Optional.ofNullable(h.getSequenceDictionary());
	    				}
	    			}
	    		case CRAM:
	    			{
    				final CramHeader cramHeader = CramIO.readCramHeader(in);
    				final Optional<SAMFileHeader> samHeader = Optional.ofNullable(
    						Container.readSAMFileHeaderContainer(cramHeader.getCRAMVersion(), in, source));
    				return samHeader.isPresent() ?Optional.ofNullable( samHeader.get().getSequenceDictionary()) :Optional.empty();
	    			}
	    		case FASTA_INDEX:
	    			{
    				try(BufferedReader br= new BufferedReader(new InputStreamReader(in))) {
	    				return Optional.of(new SAMSequenceDictionary(br.lines().
	    						map(L->CharSplitter.TAB.split(L)).
	    						map(A->new SAMSequenceRecord(A[0], Integer.parseInt(A[1]))).
	    						collect(Collectors.toList())));
	    				}
	    			}
	    		default: throw new IllegalArgumentException(type.name());
	    		}
    		}
    	catch(Exception err) {
    		return Optional.empty();
    		}
    	}
    
    
    
    
    private void notFound(String s) {
    	switch(HtsFileType.of(s).get()) {
    	case FASTA: throw new JvarkitException.FastaDictionaryMissing(s);
    	case VCF: case BCF: case COMPRESSED_VCF: throw new JvarkitException.VcfDictionaryMissing(s);
    	default:throw new JvarkitException.DictionaryMissing(s);
    	}
    }
    
    public SAMSequenceDictionary extractRequiredDictionary(final String pathOrUrl) {
    	final Optional<SAMSequenceDictionary>  opt= extractDictionary(pathOrUrl);
    	if(!opt.isPresent()) notFound(pathOrUrl);
    	return opt.get();
    	}

    public SAMSequenceDictionary extractRequiredDictionary(final Path path) {
    	final Optional<SAMSequenceDictionary>  opt= extractDictionary(path);
    	if(!opt.isPresent()) notFound(path.toString());
    	return opt.get();
    	}
    
    public SAMSequenceDictionary extractRequiredDictionary(final URL url) throws IOException {
		final Optional<SAMSequenceDictionary>  opt= extractDictionary(url);
	    if(!opt.isPresent()) notFound(url.getPath());
    	return opt.get();
		}
    
    public Optional<SAMSequenceDictionary> extractDictionary(final String pathOrurl) {
    	if(IOUtil.isUrl(pathOrurl)) {
    		try {
				return extractDictionary(new URL(pathOrurl));
				}
    		catch (MalformedURLException e) {
				throw new IllegalArgumentException(e);
				}
    		catch (IOException e) {
				throw new IllegalArgumentException(e);
				}
    		}
    	else
    		{
    		return extractDictionary(Paths.get(pathOrurl));
    		}
    	}
    
    public Optional<SAMSequenceDictionary> extractDictionary(final Path path) {
    	final Optional<HtsFileType> ot = HtsFileType.of(path);
    	if(!ot.isPresent()) return Optional.empty();
	    try {
	    	switch(ot.get()) {
		    	case SAM:
		    	case CRAM:
		    	case BAM:
	    		case INTERVAL_LIST:
	    		case COMPRESSED_INTERVAL_LIST:
	    		case DICT:
	    		case FASTA: {
	    			return Optional.of(SAMSequenceDictionaryExtractor.extractDictionary(path));
	    			}
	    		case FASTA_INDEX:
	    		case VCF:
	    		case COMPRESSED_VCF:
	    		case BCF:
	    			{
	    			try(InputStream in= Files.newInputStream(path)) {
	    				return extractDictionary(ot.get(),in, path.toString());
	    				}
	    			}
	    		default: throw new IllegalArgumentException("cannot extract dictionary from "+ot.get().name()+":"+path);
	    		}
	    	}
    	catch(MalformedURLException err ) {
    		return Optional.empty();
    		}
    	catch(IOException err) {
    		return Optional.empty();
    		}
    	}

    
    public Optional<SAMSequenceDictionary> extractDictionary(final URL url) throws IOException {
    	final Optional<HtsFileType> ot = HtsFileType.of(url);
    	if(!ot.isPresent()) return Optional.empty();
	    	try {
	    	switch(ot.get()) {
	    		case FASTA: {
	    			final String name = url.getPath();
	    	        final String ext =  FileExtensions.FASTA.stream().
	    	         	filter(name::endsWith).
	    	         	findFirst().
	    	         	get();
	    	        final URL url2=new URL(url,name.substring(0,name.length()-ext.length())+FileExtensions.DICT);
	    	        try {
	    	        	return extractDictionary(url2);
	    	        	}
	    	        catch(java.io.FileNotFoundException err) {
		    	        final URL url3=new URL(url,name+FileExtensions.FASTA_INDEX);
	    	        	return extractDictionary(url3);
	    	        	}
	    			}
	    		case FASTA_INDEX:
	    		case VCF:
	    		case BCF:
	    		case COMPRESSED_VCF:
	    		case DICT:
	    		case CRAM:
	    			{
	    			try(InputStream in=openStream(url)) {
	    				return extractDictionary(ot.get(),in, url.toString());
	    				}
	    			}
	    		case SAM:case BAM:
	    			{
	                try(SamReader sr= SamReaderFactory.makeDefault().open(SamInputResource.of(url))) {
	                	final SAMFileHeader h= sr.getFileHeader();
	    				return h==null?Optional.empty():Optional.ofNullable(h.getSequenceDictionary());
						}
	    			}
	    		case INTERVAL_LIST:
	    		case COMPRESSED_INTERVAL_LIST:
	    			{
		    		try(InputStream in=openStream(url)) {
		    			return extractDictionary(ot.get(),IOUtils.mayBeGzippedInputStream(in), url.toString());
		    			}
	    			}
	    		default: throw new IllegalArgumentException("cannot extract dictionary from "+ot.get().name()+":"+url);
	    		}
	    	}
    	catch(MalformedURLException err ) {
    		throw new IOException(err);
    		}
    	}
}
