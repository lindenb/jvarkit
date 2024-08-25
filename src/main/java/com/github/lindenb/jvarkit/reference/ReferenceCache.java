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
package com.github.lindenb.jvarkit.reference;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**
 * manage local and remote fasta references
 * @author lindenb
 *
 */
public class ReferenceCache implements Closeable {
private static Logger LOG=Logger.build(ReferenceCache.class).make();
private final File cacheDir;
private final boolean deleteOnClose;
private final List<AbstractInfo> references=new ArrayList<>();


private static abstract class AbstractInfo {
	SAMSequenceDictionary dict=null;
	Path  path;
	abstract SAMSequenceDictionary getDictionary();
	abstract Path getPath();
	abstract boolean isLocal();
	}

private static class PathInfo extends AbstractInfo {
	PathInfo(final Path p) {
		super.path=p;
		IOUtil.assertFileIsReadable(p);
		}
	@Override
	boolean isLocal() {
		return true;
		}
	@Override
	SAMSequenceDictionary getDictionary() {
		if(dict==null) {
			dict=SAMSequenceDictionaryExtractor.extractDictionary(this.path);
			if(dict==null) throw new JvarkitException.FastaDictionaryMissing(this.path);
			}
		return dict;
		}
	@Override
	Path getPath() {
		return super.path;
		}
	}

private class URLInfo extends AbstractInfo {
	final URL url;
	URLInfo(URL url) {
		this.url=url;
		}
	@Override
	boolean isLocal() {
		return false;
		}
	@Override
	SAMSequenceDictionary getDictionary() {
		if(dict==null) {
			try {
				dict=new SequenceDictionaryExtractor().extractRequiredDictionary(this.url);
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		return dict;
		}
	@Override
	Path getPath() {
		Objects.requireNonNull(this.dict);
		if(this.path!=null) return this.path;
		final String md5= StringUtils.md5(this.url.toString());
		try {
			final Path fasta= new File(cacheDir,md5+".fa").toPath();
			if(Files.exists(fasta)) throw new IllegalStateException("file exists "+fasta);
			LOG.info("download "+this.url+" to "+fasta);

			try(InputStream in=IOUtils.mayBeGzippedInputStream(ParsingUtils.getURLHelper(url).openInputStream())) {
				IOUtils.copyTo(in, fasta);
				}
			catch(IOException err) {
				Files.deleteIfExists(fasta);
				throw err;
				}
			
			LOG.info("create fai for "+fasta);
			FastaSequenceIndexCreator.buildFromFasta(fasta);
				
			
			final Path dictf = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta);
			LOG.info("write dict for "+this.url+" to "+dictf);
	        try(Writer w=Files.newBufferedWriter(dictf)) {
		        new SAMSequenceDictionaryCodec(w).encode(this.getDictionary());
		        w.flush();
	        	}
			catch(IOException err) {
				Files.deleteIfExists(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta));
				Files.deleteIfExists(fasta);
				Files.deleteIfExists(dictf);
				throw err;
				}
	        
	        try(Writer w=Files.newBufferedWriter(new File(cacheDir,md5+".readme").toPath())) {
	        	w.write(fasta.toString()+" downloaded from "+this.url+"\n");
	        	w.flush();
	        	}
	        
	        
	        this.path=fasta;
	        PathInfo pi=new PathInfo(this.path);
	        pi.dict=this.dict;
	        references.add(pi);
	        references.remove(this);
			return this.path;
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	}


public ReferenceCache() {
	this.cacheDir= IOUtil.createTempDir("ref.cache.", ".dir");
	this.deleteOnClose=true;
	}

public ReferenceCache(final File dir) {
	this.cacheDir= dir;
	this.deleteOnClose=false;
	for(File fa: dir.listFiles(new FileFilter() {
		@Override
		public boolean accept(File f) {
			return FileExtensions.FASTA.stream().anyMatch(ext->f.exists() && f.isFile() && f.getName().endsWith(ext));
			}
		}))
		{
		this.references.add(new PathInfo(fa.toPath()));
		}
	
	}

public ReferenceCache loadTable(Path p) throws IOException {
	try(BufferedReader br=IOUtils.openPathForBufferedReading(p)) {
		String line=br.readLine();
		if(line==null) throw new IOException("cannot read first line of "+p);
		final FileHeader fh=new FileHeader(line, CharSplitter.TAB);
		fh.assertColumnExists("fasta");
		while((line=br.readLine())!=null) {
			final FileHeader.RowMap row= fh.toMap(line);
			final String fasta=row.get("fasta");
			if(IOUtil.isUrl(fasta)) {
				this.references.add(new URLInfo(new URL(fasta)));
				}
			else
				{
				this.references.add(new PathInfo(Paths.get(fasta)));
				}
			}
		}
	return this;
	}


protected List<AbstractInfo> getReferences() {
	return this.references;
	}

public synchronized Path getRequiredReferenceByDict(final SAMSequenceDictionary dict) {
	final Optional<Path> opt = getReferenceByDict(dict);
	if(!opt.isPresent()) throw new SAMException("Cannot find dictionary for dictionary:\n"+dict.getSequences().stream().limit(3L).
			map(SSR->"\t*  $"+SSR.getSequenceIndex()+" "+SSR.getSequenceName()+" "+SSR.getSequenceLength()).collect(Collectors.joining("\n"))+"\n\t *  ...");
	return opt.get();
	}


public synchronized Optional<Path> getReferenceByDict(final SAMSequenceDictionary dict) {
	Objects.requireNonNull(dict);
	/* first we search local path */
	Optional<Path> local= getReferences().
			stream().
			filter(it->it.isLocal()).
			filter(it->SequenceUtil.areSequenceDictionariesEqual(it.getDictionary(), dict)).
			map(it->it.getPath()).
			findFirst()
			;
	if(local.isPresent()) return local;
	/* second we search remote path */
	return getReferences().
		stream().
		filter(it->!it.isLocal()).
		filter(it->SequenceUtil.areSequenceDictionariesEqual(it.getDictionary(), dict)).
		map(it->it.getPath()).
		findFirst()
		;
	}

@Override
public void close() {
	if(!deleteOnClose) IOUtil.deleteDirectoryTree(cacheDir);
	}
}
