package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.Closeable;
import java.io.File;
import java.util.AbstractList;
import java.util.NoSuchElementException;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

public class IndexedGenome
extends AbstractList<GenomicSequence>
implements Closeable {

private static final Logger LOG=Logger.build().prefix("IndexedGenome").make();
public static final String ENV_NAME="REF_PATH";
private final IndexedFastaSequenceFile indexedFastaSequenceFile;
private final SAMSequenceDictionary dict;
private GenomicSequence cached=null;

private IndexedGenome(final IndexedFastaSequenceFile indexedFastaSequenceFile) {
	this.indexedFastaSequenceFile=indexedFastaSequenceFile;
	this.dict=this.indexedFastaSequenceFile.getSequenceDictionary();
	}
public GenomicSequence  get(int i)
	{
	return _get(this.dict.getSequence(i).getSequenceName());
	}

private GenomicSequence  _get(final String contigName)
	{
	final SAMSequenceRecord rec=this.dict.getSequence(contigName);
	if(rec==null) throw new NoSuchElementException(contigName);
	if(this.cached!=null && this.cached.getChrom().equals(contigName)) {
		return this.cached;
		}
	this.cached= new GenomicSequence(this.indexedFastaSequenceFile,rec.getSequenceName());
	return this.cached;
	}

@Override
public int size() {
	return this.dict.size();
	}

@Override
public void close() {
	CloserUtil.close(this.indexedFastaSequenceFile);
	}
public static class Builder
	{
	private File genomeFile;
	private String envName=ENV_NAME;
	
	public Builder setGenomeFile(final File genomeFile) {
		this.genomeFile = genomeFile;
		return this;
		}
	
	public Builder setEnvName(final String envName) {
		this.envName = envName;
		return this;
		}
	
	public IndexedGenome make()
		{
		File f = this.genomeFile;
		if(f==null && envName!=null && !envName.trim().isEmpty()) {
			LOG.info("searching reference genome using getenv(\""+envName+"\")");
			String ref=System.getenv(envName);
			if(ref!=null) 
				{
				f=new File(ref);
				}
			}
		if(f==null) throw new JvarkitException.UserError("REF genome is undefined");
		IOUtil.assertFileIsWritable(f);
		
		IndexedFastaSequenceFile idx=null;
		try
			{
			idx = new IndexedFastaSequenceFile(f);
			}
		catch(final Exception err)
			{
			LOG.severe(err);
			throw new RuntimeIOException(err);
			}
		final SAMSequenceDictionary dict=idx.getSequenceDictionary();
		if(dict==null)
			{
			CloserUtil.close(idx);
			throw new JvarkitException.FastaDictionaryMissing(f.getPath());
			}

		final IndexedGenome g=new IndexedGenome(idx);
		return g;
		}
	}

public static Builder build() {
	return new Builder();
	}
}
