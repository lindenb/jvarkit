package com.github.lindenb.jvarkit.util.picard;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.Closeable;
import java.util.Iterator;




/**
 * the original picard FastqReader doesn't allow empty lines...
 */
public interface FastqReader extends Iterator<FastqRecord>, Closeable
	{
    public void setValidationStringency( ValidationStringency validationStringency);
    public ValidationStringency getValidationStringency();
	}
