package com.github.lindenb.jvarkit.util.picard;

import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.picard.fastq.FastqRecord;

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
