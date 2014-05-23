package com.github.lindenb.jvarkit.util.picard;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.logging.Logger;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SamFileReaderFactory
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private SamReaderFactory delegate;
	private static ValidationStringency DEFAULT_VALIDATION_STRINGENCY=ValidationStringency.LENIENT;
	public SamFileReaderFactory()
		{
		this.delegate=SamReaderFactory.make();
		this.delegate.disable(SamReaderFactory.Option.values());
		this.delegate.validationStringency(SamFileReaderFactory.DEFAULT_VALIDATION_STRINGENCY);
		}
	
	public static ValidationStringency getDefaultValidationStringency() {
		return DEFAULT_VALIDATION_STRINGENCY;
		}
	public static void setDefaultValidationStringency(
			ValidationStringency defaultValidationStringency) {
		SamFileReaderFactory.DEFAULT_VALIDATION_STRINGENCY = defaultValidationStringency;
		}	
	
	public static SamFileReaderFactory mewInstance()
		{
		return new SamFileReaderFactory();
		}
	
	public SamFileReaderFactory stringency(ValidationStringency validationStringency)
		{
		getDelegate().validationStringency(validationStringency);
		return this;
		}
	

	protected SamReaderFactory getDelegate() {
		return delegate;
		}

	/**
	 *  s can be a 
	 * + URL
	 * + File
	 */
	public SamReader open(String s)
		{
		if(IOUtils.isRemoteURI(s))
			{
			try
				{
				LOG.info("opening url: "+s);
				return getDelegate().open(SamInputResource.of(new URL(s)));
				}
			catch(IOException err)
				{
				throw new PicardException("Cannot open URL \""+s+"\"", err);
				}
			}
		else
			{
			return open(new File(s));
			}
		}

	
	public SamReader open(File samFile)
		{
		LOG.info("opening file "+samFile);
		return getDelegate().open(samFile);
		}
	
	public SamReader open(InputStream in)
		{
		LOG.info("opening stream");
		return getDelegate().open(SamInputResource.of(in));
		}
	
	public SamReader openStdin()
		{
		return open(System.in);
		}
	
	
	}
