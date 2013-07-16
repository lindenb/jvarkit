package com.github.lindenb.jvarkit.util.vcf;

import java.io.IOException;

import javax.xml.stream.XMLEventReader;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.variantcontext.VariantContext;

public class XMLVcfReader
	{
	private XMLEventReader reader;
	
	public FeatureCodecHeader readHeader()
			throws IOException
		{
		return null;
		}
	public VariantContext next() throws IOException
		{
		return null;
		}
	}
