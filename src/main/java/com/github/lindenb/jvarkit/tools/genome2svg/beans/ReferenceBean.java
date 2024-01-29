package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.lang.ref.SoftReference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class ReferenceBean extends PathBean {
	private SoftReference<SAMSequenceDictionary> softDicr = new SoftReference<>(null);

	private String ucscName  = null;
	public ReferenceBean() {
		}
	
	@Override
	public SAMSequenceDictionary getSAMSequenceDictionary() {
		SAMSequenceDictionary dict = softDicr.get();
		if(dict==null) {
			dict=SAMSequenceDictionaryExtractor.extractDictionary(asPath());
			if(dict==null) throw new NullPointerException("cannot get dict from "+getPath());
			softDicr = new SoftReference<>(dict);
			}
		return dict;
		}

	
	public String getHyperlink(Locatable loc) {
		return null;
		}
	}
