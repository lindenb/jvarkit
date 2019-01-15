/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

package com.github.lindenb.jvarkit.util.samtools;

import java.util.Collection;
import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;

/** how to group a read */
public enum SAMRecordPartition
	implements Function<SAMReadGroupRecord,String>
	{
	readgroup, 
	sample, 
	library, 
	platform, 
	center, 
	sample_by_platform,
	sample_by_center,
	sample_by_platform_by_center,
	any
	;
	
	public static final String OPT_DESC="Data partitioning using the SAM Read Group (see https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can be any combination of sample, library.... ";
	
	
	/** return the label of the partition for this SAMRecord, result can be null */
	public String getPartion(final SAMRecord rec,final String defaultValue)	{
		if(this.equals(any)) return "all";
		return rec==null?null:this.apply(rec.getReadGroup(),defaultValue);
	}	

	/** return the label of the partition for this SAMRecord, result can be null */
	public final String getPartion(final SAMRecord rec)	{
		return getPartion(rec,null);
	}	
	
	/** return the label of the partition for this read-group, result is defaultValue if partition is null or empty*/
	public final String apply(final SAMReadGroupRecord rg,final String defaultValue)
		{
		final String p = this.apply(rg);
		return StringUtil.isBlank(p)?defaultValue:p;
		}
	
	/** return the label of the partition for this read-group, result can be null */
	@Override
	public String apply(final SAMReadGroupRecord rg)
		{
		if(this.equals(any)) return "all";
		if( rg == null) return null;
		switch(this)
			{
			case readgroup : return rg.getId();
			case sample : return  rg.getSample();
			case library : return rg.getLibrary();
			case platform : return rg.getPlatform();
			case center : return rg.getSequencingCenter();
			case sample_by_platform : return rg.getSample()!=null && rg.getPlatform()!=null?
					String.join("_",rg.getSample(),rg.getPlatform()):null;
			case sample_by_center : return rg.getSample()!=null && rg.getSequencingCenter()!=null?
					String.join("_",rg.getSample(),rg.getSequencingCenter()):null;
			case sample_by_platform_by_center : return rg.getSample()!=null && rg.getPlatform()!=null && rg.getSequencingCenter()!=null?
					String.join("_",rg.getSample(),rg.getPlatform(),rg.getSequencingCenter()):null;
			default: throw new IllegalStateException(this.name());
			}
		}
	/** extract all distinct non-null/non-empty labels from a set of read-groups */
	public Set<String> getPartitions(final Collection<SAMReadGroupRecord> rgs) {
		if(rgs==null || rgs.isEmpty()) return Collections.emptySet();
		final Set<String> set = new TreeSet<>();
		for(final SAMReadGroupRecord rg:rgs)
			{
			final String partition = this.apply(rg);
			if(StringUtil.isBlank(partition)) continue;
			set.add(partition);
			}
		return set;
		}
	/** extract all distinct non-null labels from a Sam File header */
	public Set<String> getPartitions(final SAMFileHeader header) {
		return getPartitions(header==null?Collections.emptyList():header.getReadGroups());
		}

	}
