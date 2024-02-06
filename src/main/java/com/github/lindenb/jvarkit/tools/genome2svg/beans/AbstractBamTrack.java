package com.github.lindenb.jvarkit.tools.genome2svg.beans;

import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public abstract class AbstractBamTrack extends Track {
	private Predicate<SAMRecord> samFilter = R->true;
	private List<BamBean> bamList  = Collections.emptyList();
	
	public Predicate<SAMRecord> getSamFilter() {
		return this.samFilter;
		}
	public void setSamFilter(Predicate<SAMRecord> samFilter) {
		this.samFilter = samFilter;
		}
		

	}
