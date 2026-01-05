/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.samtools;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.SequenceUtil;

/**
 * generate a Sample and BAM list
 */
public class SampleAndBamFactory {
private ValidationStringency validationStringency = ValidationStringency.LENIENT;
private Function<SAMReadGroupRecord,String> sampleExtractor = RG->RG.getSample();
private boolean enable_multiple_samples_in_one_bam = false;
private boolean enable_filename_if_no_sample = false;
private boolean enable_duplicate_samples = false;
private boolean enable_non_indexed_bam = false;
private boolean enable_no_associated_reference = false;
private final List<RefAndDict> referenceAndDictionaries  = new ArrayList<>();

public static interface SampleAndBamRecord {
	/** returns the sample associated to the bam */
	public String getSampleName();
	/** the Bam for this record */
	public Path getPath();
	/** the Reference */
	public Optional<Path> getReference();
	}

/** enable multiple bam sharing the same sample name */
public SampleAndBamFactory setEnableDuplicateSamples(boolean enable_duplicate_samples) {
	this.enable_duplicate_samples = enable_duplicate_samples;
	return this;
	}
/** enable filename as the sample name if no sample was found in RG */
public SampleAndBamFactory setEnableFilenameIfNoSample(boolean enable_filename_if_no_sample) {
	this.enable_filename_if_no_sample = enable_filename_if_no_sample;
	return this;
	}
/** enable multiple sample in one BAM */
public SampleAndBamFactory setEnableMultipleSamplesInOneBam(boolean enable_multiple_samples_in_one_bam) {
	this.enable_multiple_samples_in_one_bam = enable_multiple_samples_in_one_bam;
	return this;
	}
/** do not check for the existence of an index */
public SampleAndBamFactory setEnableNonIndexedBam(boolean enable_non_indexed_bam) {
	this.enable_non_indexed_bam = enable_non_indexed_bam;
	return this;
	}

/** set sample name extractor */
public SampleAndBamFactory setSampleExtractor(Function<SAMReadGroupRecord, String> sampleExtractor) {
	this.sampleExtractor = sampleExtractor;
	return this;
	}

/** do not throw error if no fasta can be associated to a BAM */
public SampleAndBamFactory setEnableWithoutReference(boolean enable_no_associated_reference) {
	this.enable_no_associated_reference = enable_no_associated_reference;
	return this;
	}

private static class RefAndDict {
	final Path reference;
	final SAMSequenceDictionary dict;
	RefAndDict(final Path reference,final SAMSequenceDictionary dict) {
		this.reference = reference;
		this.dict = dict;
		}
	
	}
/** add reference to collection of known references */
public SampleAndBamFactory reference(final Path fasta) {
	final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(fasta);
	this.referenceAndDictionaries.add(new RefAndDict(fasta, dict));
	return this;
	}

private static class SampleAndBamRecordImpl implements SampleAndBamRecord {
	final String sn;
	final Path path;
	final Path reference;
	SampleAndBamRecordImpl(final String sn,final Path path,final Path reference) {
		this.sn = sn;
		this.path = path;
		this.reference = reference;
		}
	@Override
	public int hashCode() {
		return this.sn.hashCode()*31+ this.path.hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null) return false;
		if(!(obj instanceof SampleAndBamRecord)) return false;
		final SampleAndBamRecord x=SampleAndBamRecord.class.cast(obj);
		return this.getPath().equals(x.getPath()) &&
				this.getSampleName().equals(x.getSampleName());
		}
	@Override
	public Optional<Path> getReference() {
		return Optional.ofNullable(reference);
		}
	@Override
	public String getSampleName() {
		return this.sn;
		}
	@Override
	public Path getPath() {
		return this.path;
		}
	@Override
	public String toString() {
		return getSampleName()+":"+getPath();
		}
	}

private SamReaderFactory createSamReaderFactory() {
	return  SamReaderFactory.make().
			validationStringency(this.validationStringency);
	}

public List<SampleAndBamRecord> parse(final List<Path> paths) throws IOException {
	if(paths.isEmpty()) return Collections.emptyList();
	if(paths.size()==1 && paths.get(0).getFileName().toString().endsWith(".list")) {
		return parse(IOUtils.unrollPath(paths.get(0)));
		}
	
	final List<SampleAndBamRecord> L = new ArrayList<>(paths.size());
	final SamReaderFactory srf = createSamReaderFactory();
	
	for(final Path path:paths) {
		if(L.stream().anyMatch(X->X.getPath().equals(path))) {
			continue;
			}
		Path reference = null;
		try(SamReader sr = srf.open(path)) {
			if(!this.enable_non_indexed_bam) {
				if(!sr.hasIndex()) {
					throw new IOException("Bam file "+path+" is not indexed.");
					}
				}
			final SAMFileHeader header = sr.getFileHeader();

			if(!enable_no_associated_reference) {
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				reference  = this.referenceAndDictionaries.stream().
						filter(SR->SequenceUtil.areSequenceDictionariesEqual(SR.dict, dict)).
						map(SR->SR.reference).
						findFirst().orElseThrow(()->new IOException("Cannot find REFerence genome for "+path));
				}
			
			
			final Set<String> samples = header.getReadGroups().
				stream().
				map(this.sampleExtractor).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toCollection(HashSet::new));
			
			if(samples.size()>1 && !enable_multiple_samples_in_one_bam) {
				throw new IOException("multiple sample "+String.join(",", samples)+" in "+path);
				}
			
			else if(samples.isEmpty()) {
				if(enable_filename_if_no_sample) {
					samples.add(IOUtils.getFilenameWithoutCommonSuffixes(path));
					}
				else
					{
					throw new IOException("no sample in "+path);
					}
				}
			for(final String sampleName: samples) {
				final Optional<SampleAndBamRecord> other =  L.stream().
						filter(SR->SR.getSampleName().equals(sampleName)).
						findAny();
				if(other.isPresent() && !enable_duplicate_samples) {
					throw new IOException("duplicate sample "+ sampleName +" in "+path+" and "+other.get().getPath());
					}
				L.add(new SampleAndBamRecordImpl(sampleName,path,reference));
				}
			}
		}
	Collections.sort(L, (A,B)->A.getSampleName().compareTo(B.getSampleName()));
	return L;
	}
}
