package com.github.lindenb.jvarkit.tools.genome2svg.beans;


import java.lang.ref.SoftReference;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.RuntimeIOException;

public class BamBean extends PathBean {
	private ReferenceBean reference;
	private SoftReference<SAMFileHeader> samFileHeader = new SoftReference<>(null);
	
	public void setReference(ReferenceBean reference) {
		this.reference = reference;
		}
	public ReferenceBean getReference() {
		return reference;
		}
	
	@Override
	public SAMSequenceDictionary getSAMSequenceDictionary() {
		return getSAMFileHeader().getSequenceDictionary();
		}

	public SAMFileHeader getSAMFileHeader() {
		SAMFileHeader hdr = this.samFileHeader.get();
		if(hdr==null) {
			try(SamReader sr =  openSamReader()) {
				hdr = sr.getFileHeader();
				if(hdr==null) throw new IllegalStateException("Cannot get Header from "+getPath());
				this.samFileHeader =  new SoftReference<SAMFileHeader>(hdr);
				}
			catch(Throwable err) {
				throw new RuntimeIOException(err);
				}
			}
		return hdr;
		}

	
	@Override
	public String resolveContig(String ctg) {
		if(StringUtils.isBlank(ctg)) return null;
		return ContigNameConverter.fromOneDictionary(getSAMSequenceDictionary()).apply(ctg);
		}
	
	public SamReader openSamReader() {
		final SamReaderFactory srf = SamReaderFactory.make();
		srf.validationStringency(ValidationStringency.LENIENT);
		if(getReference()!=null) srf.referenceSequence(getReference().asPath());
		final SamReader sr =  srf.open(SamInputResource.of(getPath()));
		if(samFileHeader.get()==null) {
			samFileHeader = new SoftReference<>(sr.getFileHeader());
			}
		return sr;
		}
	
	public String getSampleName() {
		return  getSAMFileHeader().getReadGroups().
				stream().
				map(RG->RG.getSample()).
				filter(S->!StringUtils.isBlank(S)).
				findFirst().
				orElse(IOUtils.getFilenameWithoutCommonSuffixes(asPath()));
		}

}
