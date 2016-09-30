/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

import com.github.lindenb.jvarkit.io.IOUtils;

public class Biostar214299 extends AbstractBiostar214299
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Biostar214299.class);
	
	
	private static class Position
		{
		//String contig;
		int refpos;
		Map<Character,String> base2sample = new HashMap<>();
		}
	
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		if(super.positionFile==null) {
			return wrapException("option -"+OPTION_POSITIONFILE+ " is not defined.");
			}
		final String UNAFFECTED_SAMPLE="UNAFFECTED";
		final String AMBIGOUS_SAMPLE="AMBIGOUS";
		final String UNMAPPED="UNMAPPED";
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		final IntervalTreeMap<Position> positionsTreeMap = new IntervalTreeMap<>();
		final Set<String> samples = new HashSet<>();
		try
			{
			sfr = openSamReader(inputName);
			final SAMFileHeader header=sfr.getFileHeader();
			final SAMSequenceDictionary dict =header.getSequenceDictionary();
			if(dict==null) return wrapException("Dictionary missing in input sam");
			
			
			try ( BufferedReader br = IOUtils.openFileForBufferedReading(super.positionFile)) {
				String line;
				while((line=br.readLine())!=null) {
					if(line.trim().isEmpty() || line.startsWith("#")) continue;
					final String tokens[]=line.split("[\t]");
					if(tokens.length<4) return wrapException("Not enough columns in "+line);
					
					final String contig = tokens[0];
					if(dict.getSequence(contig)==null) 
						return wrapException("No such contig in input's sam dictionary: "+contig);
					final int refpos = Integer.parseInt(tokens[1]);
					final Interval interval = new Interval(contig, refpos, refpos);
					Position position = positionsTreeMap.get(interval);
					if(position==null) {
						position = new Position();
						//position.contig = contig;
						position.refpos = refpos;
						 positionsTreeMap.put(interval, position);
					}
					
					final String bases = tokens[2].toUpperCase();
					if(bases.length()!=1 || !bases.matches("[ATGC]"))
						return wrapException("in "+line+" bases should be one letter an ATGC");
					if(position.base2sample.containsKey(bases)) {
						return wrapException("in "+line+" bases already defined for this position");
					}
					
					final String sampleName = tokens[3].trim();
					if(sampleName.isEmpty())
						return wrapException("sample name cannot be empty");
					samples.add(sampleName);
					position.base2sample.put(bases.charAt(0), sampleName);
					
					}
			} catch (final IOException err) {
				return wrapException(err);
			}
			
			if(samples.contains(UNAFFECTED_SAMPLE)) 
				return wrapException("Sample cannot be named "+UNAFFECTED_SAMPLE);
			if(samples.contains(AMBIGOUS_SAMPLE)) 
				return wrapException("Sample cannot be named "+AMBIGOUS_SAMPLE);
			if(samples.contains(UNMAPPED)) 
				return wrapException("Sample cannot be named "+UNMAPPED);

			samples.add(UNAFFECTED_SAMPLE);
			samples.add(AMBIGOUS_SAMPLE);
			samples.add(UNMAPPED);
			
			
			
			final SAMFileHeader newHeader = new SAMFileHeader();
			newHeader.setSortOrder(header.getSortOrder());
			newHeader.setSequenceDictionary(dict);
			newHeader.addComment("generated with "+getName()+" "+getVersion()+" "+getAuthorName()+": "+getProgramCommandLine());
			/* create groups */
			for(final String sample: samples) {
				final SAMReadGroupRecord rg = new SAMReadGroupRecord(sample);
				rg.setSample(sample);
				rg.setLibrary(sample);
				newHeader.addReadGroup(rg);
			}
			
			
			sfw = openSAMFileWriter(newHeader, true);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			final SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				rec.setAttribute("RG",null);
				if(rec.getReadUnmappedFlag()) {
					rec.setAttribute("RG",UNMAPPED);
					sfw.addAlignment(rec);
					continue;
				}
				
				final Cigar cigar=rec.getCigar();
				final Collection<Position> snps =  positionsTreeMap.getContained(new Interval(rec.getContig(),rec.getUnclippedStart(),rec.getUnclippedEnd()));
				if(snps== null || snps.isEmpty()) {
					rec.setAttribute("RG",UNAFFECTED_SAMPLE);
					sfw.addAlignment(rec);
					continue;
				}
				final Map<Integer,Position> index2pos= snps.stream().collect(Collectors.toMap(P->P.refpos,P->P));
				final Set<String> selectedSamples = new HashSet<>();
				
				final byte bases[] =rec.getReadBases();
				if(bases==null) return wrapException("Bases missing in read "+rec);
				
				int refPos1=rec.getUnclippedStart();
				int readPos0=0;
			
				for(final CigarElement ce:cigar.getCigarElements())
					{
					final CigarOperator op = ce.getOperator();
					if(op.consumesReferenceBases() && op.consumesReferenceBases()) 
						{
						for(int i=0;i< ce.getLength();++i){
							final int nowRefPos1 = (refPos1+i);
							final int nowReadPos0 = (readPos0+i);
							final Position position = index2pos.get(nowRefPos1);
							if(position==null) continue;
							final char base = (char)Character.toUpperCase(bases[nowReadPos0]);
							final String sample = position.base2sample.get(base);
							if(sample==null) continue;
							selectedSamples.add(sample);
							
							index2pos.remove(nowRefPos1);
							if(index2pos.isEmpty()) break;
							
							}
						}
					if(op.consumesReferenceBases()) refPos1+=ce.getLength();
					if(op.consumesReadBases()) readPos0+=ce.getLength();
					}
				if(selectedSamples.isEmpty())  {
					rec.setAttribute("RG",UNAFFECTED_SAMPLE);
				} else if(selectedSamples.size()==1) {
					rec.setAttribute("RG",selectedSamples.iterator().next());
				} else
				{
					rec.setAttribute("RG",AMBIGOUS_SAMPLE);
				}
				
				sfw.addAlignment(rec);
				}
			progress.finish();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar214299().instanceMainWithExit(args);

	}

}
