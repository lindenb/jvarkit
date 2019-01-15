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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.List;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
BEGIN_DOC

## Cited In

  * Suppl. material of: "The Mobile Element Locator Tool (MELT): population-scale mobile element discovery and biology". [http://genome.cshlp.org/content/27/11/1916.short](http://genome.cshlp.org/content/27/11/1916.short) Gardner et al. Genome Res. 2017. 27: 1916-1929 
  * Molly M McDonough, Lillian D Parker, Nancy Rotzel McInerney, Michael G Campana, Jesús E Maldonado; Performance of commonly requested destructive museum samples for mammalian genomic studies, Journal of Mammalogy, , gyy080, https://doi.org/10.1093/jmammal/gyy080

## Example

```bash
$  java -jar dist/biostar84452.jar samtools-0.1.18/examples/toy.sam > out.sam

@HD	VN:1.4	SO:unsorted
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
@PG	ID:0	PN:com.github.lindenb.jvarkit.tools.biostar.Biostar84452	VN:b5ebf67dd2926d8a6afadb4d1e36a4959508057f	CL:samtools-0.1.18/examples/toy.sam
(...)
r002	0	ref	9	0	2I6M1P1I1P1I4M2I	*	0	0	AAAGATAAGGGATAAA	*
(...)


$ grep r002 samtools-0.1.18/examples/toy.sam
r002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*

```
## See also

* https://twitter.com/EugenomeUK/status/938031803612491776

END_DOC




 */
@Program(name="biostar84452",
	biostars=84452,
	description="remove clipped bases from a BAM file",
	keywords={"sam","bam","clip"}
	)
public class Biostar84452 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar84452.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@ParametersDelegate
	private WritingBamArgs WritingBamArgs=new WritingBamArgs();
	@Parameter(names={"-t","--tag"},description="tag to flag samrecord as processed")
	private String customTag=null;
	
	@Override
	public int doWork(final List<String> args) {
		if(!StringUtil.isBlank(this.customTag))
			{
			if(customTag.length()!=2 || !customTag.startsWith("X"))
				{
				LOG.error("Bad tag: expect length=2 && start with 'X'");
				return -1;
				}
			}
		SAMFileWriter sfw=null;
		SamReader sfr=null;
		try
			{
			sfr = super.openSamReader(oneFileOrNull(args));
			SAMFileHeader header=sfr.getFileHeader();

			sfw = this.WritingBamArgs.openSAMFileWriter(outputFile, header, true);
			
			
			long nChanged=0L;
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progress.watch(iter.next());
				if(rec.getReadUnmappedFlag())
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				final Cigar cigar=rec.getCigar();
				if(cigar==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				final String originalCigarSting = rec.getCigarString();
				final byte bases[]= rec.getReadBases();
				if(bases==null)
					{
					sfw.addAlignment(rec);
					continue;
					}
				
				final ArrayList<CigarElement> L=new ArrayList<CigarElement>();
				final ByteArrayOutputStream nseq=new ByteArrayOutputStream();
				final ByteArrayOutputStream nqual=new ByteArrayOutputStream();
				
				final byte quals[]= rec.getBaseQualities();
				int indexBases=0;
				for(final CigarElement ce:cigar.getCigarElements())
					{
					switch(ce.getOperator())
						{
						case S: indexBases+=ce.getLength(); break;
						case H: //cont
						case P: //cont
						case N: //cont
						case D:
							{
							L.add(ce);
							break;
							}
							
						case I:
						case EQ:
						case X:
						case M:
							{
							L.add(ce);
							nseq.write(bases,indexBases,ce.getLength());
							if(quals.length!=0) nqual.write(quals,indexBases,ce.getLength());
							indexBases+=ce.getLength(); 
							break;
							}
						default:
							{
							throw new SAMException("Unsupported Cigar opertator:"+ce.getOperator());
							}
						}
					
					}
				if(indexBases!=bases.length)
					{
					throw new SAMException("ERRROR "+rec.getCigarString());
					}
				if(L.size()==cigar.numCigarElements())
					{
					sfw.addAlignment(rec);
					continue;
					}
				++nChanged;
				if(!StringUtil.isBlank(this.customTag)) rec.setAttribute(this.customTag,originalCigarSting);
				rec.setCigar(new Cigar(L));
				rec.setReadBases(nseq.toByteArray());
				if(quals.length!=0)  rec.setBaseQualities(nqual.toByteArray());
				sfw.addAlignment(rec);
				}
			progress.finish();
			iter.close();
			LOG.info("Num records changed:"+nChanged);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		return 0;
		}
	
	public static void main(final String[] args)
		{
		new Biostar84452().instanceMainWithExit(args);
		}

	}
