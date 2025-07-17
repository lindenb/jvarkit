/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;

/**
BEGIN_DOC

## Cited In

  * Suppl. material of: "The Mobile Element Locator Tool (MELT): population-scale mobile element discovery and biology". [http://genome.cshlp.org/content/27/11/1916.short](http://genome.cshlp.org/content/27/11/1916.short) Gardner et al. Genome Res. 2017. 27: 1916-1929 
  * Molly M McDonough, Lillian D Parker, Nancy Rotzel McInerney, Michael G Campana, Jesus E Maldonado; Performance of commonly requested destructive museum samples for mammalian genomic studies, Journal of Mammalogy, , gyy080, https://doi.org/10.1093/jmammal/gyy080
  * A long noncoding RNA promotes parasite differentiation in African trypanosomes. Fabien GueganK & al. Sci. Adv., 8 (24), eabn2706. . DOI: 10.1126/sciadv.abn2706
  * Baaijens, J.A., Zulli, A., Ott, I.M. et al. Lineage abundance estimation for SARS-CoV-2 in wastewater using transcriptome quantification techniques. Genome Biol 23, 236 (2022). https://doi.org/10.1186/s13059-022-02805-9
  * The Mitoribosome. Methods in Molecular Biology. https://doi.org/10.1007/978-1-0716-3171-3
  * Makamure, C.E., Justinen, S., Martínez, D.E. et al. Cool temperature inhibits binary fission and results in phenotypic and transcriptomic changes that suggest inducible aging in Diadumene lineata. BMC Res Notes 18, 293 (2025). https://doi.org/10.1186/s13104-025-07378-x

## Example

```bash
$  java -jar dist/jvarkit.jar biostar84452 samtools-0.1.18/examples/toy.sam > out.sam

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
	keywords={"sam","bam","clip"},
	modificationDate = "20250408",
	jvarkit_amalgamion =  true,
	menu="Biostars"
	)
public class Biostar84452 extends OnePassBamLauncher
	{
	private static final Logger LOG = Logger.of(Biostar84452.class);
	
	@Parameter(names={"-t","--tag"},description="tag to flag samrecord as processed")
	private String customTag=null;
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam() {
		if(!StringUtil.isBlank(this.customTag))
			{
			if(customTag.length()!=2 || !customTag.startsWith("X"))
				{
				LOG.error("Bad tag: expect length=2 && start with 'X'");
				return -1;
				}
			}
		return super.beforeSam();
		}
	
	
	@Override
	protected Function<SAMRecord, List<SAMRecord>> createSAMRecordFunction() {
		return rec->{
			if(rec.getReadUnmappedFlag())
				{
				return Collections.singletonList(rec);
				}
			
			final Cigar cigar=rec.getCigar();
			if(cigar==null)
				{
				return Collections.singletonList(rec);
				}
			final String originalCigarSting = rec.getCigarString();
			final byte bases[]= rec.getReadBases();
			if(bases==null || SAMRecord.NULL_SEQUENCE.equals(bases))
				{
				return Collections.singletonList(rec);
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
				return Collections.singletonList(rec);
				}
			if(!StringUtil.isBlank(this.customTag)) rec.setAttribute(this.customTag,originalCigarSting);
			rec.setCigar(new Cigar(L));
			rec.setReadBases(nseq.toByteArray());
			if(quals.length!=0)  rec.setBaseQualities(nqual.toByteArray());
			
			return Collections.singletonList(rec);
			};
		}
		
	public static void main(final String[] args)
		{
		new Biostar84452().instanceMainWithExit(args);
		}

	}
