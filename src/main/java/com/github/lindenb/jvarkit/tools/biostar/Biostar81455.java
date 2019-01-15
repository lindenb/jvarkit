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
import java.io.BufferedReader;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Input

tab delimited file. 2 columns: CHROM and POS

## Example

```bash
echo -e "chr22\t41258261\nchr22\t52000000\nchr22\t0" |\
	java   dist/biostar81455.jar 

chr22	41258261	uc003azg.2	41253084	41258785	POSITIVE	Exon 2	41257621	41258785	0
chr22	41258261	uc011aox.2	41253084	41305239	POSITIVE	Exon 1	41253084	41253249	-5012
chr22	41258261	uc003azi.3	41253084	41328823	POSITIVE	Exon 1	41253084	41253249	-5012
chr22	41258261	uc003azj.3	41255553	41258130	NEGATIVE	Exon 1	41255553	41258130	-131
chr22	41258261	uc010gyh.1	41258260	41282519	POSITIVE	Exon 1	41258260	41258683	0
chr22	41258261	uc011aoy.1	41258260	41363888	POSITIVE	Exon 1	41258260	41258683	0
chr22	52000000	uc011asd.2	51195513	51227614	POSITIVE	Exon 4	51227177	51227614	-772386
chr22	52000000	uc003bni.3	51195513	51238065	POSITIVE	Exon 4	51237082	51238065	-761935
chr22	52000000	uc011ase.1	51205919	51220775	NEGATIVE	Exon 1	51220615	51220775	-779225
chr22	52000000	uc003bnl.1	51205919	51222087	NEGATIVE	Exon 1	51221928	51222087	-777913
chr22	52000000	uc003bns.3	51222156	51238065	POSITIVE	Exon 3	51237082	51238065	-761935
chr22	52000000	uc003bnq.1	51222224	51227600	POSITIVE	Exon 4	51227322	51227600	-772400
chr22	52000000	uc003bnr.1	51222224	51227781	POSITIVE	Exon 4	51227319	51227781	-772219
chr22	52000000	uc010hbj.3	51222224	51238065	POSITIVE	Exon 3	51237082	51238065	-761935
chr22	0	uc002zks.4	16150259	16193004	NEGATIVE	Exon 8	16150259	16151821	16150259
chr22	0	uc002zkt.3	16162065	16172265	POSITIVE	Exon 1	16162065	16162388	16162065
chr22	0	uc002zku.3	16179617	16181004	NEGATIVE	Exon 1	16179617	16181004	16179617
chr22	0	uc002zkv.3	16187164	16193004	NEGATIVE	Exon 5	16187164	16187302	16187164	
```

END_DOC
 */
@Program(name="biostar81455",
	biostars=81455,
	keywords={"bed","gene","knownGene","ucsc"},
	description="Defining precisely the exonic genomic context based on a position ."
	)
public class Biostar81455 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar81455.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-1"},description="The coordinate are one-based. The default is zero based.")
	private  boolean one_based=false;
	@Parameter(names={"-KG","--knownGene"},description=KnownGene.OPT_KNOWNGENE_DESC,required=true)
	private String kgUri = KnownGene.getDefaultUri();
	
	private IntervalTreeMap<List<KnownGene>> kgMap=null;
	
    private static int distance(int pos0,KnownGene kg)
		{
		if(pos0<kg.getTxStart()) return kg.getTxStart()-pos0;
		if(pos0>kg.getTxEnd())  return kg.getTxEnd()-pos0;
		return 0;
		}
    
    private static int distance(int pos0,KnownGene.Exon exon)
    	{
    	if(pos0<exon.getStart()) return exon.getStart()-pos0;
    	if(pos0>exon.getEnd())  return exon.getEnd()-pos0;
    	return 0;
    	}
    
    @Override
    public int doWork(final List<String> args) {
		BufferedReader r=null;
		String line;
		PrintStream out=null;
		final Pattern tab=Pattern.compile("[\t]");
		
		
    	try
			{
    		if(StringUtil.isBlank(this.kgUri))
    			{
    			LOG.error("undefined kguri");
    			return -1;
    			}
    		LOG.info("readubf "+kgUri);
    		this.kgMap = KnownGene.loadUriAsIntervalTreeMap(this.kgUri,(kg)->kg.getExonCount()!=0);
			}
    	catch(final Exception err)
    		{
			LOG.error(err);
			return -1;
    		}
    	finally {
			CloserUtil.close(r);
			}
		final ContigNameConverter contigNameConverter = ContigNameConverter.fromIntervalTreeMap(this.kgMap);
		try
    		{
    		r = super.openBufferedReader(oneFileOrNull(args));
			out = openFileOrStdoutAsPrintStream(this.outputFile);
			while((line=r.readLine())!=null)
				{
				if(line.startsWith("#"))
					{
					out.println(line);
					continue;
					}
				boolean found=false;
				final String tokens[]=tab.split(line);
				final int pos0=Integer.parseInt(tokens[1]) - (this.one_based?1:0);
				final String convertCtg = contigNameConverter.apply(tokens[0]);
				if(convertCtg==null) {
					LOG.error("CANNOT FIND contig "+tokens[0]);
					out.println("##UNKNOWN CONTIG "+line);
					continue;
					}
				final IntervalTree<List<KnownGene>> kgs;
			    if(StringUtil.isBlank(convertCtg))
			    	{
			    	kgs = null;
			    	}
			    else
				    {
					kgs = kgMap.debugGetTree(tokens[0]);
				    }
			    
			    if(kgs==null)
			    	{
			    	LOG.info("no gene found in chromosome "+tokens[0]+" (check chrom prefix?)");
			    	}
			    else
					{
			    	KnownGene bestGene=null;
					
					for(final Iterator<IntervalTree.Node<List<KnownGene>>> iter=kgs.iterator();
						iter.hasNext();)
						{
						for(final KnownGene kg:iter.next().getValue())
							{
							if(bestGene==null|| Math.abs(distance(pos0,kg))< Math.abs(distance(pos0,bestGene)))
								{
								bestGene=kg;
								}
							}
						}
					if(bestGene!=null)
						{
						//get all transcritpts
						final List<KnownGene> overlapKg =new ArrayList<>();
						for(final List<KnownGene> lkg:kgMap.getOverlapping(new Interval(
								bestGene.getChromosome(),
								bestGene.getTxStart()+1,
								bestGene.getTxEnd()
								)))
								{
								overlapKg.addAll(lkg);
								}
						
						for(final KnownGene kg:overlapKg)
							{
							KnownGene.Exon bestExon=null;
							for(KnownGene.Exon exon:kg.getExons())
								{
								if(bestExon==null || Math.abs(distance(pos0, exon))< Math.abs(distance(pos0,bestExon)))
									{
									bestExon=exon;
									}
								}
							if(bestExon!=null)
								{
								out.println(
										line+"\t"+
										kg.getName()+"\t"+
										kg.getTxStart()+"\t"+kg.getTxEnd()+"\t"+kg.getStrand()+"\t"+
										bestExon.getName()+"\t"+bestExon.getStart()+"\t"+bestExon.getEnd()+"\t"+
										distance(pos0,bestExon)
										);
								found=true;
								}
							}
						}
					
					}
				if(!found)
					{
					out.println(line+"\tNULL");
					}
				}
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
    	finally
    		{
    		CloserUtil.close(r);
    		CloserUtil.close(out);
    		}
		}
	
	public static void main(final String[] args)throws Exception
		{
		new Biostar81455().instanceMainWithExit(args);
		}
	}
