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
package com.github.lindenb.jvarkit.tools.translategff3;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.tools.kg2fa.KnownGeneToFasta;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

/**

BEGIN_DOC

Experimental.
Translates the output of bcftools consensus.
Output is a BED file containing the cDNA, the translated peptide and the mean codon usage .

## Example
```
samtools faidx ref.fa 8:11870-11890 |\
	bcftools consensus in.vcf.gz  -s SAMPLE -H A --include 'TYPE="snp"' |\
	java -jar dist/jvarkit.jar translategff3 --gff3 annot.gff3
```
END_DOC

*/
@Program(name="translategff3",
description="translates the output of bcftools consensus",
keywords={"gff","gff3","fasta","genpred","peptide","protein"},
creationDate="20241003",
modificationDate="20241003",
jvarkit_amalgamion = true
)
public class TranslateGff3 extends Launcher {
	private static final Logger LOG = Logger.of(KnownGeneToFasta.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--gff","--gff3"},description="Gff3 file",required=true)
	private Path gff3path = null;

	/** codon usage for Hs https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606 */
	@SuppressWarnings("serial")
	private final Map<String,Float> codonUsage=new HashMap<String,Float>() {{{
		put("UUU",17.6f);  put("UCU",15.2f);  put("UAU",12.2f);  put("UGU",10.6f);
		put("UUC",20.3f);  put("UCC",17.7f);  put("UAC",15.3f);  put("UGC",12.6f);
		put("UUA",7.7f);   put("UCA",12.2f);  put("UAA",1.0f);   put("UGA",1.6f);
		put("UUG",12.9f);  put("UCG",4.4f);   put("UAG",0.8f);   put("UGG",13.2f);
		put("CUU",13.2f);  put("CCU",17.5f);  put("CAU",10.9f);  put("CGU",4.5f);
		put("CUC",19.6f);  put("CCC",19.8f);  put("CAC",15.1f);  put("CGC",10.4f);
		put("CUA",7.2f);   put("CCA",16.9f);  put("CAA",12.3f);  put("CGA",6.2f);
		put("CUG",39.6f);  put("CCG",6.9f);   put("CAG",34.2f);  put("CGG",11.4f);
		put("AUU",16.0f);  put("ACU",13.1f);  put("AAU",17.0f);  put("AGU",12.1f);
		put("AUC",20.8f);  put("ACC",18.9f);  put("AAC",19.1f);  put("AGC",19.5f);
		put("AUA",7.5f);   put("ACA",15.1f);  put("AAA",24.4f);  put("AGA",12.2f);
		put("AUG",22.0f);  put("ACG",6.1f);   put("AAG",31.9f);  put("AGG",12.0f);
		put("GUU",11.0f);  put("GCU",18.4f);  put("GAU",21.8f);  put("GGU",10.8f);
		put("GUC",14.5f);  put("GCC",27.7f);  put("GAC",25.1f);  put("GGC",22.2f);
		put("GUA",7.1f);   put("GCA",15.8f);  put("GAA",29.0f);  put("GGA",16.5f);
		put("GUG",28.1f);  put("GCG",7.4f);   put("GAG",39.6f);  put("GGG",16.5f);
		}}};
	
	private void translateTranscript(PrintWriter w, byte[] sequence,Locatable interval,GeneticCode geneticCode, Gff3Feature gene,Gff3Feature transcript)
		{
		//LOG.info("translate "+transcript.getID()+" "+interval);
		if(!interval.contains(transcript)) return;
		final List<? extends Locatable> cdss=transcript.getChildren().stream().
			filter(G->G.getType().equalsIgnoreCase("cds")).
			map(CDS->new SimpleInterval(CDS)).
			sorted((A,B)->Integer.compare(A.getStart(),B.getStart())).
			collect(Collectors.toList());
		if(cdss.isEmpty()) {
			//LOG.info("no CDS for "+transcript.getID());
			return ;
			}
		
		final CharSequence cdna;
		
		if(transcript.getStrand()==Strand.NEGATIVE) {
			final StringBuilder sb=new StringBuilder(cdss.stream().mapToInt(S->S.getLengthOnReference()).sum());
			for(int y=cdss.size()-1;y>=0;y--) {
				final Locatable cds = cdss.get(y);
				for(int x=cds.getEnd();x>= cds.getStart();--x) {
					sb.append( AcidNucleics.complement((char)sequence[x-interval.getStart()]));
					}
				}
			cdna=sb;
			}
		else
			{
			final StringBuilder sb=new StringBuilder(cdss.stream().mapToInt(S->S.getLengthOnReference()).sum());
			for(Locatable cds: cdss) {
				for(int x=cds.getStart();x<= cds.getEnd();++x) {
					sb.append((char)sequence[x-interval.getStart()]);
					}
				}
			cdna=sb;
			}
		w.print(transcript.getContig());
		w.print("\t");
		w.print(transcript.getStart()-1);
		w.print("\t");
		w.print(transcript.getEnd());
		w.print("\t");
		w.print(transcript.getID());
		w.print("\t");
		w.print(transcript.getStrand()==Strand.NEGATIVE?'-':'+');
		w.print("\t");
		w.print(transcript.getName());
		w.print("\t");
		w.print(gene.getID());
		w.print("\t");
		w.print(gene.getName());
		w.print("\t");
		w.print(cdss.size());
		w.print("\t");
		w.print(cdss.stream().map(L->String.valueOf(L.getStart()-1)).collect(Collectors.joining(",")));
		w.print("\t");
		w.print(cdss.stream().map(L->String.valueOf(L.getEnd())).collect(Collectors.joining(",")));
		w.print("\t");
		w.print(cdna);
		w.print("\t");
		float sum=0f;
		int count=0;
		for(int x=0;x+2< cdna.length();x+=3) {
			String codon = (""+cdna.charAt(x+0) + cdna.charAt(x+1) + cdna.charAt(x+2)).toUpperCase();
			char aa = geneticCode.translate(codon.charAt(0), codon.charAt(1), codon.charAt(2));
			Float frequency_per_thousand = codonUsage.get(codon);
			if(frequency_per_thousand!=null) {
				sum+=frequency_per_thousand;
				count++;
				}
			w.print(aa);
			}
		w.print("\t");
		if(count==0) {
			w.print(".");
			}
		else
			{
			w.print(sum/count);
			}
		
		w.println();
		}
	
	private void translateGene(PrintWriter w, byte[] sequence,Locatable interval,GeneticCode geneticCode,Gff3Feature gene)
		{
		//LOG.info("translate "+gene.getID());
		for(Gff3Feature transcript:gene.getChildren()) {
			translateTranscript(w,sequence,interval,geneticCode,gene,transcript);
			}
		}

	private void translateGenes(PrintWriter w, byte[] sequence,Locatable interval,GeneticCode geneticCode,Collection<Gff3Feature> genes)
		{
		for(Gff3Feature gene: genes) {
			translateGene(w,sequence,interval,geneticCode,gene);
			}
		}
	
	private void recursive(Gff3Feature feat,final IntervalTreeMap<Gff3Feature> geneMap) {
		if(feat.getType().equals("gene")) {
			geneMap.put(new Interval(feat), feat);
			}
		else
			{
			for(Gff3Feature c:feat.getChildren()) {
				recursive(c,geneMap);
				}
			}
		}

	
	@Override
	public int doWork(List<String> args) {
		
		try {
		final String input = oneFileOrNull(args);
		final IntervalTreeMap<Gff3Feature> geneMap = new IntervalTreeMap<>();
		final Set<String> keepAnnots  = new HashSet<>(Arrays.asList(
				Gff3Constants.ID_ATTRIBUTE_KEY,
				Gff3Constants.PARENT_ATTRIBUTE_KEY,
				Gff3Constants.NAME_ATTRIBUTE_KEY,
				"gene_name","gene_id","transcript_id","gene_biotype"
				));
	    final Gff3Codec codec = new Gff3Codec(DecodeDepth.DEEP,ATT->!keepAnnots.contains(ATT));
	    LineIterator lr = IOUtils.openPathForLineIterator(this.gff3path);
	    codec.readHeader(lr);
        while(!codec.isDone(lr)) {
                final Gff3Feature feat = codec.decode(lr);
                if(feat==null) continue;
                recursive(feat,geneMap);
                }
        codec.close(lr);
        if(geneMap.isEmpty()) {
        	LOG.error("no gene was found in "+this.gff3path);
        	return -1;
        	}
        
	    try(BufferedReader br=super.openBufferedReader(input)) {
	    	Collection<Gff3Feature> genes = Collections.emptyList();
	    	final ByteArrayOutputStream baos = new ByteArrayOutputStream();
	    	try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
	    		Locatable interval=null;
		    	for(;;) {
		    		String line=br.readLine();
		    		
		    		if(line==null || line.startsWith(">")) {
		    			if(baos.size()>0 && ! genes.isEmpty()) {
		    				translateGenes(pw, baos.toByteArray(),interval,GeneticCode.getStandard(),genes);
		    				}
		    			
		    			if(line==null) break;
		    			baos.reset();
		    			final String seqName=line.substring(1).trim();
		    			final int hyphen = seqName.lastIndexOf('-');
		    			final int colon = seqName.lastIndexOf(':');
		    			if(hyphen>0 && colon < hyphen) {
		    				interval = new SimpleInterval(
		    						seqName.substring(0,colon),
		    						Integer.parseInt(seqName.substring(colon+1,hyphen)),
		    						Integer.parseInt(seqName.substring(hyphen+1))
		    						);
		    				if(interval.getStart()> interval.getEnd()) throw new IllegalArgumentException("bad interval "+line);
		    				}
		    			else
		    				{
		    				interval = new SimpleInterval(seqName,1,(int) Math.pow(2, 29)/* max tabix */);
		    				}
		    			genes = geneMap.getOverlapping(interval);
		    			if(genes.isEmpty()) {
		    				//LOG.warn("no gene for "+interval);
		    				}
		    			}
		    		else if(!genes.isEmpty())
		    			{
		    			for(int i=0;i< line.length();i++) {
		    				baos.write(line.charAt(i));
		    				}
		    			}
		    		}
	    		pw.flush();
	    		}
	    	}
	    
	    return 0;
		} catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new TranslateGff3().instanceMainWithExit(args);
	}

}
