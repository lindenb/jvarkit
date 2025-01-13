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
package com.github.lindenb.jvarkit.tools.gff2fa;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.LocatableUtils;
import com.github.lindenb.jvarkit.tools.kg2fa.KnownGeneToFasta;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Constants;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

/**

BEGIN_DOC

convert GFF3 to fasta/xml/bed

## Example
```
$ java -jar dist/jvarkit.jar gff2fasta -R ref.fasta in.gff3 --output-format xml | xmllint --format -
(...)
  <entry>
    <attributes>
      <attribute key="gene_name">SCN5A</attribute>
      <attribute key="gff_type">transcript</attribute>
      <attribute key="start">38674799</attribute>
      <attribute key="length">138</attribute>
      <attribute key="type">UTR5</attribute>
      <attribute key="chrom">chr3</attribute>
      <attribute key="transcript_id">ENST00000327956.6</attribute>
      <attribute key="strand">negative</attribute>
      <attribute key="protein_id">ENSP00000333674.6</attribute>
      <attribute key="build">GRCh37</attribute>
      <attribute key="end">38687267</attribute>
      <attribute key="location">chr3:38674799-38687267</attribute>
      <attribute key="transcript_type">protein_coding</attribute>
      <attribute key="gene_id">ENSG00000183873.11</attribute>
    </attributes>
    <sequence>AGTGGACACTGTGGGCATGCGTGTCCGTGCAGCACATCGCCATGCAGGAGCCTGGGGGAAGGCCTTTCCCTCAGTGAGGGCTGCAGCTTCCCCACAGGCAACGTGAGGAGAGCCTGTGCCCAGAAGCAGGATGAGAAG</sequence>
  </entry>
</gff2fasta>

chrom  start     end       strand    location                gene_id             gene_name  gff_type    length  transcript_id  transcript_name  tran
chr3   38589547  38674840  negative  chr3:38589548-38674840  ENSG00000183873.11  SCN5A      transcript  8303    ENST                            
chr3   38591811  38674798  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  5997    ENST                            
chr3   38591811  38674798  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  1999    ENST                            
chr3   38589547  38591811  negative  chr3:38589548-38591811  ENSG00000183873.11  SCN5A      transcript  2264    ENST                            
chr3   38674798  38674840  negative  chr3:38674799-38674840  ENSG00000183873.11  SCN5A      transcript  42      ENST00                          
chr3   38589552  38674853  negative  chr3:38589553-38674853  ENSG00000183873.11  SCN5A      transcript  8362    ENST                            
chr3   38591811  38674798  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  6048    ENST                            
chr3   38591811  38674798  negative  chr3:38591812-38674798  ENSG00000183873.11  SCN5A      transcript  2016    ENST                            
chr3   38589552  38591811  negative  chr3:38589553-38591811  ENSG00000183873.11  SCN5A      transcript  2259    ENST                            

$ java -jar dist/jvarkit.jar gff2fasta -R  in.fasta in.gff3 --hide-mrna --hide-peptide --hide-cdna

(...)
>build:GRCh37|chrom:chr3|end:38674807|gene_id:ENSG00000183873.11|gene_name:SCN5A|gff_type:transcript|length:9|location:chr3:38674799-38674807|protein_id:ENSP00000413996.2|start:38674799|strand:negative|transcript_id:ENST00000449557.2|transcript_type:protein_coding|type:UTR5
GATGAGAAG
>build:GRCh37|chrom:chr3|end:38687267|gene_id:ENSG00000183873.11|gene_name:SCN5A|gff_type:transcript|length:138|location:chr3:38674799-38687267|protein_id:ENSP00000333674.6|start:38674799|strand:negative|transcript_id:ENST00000327956.6|transcript_type:protein_coding|type:UTR5
AGTGGACACTGTGGGCATGCGTGTCCGTGCAGCACATCGCCATGCAGGAGCCTGGGGGAA
GGCCTTTCCCTCAGTGAGGGCTGCAGCTTCCCCACAGGCAACGTGAGGAGAGCCTGTGCC
CAGAAGCAGGATGAGAAG

```
END_DOC

*/
@Program(name="gff2fasta",
description="extract fasta from gtf",
keywords={"gff","gff3","fasta"},
creationDate="20241016",
modificationDate="20241017",
jvarkit_amalgamion = true
)
public class Gff3ToFasta extends Launcher {
	private static final Logger LOG = Logger.build(KnownGeneToFasta.class).make();
	private enum output_format {fasta,xml,bed};
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referencePath = null;
	@Parameter(names={"--hide-mrna"},description="Hide mRNA")
	private boolean hide_mrna=false;
	@Parameter(names={"--hide-cdna"},description="Hide cDNA")
	private boolean hide_cdna=false;
	@Parameter(names={"--hide-peptide"},description="Hide peptide")
	private boolean hide_peptide=false;
	@Parameter(names={"--hide-utr5"},description="Hide UTR5")
	private boolean hide_utr5=false;
	@Parameter(names={"--hide-utr3"},description="Hide UTR3")
	private boolean hide_utr3=false;
	@Parameter(names={"--region","--interval","--regions"},description="Limit to gene overlaping that region chr-end")
	private String regionStr=null;
	@Parameter(names={"--fasta-length"},description="length of fasta lines")
	private int fasta_line_length=60;
	@Parameter(names={"--output-format"},description="output format")
	private output_format outFormat = output_format.fasta;


	
	private static String getProperty(Gff3Feature feat,String key) {
		List<String> L=feat.getAttribute(key);
		return L==null||L.isEmpty()?null:L.get(0);
		}
	
	private abstract class ResultWriter implements Closeable {
		final PrintWriter w;
		ResultWriter(PrintWriter w) {
			this.w= w;
			}
		abstract  void write(FunctionalMap<String, String> props,CharSequence seq) ;
		@Override
		public void close() throws IOException {
			this.w.flush();
			this.w.close();
			}
		}
	private class FastaWriter extends ResultWriter {
		FastaWriter(PrintWriter w) {
			super(w);
			}
	
		 @Override
		 void write(FunctionalMap<String, String> props,CharSequence seq) {
			 w.print(">");
			w.print(props.entrySet().stream().
					sorted((A,B)->A.getKey().compareTo(B.getKey())).
					map(KV->KV.getKey()+":"+KV.getValue()).collect(Collectors.joining("|")));
			for(int i=0;i< seq.length();i++) {
				if(i%Gff3ToFasta.this.fasta_line_length==0) w.println();
				w.append(seq.charAt(i));
				}
			w.println();
		    }
		}

	private class TsvWriter extends ResultWriter {
		final String[] cols=new String[] {"chrom","start","end","strand","location","gene_id","gene_name","gff_type","length",
				"transcript_id","transcript_name","transcript_type","type",
				"protein_id"};
		TsvWriter(PrintWriter w) {
			super(w);
			w.print(String.join("\t", cols));
			w.println("\tsequence");
			}
	
		 @Override
		 void write(FunctionalMap<String, String> props,CharSequence seq) {
			 for(String col:cols) {
				 if(col.equals("start") && props.containsKey(col)) {
					//convert to bed
					  w.append(String.valueOf(Integer.parseInt(props.get(col))-1));
				 	}
				 else
				 {
			      w.append(props.containsKey(col)?props.get(col):".");
				 }
				w.append("\t");	
			 	}
			 w.append(seq);
			 w.append("\n");
			 }
		}
	
	private class XmlWriter extends ResultWriter {
		XMLStreamWriter xw;
		XmlWriter(PrintWriter w) {
			super(w);
			XMLOutputFactory xof = XMLOutputFactory.newFactory();
			try {
				this.xw=xof.createXMLStreamWriter(w);
				this.xw.writeStartDocument("1.0");
				this.xw.writeStartElement("gff2fasta");
				this.xw.writeAttribute("version", getVersion());
				this.xw.writeCharacters("\n");
				}
			catch(XMLStreamException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close() throws IOException {
			try {
				this.xw.writeEndElement();
				this.xw.writeEndDocument();
				this.xw.close();
				}
			catch(XMLStreamException err) {
				throw new IOException(err);
				}
			super.close();
			}
	
		 @Override
		 void write(FunctionalMap<String, String> props,CharSequence seq) {
			 try {
				 this.xw.writeStartElement("entry");
				 this.xw.writeStartElement("attributes");
				 for(String col:props.keySet()) {
					 this.xw.writeStartElement("attribute");
					 this.xw.writeAttribute("key", col);
					 this.xw.writeCharacters(props.get(col));
					 this.xw.writeEndElement();
				 	}
				 this.xw.writeEndElement();
				 this.xw.writeStartElement("sequence");
				 this.xw.writeCharacters(seq.toString());
				 this.xw.writeEndElement();
				 this.xw.writeEndElement();
				 this.xw.writeCharacters("\n");
				}
			catch(XMLStreamException err) {
				throw new RuntimeIOException(err);
				}
			 }
		}
	
	private void echo(ResultWriter w, boolean reverse,FunctionalMap<String, String> props,CharSequence seq) {
		props=props.plus("length",String.valueOf(seq.length()));
		if(props.containsKey("chrom") && props.containsKey("start") && props.containsKey("end")) {
			props=props.plus("location",props.get("chrom")+":"+props.get("start")+"-"+props.get("end"));
			}
		if(reverse) seq=AcidNucleics.reverseComplement(seq);

		w.write(props, seq);
		}
	
	
	
	private void recursive_transcript(
			final Gff3Feature gene,
			final Gff3Feature transcript,
			final ResultWriter pw,
			final ReferenceSequenceFile ref,
			FunctionalMap<String, String> props
			) {
		boolean reverse= transcript.getStrand().equals(Strand.NEGATIVE);
		props = props.plus("strand",transcript.getStrand().name().toLowerCase());
		String s = getProperty(transcript,"transcript_id");
		if(!StringUtils.isBlank(s)) props= props.plus("transcript_id",s);
		s = getProperty(transcript,"transcript_type");
		if(!StringUtils.isBlank(s)) props= props.plus("transcript_type",s);
		s = getProperty(transcript,"transcript_name");
		if(!StringUtils.isBlank(s)) props= props.plus("transcript_name",s);
		s = getProperty(transcript,"protein_id");
		if(!StringUtils.isBlank(s)) props= props.plus("protein_id",s);

		
		props= props.plus("gff_type",transcript.getType());

		
		final List<Gff3Feature> exons = new ArrayList<>();
		final List<Gff3Feature> cdss = new ArrayList<>();
		for(Gff3Feature c:transcript.getChildren()) {
			if(c.getType().equals("exon")) exons.add(c);
			else if(c.getType().equals("CDS")) cdss.add(c);
			}
		Collections.sort(exons,LocatableUtils::compareTo);
		Collections.sort(cdss,LocatableUtils::compareTo);
		
		for(int i=0;i+1< exons.size();i++) {
			final Locatable loc1 = exons.get(i);
			final Locatable loc2 = exons.get(i+1);
			if(!loc1.contigsMatch(loc2)) throw new IllegalStateException("two exons not on the same contig ?? "+loc1+" "+loc2);
			if(loc1.overlaps(loc2)) throw new IllegalStateException("two exons overlap ?? "+loc1+" "+loc2);
			}
		for(int i=0;i+1< cdss.size();i++) {
			final Locatable loc1 = cdss.get(i);
			final Locatable loc2 = cdss.get(i+1);
			if(!loc1.contigsMatch(loc2)) throw new IllegalStateException("two cds not on the same contig ?? "+loc1+" "+loc2);
			if(loc1.overlaps(loc2)) throw new IllegalStateException("two cds overlap ?? "+loc1+" "+loc2);
			}

		
		
		
		if(!this.hide_mrna && !exons.isEmpty()) {
			final StringBuilder mRNA=new StringBuilder(exons.stream().mapToInt(E->E.getLengthOnReference()).sum());

			for(final Gff3Feature exon: exons) {
				mRNA.append(
						ref.getSubsequenceAt(exon.getContig(), exon.getStart(), exon.getEnd()).
						getBaseString());
				}
	
			echo(pw,reverse,
					props.
						plus("type","mRNA").
						plus("start",String.valueOf(exons.get(0).getStart())).
						plus("end",String.valueOf(exons.get(exons.size()-1).getEnd())).
						plus("n-exons",String.valueOf(exons.size())),
					mRNA);
			}
		
		if(!(this.hide_cdna && this.hide_peptide) && !cdss.isEmpty()) {
			final StringBuilder sb=new StringBuilder(cdss.stream().mapToInt(E->E.getLengthOnReference()).sum());
			
			
			for(final Gff3Feature cds: cdss) {
				sb.append(ref.getSubsequenceAt(cds.getContig(), cds.getStart(), cds.getEnd()).getBaseString());
				}
			final CharSequence cdna = (reverse?AcidNucleics.reverseComplement(sb):sb);
			
			if(!this.hide_cdna) {
				echo(pw,false/* already reverse compl*/,
					props.
					plus("type","cDNA").
						plus("start",String.valueOf(cdss.get(0).getStart())).
						plus("end",String.valueOf(cdss.get(cdss.size()-1).getEnd())),
					cdna);
				}
			if(!this.hide_peptide && cdna.length()>2) {
				
				final GeneticCode code=GeneticCode.getStandard();
				final StringBuilder pep=new StringBuilder(cdna.length()/3);
				for(int x=0;x+2< cdna.length();x+=3) {
					pep.append(code.translate(
							cdna.charAt(x+0),
							cdna.charAt(x+1),
							cdna.charAt(x+2))
							);	
					}
				echo(pw,false,props.
						plus("type","peptide").
						plus("start",String.valueOf(cdss.get(0).getStart())).
						plus("end",String.valueOf(cdss.get(cdss.size()-1).getEnd())),
					pep);
				}
			}
		if(!cdss.isEmpty() && !(hide_utr5 && hide_utr3)) {
			
			// [EX   ]===[    EXON  ]======
			// =============[ CDS   ]======
			// [UTR  ]=================
			if(!((hide_utr5 && !reverse) || (hide_utr3 && reverse) )) {
				final Locatable firstCDS = cdss.get(0);
				if(firstCDS.getStart() > exons.get(0).getStart()) {
					final StringBuilder leftutr=new StringBuilder();
					for(Locatable exon:exons) {
						if(exon.getStart() > firstCDS.getEnd()) break;
						leftutr.append(ref.getSubsequenceAt(transcript.getContig(),
								exon.getStart(),
								Math.min(exon.getEnd(),firstCDS.getStart()-1)
								).getBaseString());
						}
					echo(pw,reverse,
						props.
							plus("type",reverse?"UTR3":"UTR5").
							plus("start",String.valueOf(exons.get(0).getStart())).
							plus("end",String.valueOf(firstCDS.getStart()-1)),
						leftutr);
					}
				}

			// ======[    EXON  ]== [EXON]=
			// =========[ CDS ]============
			//======================[ UTR]===
			if(!((hide_utr5 && !reverse) || (hide_utr3 && reverse) )) {
				final Locatable lastCDS = cdss.get(cdss.size()-1);
				if(lastCDS.getEnd()< exons.get(exons.size()-1).getEnd()) {
					final StringBuilder rightUTR=new StringBuilder();
					for(Locatable exon:exons) {
						if(exon.getEnd() <= lastCDS.getEnd()) continue;
						rightUTR.append(ref.getSubsequenceAt(transcript.getContig(),
								Math.max(exon.getStart(),lastCDS.getEnd()+1),
								exon.getEnd()
								).getBaseString());
						}
					echo(pw,reverse,props.
							plus("type",reverse?"UTR5":"UTR3").
							plus("start",String.valueOf(lastCDS.getEnd()+1)).
							plus("end",String.valueOf(exons.get(exons.size()-1).getEnd())),
							rightUTR);
					}
				}
			}
		}
	
	private void recursive(
			final Gff3Feature feat,
			final ResultWriter pw,
			final Locatable region,
			final ReferenceSequenceFile ref
			) {
		
		if(feat.getType().equals("gene") && (region==null || region.overlaps(feat))) {
			final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(ref);
			
			final SAMSequenceRecord ssr= dict.getSequence(feat.getContig());
			if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(feat.getContig(), dict);
			if(feat.getEnd() > ssr.getEnd()) {
				LOG.error("GFF feature "+feat+" is longer than the chromosome (length:"+ssr.getSequenceLength()+")");
				return;
				}
			
			
			FunctionalMap<String,String> props  =FunctionalMap.of("chrom",feat.getContig());
			String s = SequenceDictionaryUtils.getBuildName(dict).orElse(null);
			if(!StringUtils.isBlank(s)) props= props.plus("build",s);

			
			s  = getProperty(feat,"gene_id");
			if(!StringUtils.isBlank(s)) props= props.plus("gene_id",s);
			s = getProperty(feat,"gene_name");
			if(!StringUtils.isBlank(s)) props= props.plus("gene_name",s);
			s = getProperty(feat,"gene_biotype");
			if(!StringUtils.isBlank(s)) props= props.plus("gene_biotype",s);
			
			
			for(Gff3Feature c:feat.getChildren()) {
				recursive_transcript(feat,c,pw,ref,props);
				}
			}
		else
			{
			for(Gff3Feature c:feat.getChildren()) {
				recursive(c,pw,region,ref);
				}
			}
		}

	
	@Override
	public int doWork(List<String> args) {
		
		try {
			this.fasta_line_length=Math.max(1, fasta_line_length);
			
		final String input = oneFileOrNull(args);
		final Set<String> keepAnnots  = new HashSet<>(Arrays.asList(
				Gff3Constants.ID_ATTRIBUTE_KEY,
				Gff3Constants.PARENT_ATTRIBUTE_KEY,
				Gff3Constants.NAME_ATTRIBUTE_KEY,
				"gene_name","gene_id","transcript_id","gene_biotype","protein_id","transcript_type","gene_type"
				));
		
		Locatable limit_interval=
				StringUtils.isBlank(this.regionStr)?
				null:
				LocatableUtils.parse(this.regionStr)
				;
		
		
		try(ReferenceSequenceFile fasta =ReferenceSequenceFileFactory.getReferenceSequenceFile(this.referencePath)) {
			final Gff3Codec codec = new Gff3Codec(DecodeDepth.DEEP,ATT->!keepAnnots.contains(ATT));
		    LineIterator lr = IOUtils.openURIForLineIterator(input);
		    codec.readHeader(lr);
	    	try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
	    		final ResultWriter rw;
	    		switch(this.outFormat) {
	    			case xml: rw = new XmlWriter(pw);break;
	    			case bed: rw = new TsvWriter(pw);break;
	    			default: rw = new FastaWriter(pw);break;
	    			}
		        while(!codec.isDone(lr)) {
		                final Gff3Feature feat = codec.decode(lr);
		                if(feat==null) continue;
		                recursive(feat, rw,limit_interval,fasta);
		                }
		        rw.close();
		    	}
	        codec.close(lr);
			}
	    
	    return 0;
		} catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new Gff3ToFasta().instanceMainWithExit(args);
	}

}
