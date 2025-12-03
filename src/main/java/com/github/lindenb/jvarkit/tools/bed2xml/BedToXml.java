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
package com.github.lindenb.jvarkit.tools.bed2xml;

import java.io.BufferedReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.OptionalLong;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineCodec;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.dict.DictionaryXmlSerializer;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.SmartComparator;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.MinMaxDouble;
import com.github.lindenb.jvarkit.math.MinMaxInteger;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CoordMath;

/**
BEGIN_DOC

## Motivation

convert BED to XML. The idea is to help generating a format that can be processed with XSLT


## Example

```
$ java -jar dist/jvarkit.jar bed2xml src/test/resources/toy.bed.gz -R src/test/resources/toy.dict -d 10  | xmllint --format -

<?xml version="1.0" encoding="UTF-8"?>
<bed src="src/test/resources/toy.bed.gz">
  <header>
    <dictionary md5="df8200dd2a49e25bc98df5f2c45ac36a" length="85" count="2">
      <sequence name="ref" length="45" index="0" offset="0" f1="0.0" f2="0.5294117647058824" M5="7a66cae8ab14aef8d635bc80649e730b" UR="file:/home/lindenb/src/jvarkit-git/src/tests/resources/toy.fa"/>
      <sequence name="ref2" length="40" index="1" offset="45" f1="0.5294117647058824" f2="1.0" M5="1636753510ec27476fdd109a6684680e" UR="file:/home/lindenb/src/jvarkit-git/src/tests/resources/toy.fa"/>
    </dictionary>
  </header>
  <body>
    <contig name="ref" count="2" start="1" end="13" rows="1">
      <row y="0" count="2">
        <rec contig="ref" start="1" end="2" length="1" y="0" f1="0.011764705882352941" f2="0.023529411764705882">
          <col4>Hello</col4>
          <col5>World0</col5>
        </rec>
        <rec contig="ref" start="10" end="13" length="3" y="0" f1="0.11764705882352941" f2="0.15294117647058825">
          <col4>Hello</col4>
          <col5>World1</col5>
        </rec>
      </row>
    </contig>
    <contig name="ref2" count="3" start="1" end="17" rows="2">
      <row y="0" count="2">
        <rec contig="ref2" start="1" end="2" length="1" y="0" f1="0.5411764705882353" f2="0.5529411764705883">
          <col4>Hello</col4>
          <col5>World2</col5>
        </rec>
        <rec contig="ref2" start="13" end="14" length="1" y="0" f1="0.6823529411764706" f2="0.6941176470588235">
          <col4>Hello</col4>
          <col5>World3</col5>
        </rec>
      </row>
      <row y="1" count="1">
        <rec contig="ref2" start="16" end="17" length="1" y="1" f1="0.7176470588235294" f2="0.7294117647058823">
          <col4>Hello</col4>
          <col5>World4</col5>
        </rec>
      </row>
    </contig>
  </body>
  <footer>
    <contig name="ref" start="1" end="13"/>
    <contig name="ref2" start="1" end="17"/>
  </footer>
</bed>

```


END_DOC
*/

@Program(
	name="bed2xml",
	description="Convert BED to XML",
	keywords={"bed","xml"},
	creationDate = "20251128",
	modificationDate = "20251128",
	jvarkit_amalgamion = true
	)
public class BedToXml extends Launcher {
	private enum BedType {generic,bedGraph,bed12};
	private static final Logger LOG=Logger.of(BedToXml.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--dict-out"},description="Save reduced/modified dict to this file")
	private Path outputDictFile = null;
	@Parameter(names={"--type"},description="BedType. bedGraph will print the min/max of value. bed12 will display each block ")
	private BedType bedType= BedType.generic;
	@Parameter(names={"--omit-xml-desclaration"},description="Don't print XM declaration.")
	private boolean omit_xml_decl = false;

	@Parameter(names={"--regex"},description="keep chromosomes matching that regular expression")
	private String contig_regex=null;
	@Parameter(names={"--columns"},description="Optional Comma separated list names of column starting from the 3rd column (after 'end'). Use '.' to ignore.")
	private String column_names_str="";
	@Parameter(names={"--min-contig-length"},description="keep chromosomes which length is greater than 'x'")
	private int min_contig_length=0;
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE)
	private Path faidPath=null;
	@Parameter(names={"-d","--distance"},description="if >=0, add a 'y' attribute that could be used to display the bed records in a browser, " 
					+ "	this 'y' is the graphical row where the item should be displayed. " 
					+ "	This distance is the distance between two item where there is a collision. Memory consuming . "+
					DistanceParser.OPT_DESCRIPTION,
					converter = DistanceParser.StringConverter.class, splitter = NoSplitter.class)
	private int distance=-1;

	private void writeBedRecord(
		final XMLStreamWriter w,
		final BedLine bed,
		final String normContig,
		final OptionalInt optY,
		final Map<String,Long> contig2idx,
		final OptionalLong genome_length,
		final List<String> column_names
		) throws XMLStreamException {
	w.writeStartElement("rec");
	w.writeAttribute("contig",normContig);
	w.writeAttribute("start", String.valueOf(bed.getStart()-1));
	w.writeAttribute("end", String.valueOf(bed.getEnd()));
	w.writeAttribute("length", String.valueOf(bed.getLengthOnReference()));
	if(optY.isPresent()) 	w.writeAttribute("y", String.valueOf(optY.getAsInt()));
	if(contig2idx!=null && genome_length.isPresent()) {
		final long idx = contig2idx.get(normContig);
		w.writeAttribute("f1", String.valueOf((idx+bed.getStart()-1)/(double)genome_length.getAsLong()));
		w.writeAttribute("f2", String.valueOf((idx+bed.getEnd())/(double)genome_length.getAsLong()));
		}
	
	for(int col=3;col < bed.getColumnCount();col++) {
		if(col>=9 && col< 12 && this.bedType.equals(BedType.bed12)) continue;
		final String v = bed.getOrDefault(col, "");
		
		final String tag;
		if(col>=3  && col-3 < column_names.size() && !StringUtils.isBlank(column_names.get(col-3)) && !column_names.get(col-3).equals(".")) {
			tag  = column_names.get(col-3);
			}
		else if(this.bedType.equals(BedType.bed12) && col < 9) {
			switch(col) {
				case 3: tag = "name";break;
				case 4: tag = "score";break;
				case 5: tag = "strand";break;
				case 6: tag = "thickStart";break;
				case 7: tag = "thickEnd";break;
				case 8: tag = "itemRgb";break;
				default: throw new IllegalStateException();
				}
			}
		else if(this.bedType.equals(BedType.bedGraph) && col ==3) {
			tag ="score";
			}
		else
			{
			tag = "col"+(col+1);
			}
		
		w.writeStartElement(tag);
		if(contig2idx!=null && genome_length.isPresent() && this.bedType.equals(BedType.bed12) && (col==6 || col==7)) {
			final long idx = contig2idx.get(normContig);
			final int thick_pos  = Integer.parseInt(bed.getOrDefault(col, "thick missing"));
			w.writeAttribute("f", String.valueOf((idx+thick_pos)/(double)genome_length.getAsLong()));
			}
		
		w.writeCharacters(v);
		w.writeEndElement();
		}
	if(this.bedType.equals(BedType.bed12) && bed.getColumnCount()>11) {
		int block_count = Integer.parseInt(bed.get(9));
		final int[] sizes = CharSplitter.COMMA.splitAsStringList(bed.get(10)).stream().mapToInt(S->Integer.parseInt(S)).toArray();
		final int[] starts = CharSplitter.COMMA.splitAsStringList(bed.get(11)).stream().mapToInt(S->Integer.parseInt(S)).toArray();
		w.writeStartElement("blocks");
		w.writeAttribute("count", String.valueOf(block_count));
		for(int i=0;i<block_count;i++) {
			final int block_start = bed.getStart() - 1 +starts[i];
			final int block_end   = starts[i]+sizes[i];
			w.writeEmptyElement("block");
			w.writeAttribute("start", String.valueOf(block_start));
			w.writeAttribute("end", String.valueOf(block_end));
			
			if(contig2idx!=null && genome_length.isPresent()) {
				final long idx = contig2idx.get(normContig);
				w.writeAttribute("f1", String.valueOf((idx+block_start)/(double)genome_length.getAsLong()));
				w.writeAttribute("f2", String.valueOf((idx+block_end)/(double)genome_length.getAsLong()));
				}
			}
		w.writeEndElement();
		}
	w.writeEndElement();
	w.writeCharacters("\n");
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Pattern contigRegex = (this.contig_regex==null?null:Pattern.compile(this.contig_regex));
			final MinMaxDouble minMaxBedGraph = new MinMaxDouble();
			final MinMaxDouble minMaxBedScore = new MinMaxDouble();
			final String input=oneFileOrNull(args);
			final ContigNameConverter  ctgConvert;
			final SAMSequenceDictionary dict;
			final Map<String,Long> contig2idx;
			final Map<String,MinMaxInteger> contig2minmaxpos=new HashMap<>();
			final List<String> colNames = CharSplitter.COMMA.splitAsStringList(this.column_names_str);
			final OptionalLong genome_length;
			int max_rows = -1;
			if(faidPath!=null) {
				final SAMSequenceDictionary dict0 = new SequenceDictionaryExtractor().extractRequiredDictionary(faidPath);
				
				dict = new SAMSequenceDictionary(
						dict0.getSequences().stream().
						filter(SSR->min_contig_length<=0 || SSR.getSequenceLength()>=min_contig_length).
						filter(SSR->contigRegex==null?true:contigRegex.matcher(SSR.getSequenceName()).matches()).
						collect(Collectors.toList())
						);
				if(dict.isEmpty()) {
					LOG.warn("empty dictionary");
					}
				ctgConvert = ContigNameConverter.fromOneDictionary(dict);
				contig2idx = new HashMap<>(dict.size());
				long x=0;
				for(SAMSequenceRecord ssr:dict.getSequences()) {
					contig2idx.put(ssr.getContig(), x);
					x+=ssr.getLengthOnReference();
					}
				genome_length = OptionalLong.of(dict.getReferenceLength());
				
				if(this.outputDictFile!=null && !dict.isEmpty()) {
					if(Files.isSameFile(this.outputDictFile, this.faidPath)) {
						LOG.warn("dict-out path is the same as dict-in");
						return -1;
						}
					else
						{
						 try(PrintWriter writer = IOUtils.openPathForPrintWriter(this.outputDictFile)) {
							final SAMSequenceDictionaryCodec codec = new SAMSequenceDictionaryCodec(writer);
							codec.encode(dict);
							writer.flush();
							}
						}
					}
				
				}
			else
				{
				dict=null;
				ctgConvert = ContigNameConverter.getIdentity();
				contig2idx = null;
				genome_length = OptionalLong.empty();
				}
			final BedLineCodec bedCodec = new BedLineCodec();
			try(BufferedReader br=(input==null||input.equals("-")?IOUtils.openStreamForBufferedReader(stdin()):IOUtils.openURIForBufferedReading(input))) {
				try(OutputStream out= super.openPathOrStdoutAsStream(this.outputFile)) {
					final Map<String,List<BedLine>> contig2rec = (this.distance<0?null:new HashMap<>());
					final XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
					final String encoding = "UTF-8";
					final XMLStreamWriter w= xmlfactory.createXMLStreamWriter(out, encoding);
					if(!omit_xml_decl) w.writeStartDocument(encoding, "1.0");
					w.writeStartElement("bed");
					if(input!=null) w.writeAttribute("src", input);
					w.writeAttribute("type", bedType.name());
					if(distance>0) w.writeAttribute("distance", String.valueOf(distance));
					
					w.writeComment("Generated with bed2xml . Jvarkit by Pierre Lindenbaum PhD. Version:"+JVarkitVersion.getInstance().getVersion());
					
					w.writeStartElement("header");
					if(dict!=null) {
						new DictionaryXmlSerializer().writeDictionary(w, dict);
						}
					w.writeEndElement();//header
					w.writeStartElement("body");
					for(;;) {
						String line = br.readLine();
						if(line==null) {
							if(contig2rec!=null) {
								final List<String> ordered_ctg= new ArrayList<>(contig2rec.keySet());
								if(dict!=null) {
									Collections.sort(ordered_ctg,new ContigDictComparator(dict));
									}
								else
									{
									Collections.sort(ordered_ctg,new SmartComparator());
									}
								for(String ctg: ordered_ctg) {
									final List<BedLine> records = contig2rec.get(ctg);
									if(records==null || records.isEmpty()) continue;
									w.writeStartElement("contig");
									w.writeAttribute("name", ctg);
									w.writeAttribute("count", String.valueOf(records.size()));
									w.writeAttribute("start", String.valueOf(records.stream().mapToInt(R->R.getStart()-1).min().getAsInt()));
									w.writeAttribute("end", String.valueOf(records.stream().mapToInt(R->R.getEnd()).max().getAsInt()));
									
									final Pileup<BedLine> cmp = new Pileup<>((LEFT,RIGHT)->CoordMath.getLength(LEFT.getEnd(),RIGHT.getStart())>= BedToXml.this.distance);
									cmp.addAll(records);
									final List<List<BedLine>> rows = cmp.getRows();
									max_rows = Math.max(rows.size(), max_rows);
									w.writeAttribute("rows", String.valueOf(rows.size()));
									
									for(int y=0;y< rows.size();++y) {
										w.writeStartElement("row");
										w.writeAttribute("y", String.valueOf(y));
										w.writeAttribute("count", String.valueOf( rows.get(y).size()));
										for(BedLine rec: rows.get(y)) {
											writeBedRecord(w,rec,ctg,OptionalInt.of(y),contig2idx,genome_length,colNames);
											}
										w.writeEndElement();//contig
										}
									w.writeEndElement();//contig
									}
								
								}
							break;
							}
						final BedLine rec = bedCodec.decode(line);
						if(rec==null) continue;
						final String ctg = ctgConvert.apply(rec.getContig());
						if(StringUtils.isBlank(ctg)) continue;
						
						if(dict==null && contigRegex!=null) {
							if(!contigRegex.matcher(ctg).matches()) continue; 
							}
						
						MinMaxInteger coordMinMax = contig2minmaxpos.get(ctg);
						if(coordMinMax==null) {
							coordMinMax=new MinMaxInteger();
							contig2minmaxpos.put(ctg, coordMinMax);
							}
						coordMinMax.accept(rec.getStart()-1);
						coordMinMax.accept(rec.getEnd());
						
						
						if(bedType.equals(BedType.bedGraph)) {
							minMaxBedGraph.accept(Double.parseDouble(rec.getOrDefault(3, "bedgraph: column $3 missing")));
							}
						else if(bedType.equals(BedType.bed12) && rec.getColumnCount()>4 &&  StringUtils.isInteger(rec.getOrDefault(4, ""))) {
							minMaxBedScore.accept(Integer.parseInt(rec.get(4)));
							}
						
						if(contig2rec!=null) {
							List<BedLine> records = contig2rec.get(ctg);
							if(records==null) {
								records = new ArrayList<>();
								contig2rec.put(ctg, records);
								}
							records.add(rec);
							}
						else
							{
							writeBedRecord(w,rec,ctg,OptionalInt.empty(),contig2idx,genome_length,colNames);
							}
						}
					
					
					w.writeEndElement();//body
					
					w.writeStartElement("footer");
					if(!minMaxBedGraph.isEmpty()) {
						w.writeEmptyElement("bedgraph");
						w.writeAttribute("min", String.valueOf(minMaxBedGraph.getMinAsDouble()));
						w.writeAttribute("max", String.valueOf(minMaxBedGraph.getMaxAsDouble()));
						}
					if(!minMaxBedScore.isEmpty()) {
						w.writeEmptyElement("score");
						w.writeAttribute("min", String.valueOf(minMaxBedScore.getMinAsDouble()));
						w.writeAttribute("max", String.valueOf(minMaxBedScore.getMaxAsDouble()));
						
						}
					for(String ctg: contig2minmaxpos.keySet()) {
						final MinMaxInteger mM = contig2minmaxpos.get(ctg);
						if(mM.isEmpty()) continue;
						w.writeEmptyElement("contig");
						w.writeAttribute("name", ctg);
						w.writeAttribute("start", String.valueOf(mM.getMinAsInt()));
						w.writeAttribute("end", String.valueOf(mM.getMaxAsInt()));
						}
					
					if(max_rows>0) {
						w.writeEmptyElement("max-rows");
						w.writeAttribute("count", String.valueOf(max_rows));
						}
					
					w.writeEndElement();
					
					w.writeEndElement();//bed
					w.writeEndDocument();
					w.close();
					}
				}
			return 0;
		} catch (final Throwable e) {
			LOG.error(e);
			return -1;
		}
		
		
		
		}
	
	public static void main(final String[] args) {
		new BedToXml().instanceMainWithExit(args);
	}

}
