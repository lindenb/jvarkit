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
package com.github.lindenb.jvarkit.primer3;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

/** describe a result for primer3 */
public class Primer3 extends AbstractList<Primer3.Pair> {
private String SEQUENCE_ID = "undefined";
private String SEQUENCE_TEMPLATE = null;
private int PRIMER_FIRST_BASE_INDEX=0;
private final List<Primer> primerLeft =  new ArrayList<>();
private final List<Primer> primerRight =  new ArrayList<>();
private final List<Pair> pairs =  new ArrayList<>();

private abstract class Entity implements Locatable {
	protected final int index;
	protected final Map<String,String> attributes=new HashMap<>();
	protected Entity(int index) {
		this.index= index;
	}
	public int getIndex() {
		return this.index;
	}
	
	@Override
	public String getContig() {
		return Primer3.this.getContig();
		}

	
	protected String getAttribute(final String key) {
		if(!this.attributes.containsKey(key)) throw new IllegalArgumentException("Cannot get "+key+" in "+this.attributes);
		return this.attributes.get(key);
		}
	
	void put(final String[] key_split_underscore,final String value) {
		final String key = String.join("_", Arrays.asList(key_split_underscore).subList(3, key_split_underscore.length));
		this.attributes.put(key,value);
		}
	StringBuilder writeBoulder(StringBuilder w) {
		for(final String key:this.attributes.keySet()) {
			w.append("PRIMER_");
			w.append(getEntityName());
			w.append("_");
			w.append(String.valueOf(getIndex()));
			if(!key.isEmpty()) w.append("_");
			w.append(key);
			w.append("=");
			w.append(attributes.get(key));
			w.append("\n");
			}
		return w;
		}
	protected abstract String getEntityName();
	abstract void writeXml(final XMLStreamWriter w) throws XMLStreamException;

	
}
/** describe a primer  */
public class Primer extends Entity implements CharSequence {
	private String cached_seq = null;
	private final int side;// 0 LEFT 1 RIGHT
	Primer(int index,int side) {
		super(index);
		this.side=side;
		}
	public boolean isLeft() {
		return this.side==0;
		}
	
	public boolean isRight() {
		return this.side==1;
	}
	
	public String getSequence() {
		if(cached_seq==null) {
			cached_seq =  getAttribute("SEQUENCE").toUpperCase();
			}
		return cached_seq;
		}
	public double getMeltingTemperature() {
		return Double.parseDouble(getAttribute("TM"));
		}
	public double getGcPercent() {
		return Double.parseDouble(getAttribute("GC_PERCENT"));
		}
	
	public double getSelfAny() {
		return Double.parseDouble(getAttribute("SELF_ANY"));
		}
	public double getSelfEnd() {
		return Double.parseDouble(getAttribute("SELF_END"));
		}
	public double getHairpinTh() {
		return Double.parseDouble(getAttribute("HAIRPIN_TH"));
		}
	@Override
	public char charAt(int index) {
		return getSequence().charAt(index);
		}
	@Override
	public int length() {
		return getSequence().length();
		}
	public Pair getPair() {
		return pairs.get(getIndex());
		}
	@Override
	public CharSequence subSequence(int start, int end) {
		return  getSequence().subSequence(start, end);
		}
	@Override
	public boolean equals(final Object obj) {
		if(obj==this)return true;
		if(obj==null || !(obj instanceof Primer)) return false;
		return getSequence().equals(Primer.class.cast(obj).getSequence());
		}
	@Override
	public int hashCode() {
		return getSequence().hashCode();
		}	
	protected String getEntityName() {
		return (isLeft()?"LEFT":"RIGHT");
		}
	

	@Override
	public int getStart() {
		return Integer.parseInt(CharSplitter.COMMA.split(getAttribute(""))[0])+1;
		}
	
	@Override
	public int getEnd() {
		return getStart() + this.length() -1;
		}
	
	@Override
	void writeXml(final XMLStreamWriter w) throws XMLStreamException {
		w.writeStartElement(getEntityName().toLowerCase()+"-primer");
		w.writeAttribute("index", String.valueOf(getIndex()));
		w.writeAttribute("start", String.valueOf(getStart()));
		w.writeAttribute("end", String.valueOf(getEnd()));
		w.writeAttribute("length", String.valueOf(length()));
		w.writeAttribute("tm", String.valueOf(getMeltingTemperature()));
		w.writeAttribute("gc", String.valueOf(getGcPercent()));
		w.writeAttribute("self-any", String.valueOf(getSelfAny()));
		w.writeAttribute("self-end", String.valueOf(getSelfEnd()));
		w.writeAttribute("hairpin-th", String.valueOf(getHairpinTh()));
		w.writeStartElement("sequence");
		w.writeCharacters(getSequence());
		w.writeEndElement();
		w.writeEndElement();
		}
	
	@Override
	public String toString() {
		return getSequence();
		}
	}

/** describe a pair of primers */
public class Pair extends Entity {
	Pair(int index) {
		super(index);
		}
	public Primer getLeft() {
		return primerLeft.get(getIndex());
	}
	
	public Primer getRight() {
		return primerRight.get(getIndex());
	}
	
	public double getPenalty() { return Double.parseDouble(getAttribute("PENALTY"));}
	public double getComplAnyTh() { return Double.parseDouble(getAttribute("COMPL_ANY_TH"));}
	public double getComplEndth() { return Double.parseDouble(getAttribute("COMPL_END_TH"));}
	public int getProductSize() { return Integer.parseInt(getAttribute("PRODUCT_SIZE"));}

	@Override
	public int hashCode() {
		return getLeft().hashCode()*31 + getRight().hashCode();
		}	

	
	@Override
	public boolean equals(Object obj) {
		if(obj==this)return true;
		if(obj==null || !(obj instanceof Pair)) return false;
		return getLeft().equals(Pair.class.cast(obj).getLeft()) &&
				getRight().equals(Pair.class.cast(obj).getRight());
		}
	@Override
	protected String getEntityName() {
		return "PAIR";
		}
	
	@Override
	void writeXml(final XMLStreamWriter w) throws XMLStreamException {
		w.writeStartElement("pair");
		w.writeAttribute("index", String.valueOf(getIndex()));
		w.writeAttribute("start", String.valueOf(getStart()));
		w.writeAttribute("end", String.valueOf(getEnd()));
		w.writeAttribute("length", String.valueOf(getProductSize()));
		getLeft().writeXml(w);
		getRight().writeXml(w);
		w.writeEndElement();
		}

	
	@Override
	public int getStart() {
		return getLeft().getStart();
		}
	@Override
	public int getEnd() {
		return getRight().getEnd();
		}
	
	@Override
	public String toString() {
		final StringBuilder sb=new StringBuilder();
		writeBoulder(sb);
		getLeft().writeBoulder(sb);
		getRight().writeBoulder(sb);
		return sb.toString();
		}
	}


private final CharSequence templateSequence = new AbstractCharSequence() {
	@Override
	public int length() {
		return PRIMER_FIRST_BASE_INDEX+SEQUENCE_TEMPLATE.length();
	}
	
	@Override
	public char charAt(int index) {
		if(index < PRIMER_FIRST_BASE_INDEX ) return 'N';
		index -= PRIMER_FIRST_BASE_INDEX;
		if(index >= SEQUENCE_TEMPLATE.length()) return 'N';
		return SEQUENCE_TEMPLATE.charAt(index);
	}
};

private Primer3() {
}

public String getContig() {
	return SEQUENCE_ID;
	}

private void put(final String key,final String value) {

	if(key.equals("SEQUENCE_ID")) {
		this.SEQUENCE_ID = value;
		return;
		}
	else if(key.equals("SEQUENCE_TEMPLATE")) {
		this.SEQUENCE_TEMPLATE = value;
		return;
		}
	else if(key.equals("PRIMER_FIRST_BASE_INDEX")) {
		this.PRIMER_FIRST_BASE_INDEX = Integer.parseInt(value);
		return;
		}
	final String[] tokens= CharSplitter.UNDERSCORE.split(key);
	if(tokens[0].equals("PRIMER") && tokens.length>1) {
		int index= (tokens.length>2 && StringUtils.isInteger(tokens[2]) ? Integer.parseInt(tokens[2]):-1);
		if(index < 0 ) return;
		if(tokens[1].equals("PAIR")) {
			while(index >= this.pairs.size()) pairs.add(new Pair(pairs.size()));
			pairs.get(index).put(tokens,value);
			}
		else if(tokens[1].equals("LEFT")) {
			while(index >= this.primerLeft.size()) primerLeft.add(new Primer(primerLeft.size(),0));
			primerLeft.get(index).put(tokens,value);
			}
		else if(tokens[1].equals("RIGHT")) {
			while(index >= this.primerRight.size()) primerRight.add(new Primer(primerRight.size(),1));
			primerRight.get(index).put(tokens,value);
			}
		}
	}

public CharSequence getTemplate() {
	return this.templateSequence;
	}

public Locatable getTemplateInterval() {
	return new SimpleInterval(
			this.SEQUENCE_ID,
			this.PRIMER_FIRST_BASE_INDEX+1,
			this.PRIMER_FIRST_BASE_INDEX+SEQUENCE_TEMPLATE.length()
			);
	}


@Override
public Pair get(int index) {
	return this.pairs.get(index);
	}
@Override
public int size() {
	return this.pairs.size();
	}

public void writeXml(final XMLStreamWriter w) throws XMLStreamException {
	w.writeStartElement("primer3");
	w.writeStartElement("template");
	w.writeAttribute("id", this.SEQUENCE_ID);
	w.writeAttribute("first-base-index", String.valueOf(this.PRIMER_FIRST_BASE_INDEX));
	w.writeCharacters(this.SEQUENCE_TEMPLATE);
	w.writeEndElement();
	w.writeStartElement("pairs");
	w.writeAttribute("count", String.valueOf(size()));
	for(Pair p: this.pairs) {
		p.writeXml(w);
		}
	w.writeEndElement();
	w.writeEndElement();
	}

private Primer3 validate() {
	if(StringUtils.isBlank(SEQUENCE_ID)) throw new IllegalArgumentException("SEQUENCE_ID MISSING");
	if(StringUtils.isBlank(SEQUENCE_TEMPLATE)) throw new IllegalArgumentException("SEQUENCE_TEMPLATE MISSING");
	return this;
	}

@Override
public String toString() {
	final StringBuilder sb=new StringBuilder();
	sb.append("SEQUENCE_ID=").append(this.SEQUENCE_ID).append("\n");
	sb.append("SEQUENCE_TEMPLATE=").append(this.SEQUENCE_TEMPLATE).append("\n");
	sb.append("PRIMER_FIRST_BASE_INDEX=").append(this.PRIMER_FIRST_BASE_INDEX).append("\n");
	for(Pair p :this.pairs) {
		sb.append(p.toString());
		}
	sb.append("=\n");
	return sb.toString();
	}

public static CloseableIterator<Primer3> iterator(final Reader r) throws IOException {
	return new MyIter(r);
	}

public static CloseableIterator<Primer3> iterator(final Path path) throws IOException {
	return iterator(IOUtils.openPathForBufferedReading(path));
	}

private static class MyIter extends AbstractCloseableIterator<Primer3> {
	private final BufferedReader br;
	MyIter(final Reader r) throws IOException {
		this.br = (r instanceof BufferedReader? BufferedReader.class.cast(r):new BufferedReader(r));
		}
	@Override
	protected Primer3 advance() {
		try {
			Primer3 curr = null;
			String line;
			while((line=br.readLine())!=null) {
				if(line.equals("=")) {
					if(curr!=null) return curr.validate();
					continue;
					}
				if(curr==null) curr = new Primer3();
				final int eq = line.indexOf("=");
				if(eq==-1) throw new IOException("Cannot find '=' in primer3 output "+line);
				curr.put(line.substring(0,eq),line.substring(eq+1));
				}
			return curr==null?null:curr.validate();
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	@Override
	public void close() {
		try {br.close();}
		catch(IOException err) {}
		}
	}

}

