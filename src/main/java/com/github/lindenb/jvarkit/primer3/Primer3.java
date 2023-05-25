/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;

public class Primer3 extends AbstractList<Primer3.Pair> {
private String SEQUENCE_ID = null;
private String SEQUENCE_TEMPLATE = null;
private int PRIMER_FIRST_BASE_INDEX=0;
private List<Primer> primerLeft =  new ArrayList<>();
private List<Primer> primerRight =  new ArrayList<>();
private List<Pair> pairs =  new ArrayList<>();

private abstract class Entity {
	protected final int index;
	protected final Map<String,String> attributes=new HashMap<>();
	protected Entity(int index) {
		this.index= index;
	}
	public int getIndex() {
		return this.index;
	}
	
	protected String getAttribute(String key) {
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
	
}

public class Primer extends Entity implements CharSequence, Locatable {
	private String cached_seq = null;
	private final int side;// 0 LEFT 1 RIGHT
	Primer(int index,int side) {
		super(index);
		this.side=side;
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
	public boolean equals(Object obj) {
		if(obj==this)return true;
		if(obj==null || !(obj instanceof Primer)) return false;
		return getSequence().equals(Primer.class.cast(obj).getSequence());
		}
	@Override
	public int hashCode() {
		return getSequence().hashCode();
		}	
	protected String getEntityName() {
		return (side==0?"LEFT":"RIGHT");
		}
	
	@Override
	public String getContig() {
		return Primer3.this.getContig();
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
	public String toString() {
		return getSequence();
		}
	}


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

private static class MyIter extends AbstractCloseableIterator<Primer3> {
	private final BufferedReader br;
	MyIter(final Reader r) throws IOException {
		this.br = new BufferedReader(r);
		}
	@Override
	protected Primer3 advance() {
		try {
			Primer3 curr = null;
			String line;
			while((line=br.readLine())!=null) {
				if(line.equals("=")) {
					if(curr!=null) return curr;
					continue;
					}
				if(curr==null) curr = new Primer3();
				final int eq = line.indexOf("=");
				if(eq==-1) throw new IOException("Cannot find '=' in primer3 output "+line);
				curr.put(line.substring(0,eq),line.substring(eq+1));
				}
			return curr;
			}
		catch(IOException err) {
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

