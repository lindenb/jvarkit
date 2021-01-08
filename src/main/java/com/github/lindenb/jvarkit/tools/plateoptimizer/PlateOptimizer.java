/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.plateoptimizer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.function.Function;
import java.util.function.IntSupplier;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;

/**
BEGIN_DOC

## Example

```

```

END_DOC
 */
@Program(
		name="plateoptimizer",
		creationDate="20191125",
		modificationDate="20200114",
		generate_doc=false,
		keywords= {"plate"}
		)
public class PlateOptimizer extends Launcher {
	private static final Logger LOG = Logger.build(PlateOptimizer.class).make();
	@Parameter(names= {"-o","--ouput"},description="Output prefix",required=true)
	private String output;
	@Parameter(names= {"--seed"},description="Random seed . -1 == use currentTimeMillis")
	private long seed= -1L;

	final int num_female_per_plates = 8;
	private final String watched_columns[]={"Typeind","Type_Donneur","hsf_Patient"};
	private final List<String> headers = new ArrayList<>();
	private final Map<String,Integer> column2index = new HashMap<>();
	private final List<Content> all_contents = new ArrayList<>();
	private volatile int ctrlc_catched_count=5;
	private Random random = new Random(System.currentTimeMillis());
	
	private class Content 
		{
		final int content_index ;
		final String tokens[];
		Content(final int content_index,final String tokens[]) {
			this.content_index = content_index;
			this.tokens=tokens;
			}
		
		
		@Override
		public String toString() {
			return String.join("\t", tokens);
			}
		String get(final String s) {
			final int col = getColumnByName(s);
			return col>=this.tokens.length?"NA":this.tokens[col];
			}
		String getSexe() {
			return get("Sexe");
			}
		String getTube() {
			return get("Type_Tubes");
			}
		String getBox() {
			String s= get("Pos_Box_DIV");
			int hyphen = s.indexOf("-");
			if(hyphen==-1) throw new IllegalStateException(s);
			return s.substring(0,hyphen);
			}
		
		
		String getTypeInd() {
			return get("Typeind");
		}
		String getTypeDonneur() {
			return get("Type_Donneur");
		}
		String getHsfPatient() {
			return get("hsf_Patient");
		}
		
		boolean isFemale() {
			return getSexe().equals("F");
			}
		boolean isMale() {
			return getSexe().equals("M");
			}
		String getName() {
			return get("id_KIT");
			}
		boolean isAlive() {
			return getTypeDonneur().equals("Vivant");
			}
		}
	
	/*
	private class Criteria {
		int column;
		int score=-10;
		Predicate<String> emptyTest = S->StringUtils.isBlank(S) || S.equals("NA") || S.equals("N/A");
		
		BiPredicate<String,String> comparator = (A,B)->{
			if(emptyTest.test(A) && emptyTest.test(B)) return false;
			return A.equals(B);
			};
		boolean isSameContent(final Content c1,final Content c2) {
			String s1 = c1.tokens[this.column];
			String s2 = c2.tokens[this.column];
			return this.comparator.test(s1, s2);
		}
	}*/
	
	private class Cell 
		{
		final Plate plate;
		final int x;
		final int y;
		Content content = null;
		Cell(final Plate plate,int x,int y) {
			this.plate = plate;
			this.x=x;
			this.y=y;
			}
		String getLabel() {
			return (char)((int)'A'+y)+String.format("%02d", x+1);
			}
		/** return true if cell is empty */
		boolean isEmpty() {
			return this.content == null;
		}
		
		private boolean isEmptyString(final String s) {
			return StringUtils.isBlank(s) || s.equals("NA") || s.equals("N/A");
			}
		int getScore(final Function<Content, String> cellExtract) {
			if(this.isEmpty()) return 0;
			int score=0;
			for(int x=-1;x<=1;++x)
				{
				for(int y=-1;y<=1;++y)
					{
					if(x==0 && y==0) continue;
					final Cell neighbour = this.plate.get(this.x + x, this.y + y);
					if( neighbour == null ) continue;
					if( neighbour.isEmpty()) continue;
					final String s1 = cellExtract.apply(neighbour.content);
					if(isEmptyString(s1)) continue;
					final String s2 = cellExtract.apply(this.content);
					if(isEmptyString(s2)) continue;
					if(!s1.equals(s2)) continue;
					
					final int n= 100;
					if(x==0 || y==0)
						{
						score += n;
						}
					else //diagonal
						{
						score += n/10;
						}
					
					}
				}
			return score;
			}
		}
	
	private class Plate implements Iterable<Cell>{
		final String name;
		final Cell cells[] = new Cell[96];
		Plate(final String name) {
			this.name  = name;
			int i=0;
			for(int y=0;y<8;++y)
				for(int x=0;x<12;++x)
					{
					final Cell cell=new Cell(this,x,y);
					this.cells[i++] = cell;
					}
			}
		Plate(final Plate cp) {
			this(cp.name);
			for(int i=0;i< this.cells.length;i++) this.cells[i].content=cp.cells[i].content;
			}
		
		public Stream<Cell> stream() {
			return Arrays.stream(this.cells);
			}
		
		@Override
		public Iterator<Cell> iterator() {
			return stream().iterator();
			}
		
		private OptionalDouble getRatioX(Predicate<Cell> pred1,Predicate<Cell> pred2) {
			final long count1 = this.stream().filter(C->!C.isEmpty()).filter(pred1).count();
			final long count2 = this.stream().filter(C->!C.isEmpty()).filter(pred2).count();
			if(count1+count2==0L) return OptionalDouble.empty();
			return OptionalDouble.of(count1/(double)(count1+count2));
			}
		OptionalDouble getRatioTypeInd() {
			return getRatioX(
					C->C.content.getTypeInd().equals("D"),
					C->C.content.getTypeInd().equals("R")
					);
			}
		
		OptionalDouble getRatioTypeDonneur() {
			final long countD = this.stream().filter(C->!C.isEmpty()).filter(C->C.content.isAlive()).count();
			final long countR = this.stream().filter(C->!C.isEmpty()).filter(C->!C.content.isAlive()).count();
			if(countR+countD==0L) return OptionalDouble.empty();
			return OptionalDouble.of(countD/(double)(countR+countD));
			}
		
		
		void shuffle() {
			final List<Content> contents = this.stream().
					map(C->C.content).
					collect(Collectors.toCollection(ArrayList::new));
			Collections.shuffle(contents,PlateOptimizer.this.random);
			for(int i=0;i< contents.size();i++) {
				this.cells[i].content=contents.get(i);
			}
		}
		
		Cell get(int x,int y) {
			if(x<0 || x>=12) return null;
			if(y<0 || y>=8) return null;
			return this.cells[y*12+x];
			}
		/* number of original plates */
		int getNumberOfTubes() {
			return  stream().
					filter(C->!C.isEmpty()).
					map(C->C.content.getTube()).
					collect(Collectors.toSet()).
					size();
		}
		
		int getNumberOfBoxes() {
			return stream().
					filter(C->!C.isEmpty()).
					map(C->C.content.getBox()).
					collect(Collectors.toSet()).
					size();
		}
		
		boolean hasEmptyCell() {
			return stream().anyMatch(C->C.isEmpty());
		}
		
		void writeXML(final XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("table");
			
			w.writeStartElement("caption");
			w.writeCharacters(this.name +
					" boxes:"+this.getNumberOfBoxes()+
					" ratio D/R :"+this.getRatioTypeDonneur()+
					" ratio Ind :"+this.getRatioTypeInd()+
					" n-females:"+ this.stream().filter(C->!C.isEmpty()).filter(C->C.content.isFemale()).count());
			w.writeEndElement();
			
			
			w.writeStartElement("thead");
			
			w.writeStartElement("tr");
			w.writeEmptyElement("th");
			for(int x=0;x< 12;++x) {
				w.writeStartElement("th");
				w.writeCharacters(String.format("%02d", x+1));
				w.writeEndElement();//th
				}
			w.writeEndElement();//tr
			w.writeEndElement();//thead
			
			w.writeStartElement("tbody");
			for(int y=0;y< 8;++y) {
				w.writeStartElement("tr");
				w.writeStartElement("th");
				w.writeCharacters(String.valueOf( (char)((int)'A'+y)));
				w.writeEndElement();//th
				for(int x=0;x< 12;++x) {
					w.writeStartElement("td");
					final Cell cell = this.get(x, y);
					if(!cell.isEmpty()) {
						String bckg="background:AliceBlue;";
						if(cell.content.isFemale()) bckg+="color:pink;";
						if(cell.content.isMale()) bckg+="color:blue;";
						w.writeAttribute("style", bckg);
						w.writeCharacters(cell.content.getName()+":");
						w.writeCharacters(String.join(";",cell.content.getBox(),cell.content.getTypeDonneur(),cell.content.getTypeInd()));
						}
					w.writeEndElement();//td
					}
				w.writeEndElement();//tr
				}
			
			w.writeEndElement();//tbody
			w.writeEndElement();
			}
	
	}
	
	private class Solution implements Comparable<Solution>{
		final long id;
		final List<Plate> plates = new ArrayList<>();
		int score = 0;
		Solution(final long id) {
			this.id = id;
		}
		
		
		List<Plate> getPlatesButLastOne() {
			if(this.plates.isEmpty()) throw new IllegalStateException();
			return this.plates.subList(0, this.plates.size()-1);
		}
		
		boolean hasEmptyCells() {
			//empty cell only allowed in last plate
			for(int i=0;i+1/* YES*/ < this.plates.size();i++) {
				if(this.plates.get(i).hasEmptyCell()) return true;
			}
			return false;
		}
		
		@Override
		public int compareTo(final Solution o) {
			
			
			int i1 = this.getPlatesButLastOne().
					stream().
					mapToInt(T->T.getNumberOfBoxes()).
					max().orElse(10000);
			int i2= o.getPlatesButLastOne().
					stream().
					mapToInt(T->T.getNumberOfBoxes()).
					max().orElse(10000);
			
			int d= Integer.compare(i1, i2);
			if(d!=0) return d;
			
			final Predicate<OptionalDouble> testOpt = opt->opt.isPresent() && opt.getAsDouble()>=0.4 && opt.getAsDouble() <=0.6;
			
			i1 = (int)this.getPlatesButLastOne().
					stream().
					filter(P->testOpt.test(P.getRatioTypeDonneur())).
					count();
			i2 = (int)o.getPlatesButLastOne().
					stream().
					filter(P->testOpt.test(P.getRatioTypeDonneur())).
					count();
			
			d= Integer.compare(i2, i1);//inverse
			if(d!=0) return d;
			
			
			i1 = (int)this.getPlatesButLastOne().
					stream().
					filter(P->testOpt.test(P.getRatioTypeInd())).
					count();
			i2 = (int)o.getPlatesButLastOne().
					stream().
					filter(P->testOpt.test(P.getRatioTypeInd())).
					count();
			d= Integer.compare(i2,i1);//inverse
			if(d!=0) return d;
			
			return 0;
			}
		
		
		public int compareToOld(final Solution o) {
			final Function<Content, String> xcell1 = C->C.getTypeInd();
			final Function<Content, String> xcell2 = C->C.getTypeDonneur();
			final Function<Content, String> xcell3 = C->C.getHsfPatient();
				
			
			for(final Function<Content, String> xcell: Arrays.asList(xcell1,xcell2,xcell3)) {
				int i1 = this.plates.
						stream().
						flatMap(P->P.stream()).
						mapToInt(C->C.getScore(xcell)).
						sum();
				int i2= o.plates.
						stream().
						flatMap(P->P.stream()).
						mapToInt(C->C.getScore(xcell)).
						sum();
				final int i= Integer.compare(i1, i2);
				if(i!=0) return i;
				}
			int i1 = this.plates.
					stream().
					mapToInt(T->T.getNumberOfTubes()).
					sum();
			int i2= o.plates.
					stream().
					mapToInt(T->T.getNumberOfTubes()).
					sum();
			
			return Integer.compare(i1, i2);
			}
		
		

		private void fillEmptyPlates() {
			final int n_plates = (int)Math.ceil(PlateOptimizer.this.all_contents.size()/96.0);
			if(n_plates==0) throw new IllegalStateException();
			for(int p=0;p<n_plates;++p)
				{
				final Plate plate = new Plate("PLATE_ID"+(p+1));
				this.plates.add(plate);
				}
			}
		
		Solution fill2() {
			this.fillEmptyPlates();
			final List<Content> contents = new ArrayList<>(all_contents);
			final IntSupplier nextFemaleIndex = ()->{
				for(int i=0;i< contents.size();i++)
					{
					if(contents.get(i).isFemale()) return i;
					}
				return -1;
				};
			
			//fill females
			for(final Plate p:this.plates) {
				int index=0;
				while(index < num_female_per_plates && index< p.cells.length) {
 					int i= nextFemaleIndex.getAsInt();
 					if(i==-1) break;
 					Content female = contents.remove(i);
 					p.cells[index++].content = female;
					}
				}
			//fill other cells
			for(final Plate p:this.plates) {
				int index=0;
				while(index< p.cells.length && p.cells[index].content!=null)
					{
					index++;
					}
				while(index< p.cells.length && !contents.isEmpty())
					{
					p.cells[index].content = contents.remove(0);//remove front
					index++;
					}
				p.shuffle();
				}
			if(!contents.isEmpty()) throw new IllegalStateException();
			return this;
			}
		
		
		Solution fill() {
			this.fillEmptyPlates();
			
			final List<Cell> cells = this.plates.stream().
					flatMap(P->P.stream()).
					limit(PlateOptimizer.this.all_contents.size()).//evite puits vides
					collect(Collectors.toCollection(ArrayList::new));
			Collections.shuffle(cells,PlateOptimizer.this.random);
			int i=0;
			for(final Content content: PlateOptimizer.this.all_contents) {
				cells.get(i).content  = content;
				i++;
				}
			
			return this;
			}
		
		private void mute(final List<Cell> cells) {
			if(cells.isEmpty()) return;
			int repeat = 1 +  PlateOptimizer.this.random.nextInt(5);

			int max_fail = 0;
			while(max_fail<10_000 && repeat>0) {
				int index1= PlateOptimizer.this.random.nextInt(cells.size());
				int index2= PlateOptimizer.this.random.nextInt(cells.size());
				if(index1==index2) {
					max_fail++;
					continue;
					}
				final Content c1 = cells.get(index1).content;
				final Content c2 = cells.get(index2).content;
				if(c1==null && c2==null)  {
					max_fail++;
					continue;
					}
				cells.get(index1).content = c2;
				cells.get(index2).content = c1;
				repeat--;
				}
			}
		
		Solution mute(long iter_id)
			{
			final Solution sol = new Solution(iter_id);
			
			sol.plates.addAll(
				this.plates.stream().
				map(P->new Plate(P)).
				collect(Collectors.toList())
				);
			
			boolean mute_per_plate = random.nextBoolean();
			
			if(mute_per_plate) {
				for(final Plate plate:sol.plates) {
					final List<Cell> cells = Arrays.asList(plate.cells);
					mute(cells);
					}
				}
			else
				{
				final List<Cell> cells = sol.plates.stream().
						flatMap(P->Arrays.stream(P.cells)).
						collect(Collectors.toCollection(ArrayList::new));				
				mute(cells);
				}
			return sol;
			}
		
		private void save() {
			try {
				try(PrintWriter w= new PrintWriter(output+".txt")) {
					w.print("#plate\tcell\tsrc.index");
					for(final String header:PlateOptimizer.this.headers) w.print("\t"+header);
					w.println();
					for(final Plate p:this.plates) {
						
						for(final Cell cell:p.cells) {
							w.print(p.name);
							w.print("\t");
							w.print(cell.getLabel());
							w.print("\t");
							w.print(cell.content==null?".":String.valueOf(cell.content.content_index+1));
							w.print("\t");
							if(cell.content!=null) {
								w.print(String.join("\t", cell.content.tokens));
								}
							w.println();
							}
						}
					w.flush();
					}
				final XMLOutputFactory xof = XMLOutputFactory.newFactory();
				try(PrintWriter pw= new PrintWriter(output+".html")) {
					final XMLStreamWriter w=xof.createXMLStreamWriter(pw);
					w.writeStartElement("html");
					
					w.writeStartElement("head");
					w.writeStartElement("title");
					w.writeCharacters("Iteration :"+this.id);
					w.writeEndElement();//title
					w.writeStartElement("style");
					w.writeCharacters("table{font-family:Arial,Verdana,sans-serif;color:darkgray;font-size:14px;border-collapse:collapse;border:1px solid black;padding:5px;margin:5px;}");
					w.writeEndElement();//title
					
					w.writeEndElement();//head
					
					w.writeStartElement("body");
					w.writeStartElement("h2");
					w.writeCharacters("Iteration :"+this.id);
					w.writeEndElement();//h2
					w.writeStartElement("div");
					for(final Plate p:this.plates) {
						p.writeXML(w);
					}
					w.writeEndElement();
					w.writeEmptyElement("hr");
					w.writeEndElement();//body
					w.writeEndElement();//html
					w.flush();
					w.close();
					}

				}
			catch(IOException|XMLStreamException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public String toString() {
			return String.valueOf(this.id);
			}
		}
	
	private Solution best=null;
	
	private void iteration(long iter_id) {
		Solution sol = null;
		if(iter_id == 0 || iter_id%20==0 || this.best==null)
			{
			sol = new Solution(iter_id).fill();
			}
		else
			{
			sol = this.best.mute(iter_id);
			}
		
		if(sol.hasEmptyCells()) return;
		
		if(!sol.plates.stream().
				allMatch(P->P.stream().map(C->!C.isEmpty()  && C.content.isFemale()).count()>=num_female_per_plates)) return;

		if(this.best==null || sol.compareTo(this.best)<0)
			{
			this.best=sol;
			LOG.info("["+iter_id+"]");
			sol.save();
			}
		}
	
	private int getColumnByName(final String s) {
		if(!this.column2index.containsKey(s)) {
			throw new IllegalArgumentException("not column '"+s+"' in "+String.join(" , ", this.column2index.keySet()));
			}
		return this.column2index.get(s);
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.seed!=-1) {
				this.random = new Random(this.seed);
				}
			else
				{
				this.random = new Random(System.currentTimeMillis());
				}
			final CharSplitter tab = CharSplitter.TAB;
			final String input = oneFileOrNull(args);
			try(BufferedReader br = (input==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(input)
					)) {
				String line = br.readLine();
				if(line==null) {
					LOG.error("cannot read first line");
				}
				this.headers.addAll(Arrays.asList(tab.split(line)));
				for(int i=0;i< this.headers.size();++i) this.column2index.put(this.headers.get(i),i);
				LOG.info("header: "+String.join(";", this.headers));
				while((line=br.readLine())!=null) {
					if(StringUtils.isBlank(line)) continue;
					line = line.replaceAll("<e9>", "e");
					final String tokens[]=tab.split(line);
					if(tokens.length!=this.headers.size()) throw new JvarkitException.TokenErrors(this.headers.size(), tokens);
					final Content content = new Content(this.all_contents.size(),tokens);
					this.all_contents.add(content);
				}
			}
			if(this.all_contents.isEmpty()) {
				LOG.error("no data defined");
				return -1;
			}
			
			
			Runtime.getRuntime().addShutdownHook(new Thread() {
			        public void run() {
			        	ctrlc_catched_count --;
			        	if(ctrlc_catched_count==0) {
			                LOG.warn("Shutting down ...");
				        	Runtime.getRuntime().removeShutdownHook(this);
			        		}
			        	else
			        		{
			        		LOG.info("press again "+ctrlc_catched_count+" times.");
			        		}
			        }
			    });
			
			
			long n_iteration = 0L;
			while(this.ctrlc_catched_count>0) {
				iteration(++n_iteration);
				}
			ctrlc_catched_count=0;
			LOG.info("done");
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			
			}
		}
	
	public static void main(final String[] args) {
		new PlateOptimizer().instanceMainWithExit(args);
	}
}
