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
package com.github.lindenb.jvarkit.tools.barcode;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
 * Barcode generator for EricCharp
 *
 */
@Program(
		name="barcodegenerator",
		description="Barcode generator for EricCharp",
		keywords={"barcode"},
		creationDate="20230629",
		modificationDate="20230629",
		jvarkit_amalgamion =  true
		)
public class BarcodeGenerator extends Launcher {
	private static final Logger LOG = Logger.build(BarcodeGenerator.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"--polyx"},description ="max polyx")
	private int max_stretch_inclusive = 3;
	@Parameter(names={"--size"},description ="barcode size")
	private int barcode_length = 8;
	@Parameter(names={"--number"},description ="number of barcode to generated")
	private int n_generate = 96;
	@Parameter(names={"--random"},description ="random seed. -1 == current time")
	private long random_seed=-1;
	@Parameter(names={"--min-gc"},description ="min gc")
	private double min_gc=0.2;
	@Parameter(names={"--max-gc"},description ="max gc")
	private double max_gc=0.8;
	@Parameter(names={"--max-generation"},description ="max generations. -1 : forever")
	private long max_generation=-1L;
	@Parameter(names={"--min-differences"},description ="min differences between barcodes")
	private int min_differences = 4;
	
	private Random random;

	
	private class Solution {
		final List<Barcode> barcodes = new ArrayList<>(BarcodeGenerator.this.n_generate);
		long generation = 0L;
		int getCountCommonBases() {
			int n=0;
			for(int i=0;i+1 < barcodes.size();i++) {
				final Barcode barcodei = barcodes.get(i);
				for(int j=i+1;j < barcodes.size();j++) {
					final Barcode barcodej = barcodes.get(j);
					n+= barcodei.getCountCommonBases(barcodej);
				}
			}
			return n;
		}
		int getMaxConsecutiveBases() {
			int n=0;
			for(int i=0;i+1 < barcodes.size();i++) {
				final Barcode barcodei = barcodes.get(i);
				for(int j=i+1;j < barcodes.size();j++) {
					final Barcode barcodej = barcodes.get(j);
					n= Math.max(n,barcodei.getMaxConsecutiveBases(barcodej));
				}
			}
			return n;
		}
		
		void fill() {
			while(this.barcodes.size()< BarcodeGenerator.this.n_generate) {
				final Barcode barcode = makeNewIndex();
				int i;
				for(i=0;i< this.barcodes.size();i++) {
					final Barcode barcodei=this.barcodes.get(i);
					final int consecutive = barcodei.getMaxConsecutiveBases(barcode);
					if(consecutive> 4) {
						//System.err.println("Too many consecutive bases:\n " +barcode+"\n "+barcodei+"\n "+consecutive);
						break;
						}
					final int common = barcodei.getCountCommonBases(barcode);
					final int differences = barcode.length()-common;
					if(differences < BarcodeGenerator.this.min_differences) 
						{
						break;
						}
					/*
					final double fract =  common /(double)barcode.length();
					if( fract > BarcodeGenerator.this.max_fraction_common_bases ) {
						break;
						}*/
					}
				if(i!=this.barcodes.size()) continue;
				this.barcodes.add(barcode);
				//System.err.println(solution.barcodes.size());
				}
			}
		
		Solution mute() {
			final Solution sol = new Solution();
			sol.barcodes.addAll(this.barcodes);
			sol.barcodes.remove(random.nextInt(sol.barcodes.size()));
			sol.fill();
			return sol;
		}
		
		void sort() {
			Collections.sort(barcodes,(A,B)->A.toString().compareTo(B.toString()));
		}
		
		boolean isBetterThan(final Solution other)  {
			int n1 = getMaxConsecutiveBases();
			int n2 = other.getMaxConsecutiveBases();
			int n =  Integer.compare(n1, n2);
			if(n!=0) return n<0;
			
			n1 = getCountCommonBases();
			n2 = other.getCountCommonBases();
			n =  Integer.compare(n1, n2);
			if(n!=0) return n<0;
			return false;
			}
		void print(PrintWriter out) {
			this.sort();
			out.println(">> generation "+ generation +" common:"+ this.getCountCommonBases()+" max-consecutive:"+this.getMaxConsecutiveBases());
			for(int i=0;i< this.barcodes.size();i++) {
				final Barcode barcode = this.barcodes.get(i);
				out.print(i+1);
				out.print("\t");
				out.print(barcode.bases);
				out.print("\t");
				out.print(barcode.getGCPercent());
				out.print("\t");
				out.print(barcode.getPolyX());
				out.println();
				}
			out.println();
			out.println("position\tA\tT\tG\tC");
			for(int k=0;k< barcode_length;++k) {
				final double[] count = new double[] {0,0,0,0};
				for(int i=0;i< this.barcodes.size();i++) {
					final Barcode barcode = this.barcodes.get(i);
					switch(barcode.charAt(k)) {
						case 'A': count[0]++;break;
						case 'T': count[1]++;break;
						case 'G': count[2]++;break;
						case 'C': count[3]++;break;
						default:break;
						}
					}
				out.print(k+1);
				for(int i=0;i< count.length;i++) {
					out.print("\t");
					out.print((int)((count[i]/barcodes.size())*100.0));
					out.print("%");
					}
				out.println();
				}
			out.println("<< "+ generation);
			out.flush();
			
		}
	}
	
	private class Barcode extends AbstractCharSequence {
		final char[] bases = new char[BarcodeGenerator.this.barcode_length];
		@Override
		public char charAt(int index) {
			return bases[index];
			}
		public double getGCPercent() {
			double gc=0.0;
			for(int i=0;i< length();i++) {
				switch(charAt(i)) {
					case 'G':case 'C': gc++; break;
					default: break;
					}
				}
			return gc/length();
			}
		
		public int getPolyX() {
			int m=0;
			int i=0;
			while(i< length()) {
				int n=1;
				while(i+n < length() && charAt(i)==charAt(i+n)) {
					n++;
					}
				i+=n;
				m = Math.max(n,m);
				}
			return m;
			}
		int getCountCommonBases(final Barcode other) {
			int n=0;
			for(int i=0;i< length();i++) {
				if(charAt(i)==other.charAt(i)) n++;
				}
			return n;
			}
		
		int getMaxConsecutiveBases(final Barcode other) {
			int m=0;
			int i=0;
			while(i< length()) {
				int n=0;
				while(i+n < length() && this.charAt(i+n)==other.charAt(i+n)) {
					n++;
					}
				i+=n+1;
				m = Math.max(n,m);
				}
			return m;
			}
		
		@Override
		public int length() {
			return bases.length;
			}
		}
	private Barcode makeNewIndex() {
		final String ATGC="ATGC";
		for(;;) {
			final Barcode barcode = new Barcode();
			for(int i=0;i< this.barcode_length;++i) {
				barcode.bases[i] = ATGC.charAt(this.random.nextInt(ATGC.length()));
				}
			if(barcode.getPolyX()>= this.max_stretch_inclusive) continue;
			final double gc = barcode.getGCPercent();
			if(gc < this.min_gc || gc > this.max_gc) continue;
			return barcode;
			}
		}
	
	private Solution makeSolution() {
		final Solution solution = new Solution();
		solution.fill();
		return solution;
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			this.random = new Random(random_seed==-1?System.currentTimeMillis():this.random_seed);
			Solution best = null;
			long generation= 0L;
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				for(;;) {
					generation++;
					final Solution sol = makeSolution();
					
					if(best==null || sol.isBetterThan(best)) {
						best=sol;
						sol.generation= generation;
						sol.print(pw);
						}
					
					for(;;)  {
						final Solution mute = best.mute();
						if(best==null || mute.isBetterThan(best)) {
							LOG.info("mute is cool");
							best=mute;
							sol.generation= generation;
							sol.print(pw);
							}
						else
							{
							break;
							}
						}
					
					if(max_generation>0L && max_generation==generation) break;
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			err.printStackTrace();
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new BarcodeGenerator().instanceMain(args);
	}
	
}
