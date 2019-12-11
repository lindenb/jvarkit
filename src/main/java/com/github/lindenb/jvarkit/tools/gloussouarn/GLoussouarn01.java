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

package com.github.lindenb.jvarkit.tools.gloussouarn;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;


/**
 * 2019-12-10 outil demande par G Loussouarn pour manipuler des fichiers tabulaires
 * @author lindenb
 *
 */
public class GLoussouarn01  {
	private static final Logger LOG=Logger.getLogger(GLoussouarn01.class.getName());
	private final List<XY> data =  new ArrayList<>(1_000_000);

	private static class XY{
		float x;
		float y;
		}
	
	private class Slice
		{
		final int x1;
		final int x2;
		Slice(int x1,int x2) {
			this.x1=x1;
			this.x2=x2;
			}
		String getTitle() { return "Y";}
		int getRows() { return x2-x1;}
		String getAt(int y) {
			final int idx = x1+y;
			if(idx>=x2) return ".";
			return String.valueOf(GLoussouarn01.this.data.get(idx));
			}
		}
	
	private int readData(final Path path) throws IOException {
		BufferedReader br = null;
		try {
		if(path.getFileName().endsWith(".gz")) {
			br = new BufferedReader(new InputStreamReader(Files.newInputStream(path)));
		} else
			{
			br = Files.newBufferedReader(path);
			}
		for(;;) {
			final String line= br.readLine();
			if(line==null) break;
			if(line.trim().isEmpty()) continue;
			final int tab= line.indexOf('\t');
			if(tab==-1) throw new IOException("cannot find tabulation in "+line	);
			final XY xy = new XY();
			xy.x = Float.parseFloat(line.substring(0,tab).trim());
			xy.y = Float.parseFloat(line.substring(tab+1).trim());
			this.data.add(xy);
			}
		br.close();
		} finally
			{
			if(br!=null) br.close();
			}
		return this.data.size();
		}
	private void usage(final PrintStream out) {
		out.println(this.getClass().getSimpleName());
		out.println("Usage: java -jar dist/"+this.getClass().getSimpleName().toLowerCase()+".jar [option] file");
		out.println("Options:");
		out.println("  -h help this screen.");
		out.println("  -o (file) output file (required)");
		out.println();
		}
	private int instanceMain(final String args[]) {
		try {
			Path out = null;
			int optind=0;
			while(optind< args.length) {
				if(args[optind].equals("-h") ) {
					usage(System.out);
					return 0;
					}
				if(args[optind].equals("-o") && optind+1 < args.length) {
					out = Paths.get(args[++optind]);
					}
				else if(args[optind].equals("--")) {
					++optind;
					break;
				}else if(args[optind].startsWith("-")) {
					LOG.severe("Illegal option "+args[optind]);
					return -1;
				}
				else
					{
					break;
					}
				++optind;
			}
			if(optind+1!=args.length) {
				LOG.severe("Illegal number of arguments.");
				usage(System.err);
				return -1;
				}
			if(out==null)  {
				LOG.severe("Output file is undefined.");
				usage(System.err);
				return -1;
				}
			final Path inPath = Paths.get(args[optind]);
			if(inPath.equals(out)) {
				LOG.severe("input and output are same file.");
				return -1;
			}
			if( this.readData(inPath)==0) {
				LOG.severe("no data in "+inPath);
				return -1;
			}
			
			final List<Slice> slices = new  ArrayList<>();
			int x1=0;
			while(x1< this.data.size()) {
				final XY xy1 = this.data.get(x1);
				if(xy1.x<0) { x1++; continue;}
				
				int x2 = x1 +1 ;
				
				while(x2< this.data.size() && this.data.get(x2).x==xy1.x) {
					x2++;
					}
				slices.add(new Slice(x1,x2));
				x1 = x2;
			}
			if(slices.isEmpty()) {
				LOG.severe("no slice was found.");
				return -1;
			}
			
			try(PrintWriter w= new PrintWriter(Files.newBufferedWriter(out))) {
				w.print("#X");
				for(final Slice slice: slices) {
					w.print("\t");
					w.print(slice.getTitle());
				}
				w.println();
				
				final int maxrows = slices.stream().mapToInt(S->S.getRows()).max().orElse(0);
				for(int y=0;y< maxrows;y++) {
					final Slice first =  slices.get(0);
					w.print("?");
					for(final Slice slice: slices) {
						w.print("\t");
						w.print(slice.getAt(y));
					}
					w.println();
				}
				
				
			w.flush();
			}
			return 0;
		} catch(final Throwable err ) {
			return -1;
		}
	}
	
	private void instanceMainWithExit(final String args[]) {
		int err= instanceMain(args);
		System.exit(err);
		}
	public static void main(String[] args) {
		new GLoussouarn01().instanceMainWithExit(args);
	}

}
