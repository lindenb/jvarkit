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

package com.github.lindenb.jvarkit.tools.gloussouarn;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.logging.Logger;
import java.util.prefs.Preferences;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;


/**
 * 2019-12-10 outil demande par G Loussouarn pour manipuler des fichiers tabulaires
 * 
 * le fichier de confg consitste en des paires d'interval (time-start, time-end)
 * le fichier de donnees contient 2 valeurs: (time: value)
 * 
 * on veut un tableau croise avec les intervals
 * 
 * @author lindenb
 *
 */
public class GLoussouarn01  {
	private static final Logger LOG=Logger.getLogger(GLoussouarn01.class.getName());
	private static final String SUFFIX=".convert.txt";
	private static final Pattern SPACES = Pattern.compile("[ \t]+");
	private static class XY{
		double time;
		double y;
		}

	private class Converter {
	private final List<XY> data =  new ArrayList<>(1_000_000);
	private final List<UserInterval> intervals =  new ArrayList<>(100);
	private Double deltaTime = null;

	
	
	private  class UserInterval {
		double timeStart;
		double timeEnd;
		List<XY> data= null;
		
		String getTitle() { return String.valueOf(timeStart)+":"+this.timeEnd;}
		int getRows() { return this.data.size();}
		String getAt(int y) {
			if(y>=this.data.size()) return ".";
			return String.format("%.4f",this.data.get(y).y);
			}
	}
	
		
	private BufferedReader open(final Path path) throws IOException {
		if(path.getFileName().endsWith(".gz")) {
			return new BufferedReader(new InputStreamReader(Files.newInputStream(path)));
		} else
			{
			return Files.newBufferedReader(path);
			}
		}
	
	private int readUserIntervals(final Path path) throws IOException {
		try(BufferedReader br=open(path)) {
		for(;;) {
			final String line= br.readLine();
			if(line==null) break;
			if(line.trim().isEmpty()) continue;
			final String tokens[]=SPACES.split(line);
			if(tokens.length!=2) throw new IOException("Expected two words in "+line);
			final UserInterval ui = new UserInterval();
			ui.timeStart = Double.parseDouble(tokens[0].trim());
			ui.timeEnd = Double.parseDouble(tokens[1].trim());
			if(ui.timeStart>ui.timeEnd) throw new IOException("Bad times in line "+line);
			
			ui.data = this.data.stream().
				filter(xy-> xy.time>=ui.timeStart && xy.time<=ui.timeEnd).
				sorted((A,B)->Double.compare(A.time, B.time)).
				collect(Collectors.toList()); 
				
			
			this.intervals.add(ui);
			}
		} 
		return this.intervals.size();
	}
	
	
	private int readData(final Path path) throws IOException {
		final double precision = 1E-3;
		try(BufferedReader br=open(path)) {
		for(;;) {
			final String line= br.readLine();
			if(line==null) break;
			if(line.trim().isEmpty()) continue;
			final String tokens[]=SPACES.split(line);
			if(tokens.length!=2) throw new IOException("Expected two words in "+line);
			final XY xy = new XY();
			xy.time = Double.parseDouble(tokens[0].trim());
			xy.y = Double.parseDouble(tokens[1].trim());
			if(!data.isEmpty()) {
				final XY prev=data.get(data.size()-1);
				final double dtime = xy.time - prev.time;
				if(dtime<=0) throw new IOException("Bad order by time");
				if(this.deltaTime==null ) {
					this.deltaTime = dtime;
					LOG.info("delta time: "+this.deltaTime);
				} else if(Math.abs(this.deltaTime.doubleValue()- dtime) > precision)
					throw new IOException("Bad delta times !! got "+this.deltaTime+" and then "+dtime+" line "+line);
			}
			this.data.add(xy);
			}
		} 
		return this.data.size();
		}
	Path convert(final Path intervalPath,final Path inPath) throws IOException {
		if(inPath.getFileName().toString().endsWith(SUFFIX)) throw new IOException("File "+inPath+" ends with "+SUFFIX);
		
		final SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMddHHmmss");
		
		final Path out = inPath.getParent().resolve(sdf.format(new Date())+"."+inPath.getFileName().toString()+SUFFIX);
		if(Files.exists(out)) throw new IOException("File exists "+out);
		if( this.readData(inPath)==0) {
			throw new IOException("no data in "+inPath);
			}
		
		if( this.readUserIntervals(intervalPath)==0) {
			throw new IOException("no data in "+intervalPath);
			}
		
		try(PrintWriter w= new PrintWriter(Files.newBufferedWriter(out))) {
			w.print("#t");
			for(final UserInterval ui: this.intervals) {
				w.print("\t");
				w.print(ui.getTitle());
			}
			w.println();
			
			final int maxrows = this.intervals.stream().mapToInt(S->S.getRows()).max().orElse(0);
			for(int y=0;y< maxrows;y++) {
				w.printf("%.6f",y*this.deltaTime);
				for(final UserInterval ui: this.intervals) {
					w.print("\t");
					w.print(ui.getAt(y));
				}
				w.println();
			}
			
			
		w.flush();
		}
		return out;
	}
	
	}
	
	
	
	
	 int instanceMain(final String args[]) {
		try {
			Locale.setDefault(Locale.ENGLISH);
			
			final Preferences prefs = Preferences.userNodeForPackage(GLoussouarn01.class);

			
			if(args.length!=0) {
				LOG.severe("too many args");
				return -1;
			}
			
			JFileChooser jfc = new JFileChooser();
			File directory = null;
			String dirStr =prefs.get("config.dir", null);
			if(dirStr!=null) directory=new File(dirStr);
			if(directory==null || !directory.exists() || !directory.isDirectory()) directory=null;
			
			jfc.setCurrentDirectory(directory);
			jfc.setDialogTitle("Configuration file");
			
			if(jfc.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION) {
				LOG.warning("User canceled");
				return -1;
			}
			final File configFile = jfc.getSelectedFile();
			if(configFile==null || !configFile.isFile()) {
				LOG.warning("User canceled");
				return -1;
			}
			prefs.put("config.dir",configFile.getParentFile().getPath());
			
			jfc = new JFileChooser();
			
			directory=null;
			dirStr =prefs.get("data.dir", null);
			if(dirStr!=null) directory=new File(dirStr);
			if(directory==null || !directory.exists() || !directory.isDirectory()) directory=null;

			
			jfc.setDialogTitle("Data files");
			jfc.setCurrentDirectory(directory);
			jfc.setMultiSelectionEnabled(true);
			jfc.setFileFilter(new FileFilter() {
				@Override
				public String getDescription() {
					return "Data files";
				}
				
				@Override
				public boolean accept(File f) {
					return f!=null && !f.getName().endsWith(SUFFIX);
				}
			});
			
			
			if(jfc.showOpenDialog(null)!=JFileChooser.APPROVE_OPTION) {
				LOG.warning("User canceled");
				return -1;
			}
			final File dataFiles[] = jfc.getSelectedFiles();
			if(dataFiles==null || dataFiles.length==0) {
				LOG.warning("User canceled");
				return -1;
				}
			
			int i=0;
			for(final File dataFile:dataFiles) {
				final Converter converter = new Converter();
				Path out = converter.convert(configFile.toPath(), dataFile.toPath());
				LOG.info("converted : "+out);
				i++;
				
				prefs.put("data.dir",dataFile.getParentFile().getPath());

			}
			JOptionPane.showMessageDialog(null, String.valueOf(i)+ " file(s) converted");
			
			prefs.sync();
			return 0;
		} catch(final Throwable err ) {
			err.printStackTrace();
			JOptionPane.showMessageDialog(null, String.valueOf(err.getMessage()));
			return -1;
		}
	}
	
	public static void main(final String[] args) {
		try {
		SwingUtilities.invokeAndWait(()->{
			int ret=new GLoussouarn01().instanceMain(args);
			System.exit(ret);
		});
		} catch(final Throwable err) {
			err.printStackTrace();
		}
	}

}
