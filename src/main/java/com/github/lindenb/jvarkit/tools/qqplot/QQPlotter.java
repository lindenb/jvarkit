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
package com.github.lindenb.jvarkit.tools.qqplot;

import java.awt.Color;
import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.canvas.Canvas;
import com.github.lindenb.jvarkit.canvas.CanvasFactory;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.MathUtils;
import com.github.lindenb.jvarkit.util.FunctionalMap;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC



# Example


END_DOC
**/
@Program(
		name="qqplotter",
		description="plot QQplot",
		creationDate="20250324",
		modificationDate="20250324",
		keywords= {"qqplot","gwas","statistics"},
		jvarkit_amalgamion = true
		)
public class QQPlotter extends Launcher {
	private static final Logger LOG = Logger.of(QQPlotter.class);
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-d","--delim"},description= "column separator as a java regexe")
	private String separator = "[ \t]+";
	@Parameter(names={"-w","-width"},description= "image width")
	private int image_width = 1000;
	@Parameter(names={"-p","--column-pvalue"},description= "column label for 'p-value'")
	private String pvalue_label = "P";
	@Parameter(names={"--column-label"},description= "column label for 'name of the point'")
	private String label_label = "";
	@Parameter(names={"--column-url"},description= "column label for 'url of the point'")
	private String url_label = "";
	@Parameter(names={"--disable-log10"},description= "data are already -log10(x), so do not apply -log10")
	private boolean disable_minus_log10 = false;

	
	private class Record {
		double p_value;
		String label = null;
		String url = null;
		double radius = -1;
		Color fill = Color.ORANGE;
		Color stroke = Color.YELLOW;
		}
	
	
	@Override
	public int doWork(List<String> args) {
		try {
			final Pattern ws = Pattern.compile(this.separator);
			final Function<String,List<String>> splitter= S->Arrays.asList(ws.split(S));
			final String input=oneFileOrNull(args);
			final List<Record> records =new ArrayList<>();
			try(BufferedReader br = super.openBufferedReader(input)) {
				String line = br.readLine();
				if(line==null) {
					LOG.error("Cannot read headerof "+(input==null?"stdin":input));
					}
				final FileHeader header=new FileHeader(line, splitter);
				final int pvalue_column = header.getColumnIndex(this.pvalue_label);
				final int label_column =  StringUtils.isBlank(this.label_label)?-1:header.getColumnIndex(this.label_label);
				final int url_column =  StringUtils.isBlank(this.url_label)?-1:header.getColumnIndex(this.url_label);
				for(;;) {
					line = br.readLine();
					if(line==null) break;
					final FileHeader.RowMap row = header.toMap(line);
					final Record rec = new Record();
					rec.p_value  = Double.parseDouble(row.at(pvalue_column));
					rec.radius = 4;
					if(!disable_minus_log10) {
						rec.p_value = - Math.log10(rec.p_value);
						}
					if(label_column!=-1) {
						rec.label = row.at(label_column);
						}
					if(url_column!=-1) {
						rec.url = row.at(url_column);
						}
					records.add(rec);
					}
				}
			if(records.isEmpty()) {
				LOG.info("no data");
				return -1;
				}
			Collections.sort(records,(A,B)->Double.compare(A.p_value, B.p_value));
			final double min_y = records.stream().mapToDouble(R->R.p_value).min().orElse(0.0);
			final double max_y = records.stream().mapToDouble(R->R.p_value).max().orElse(1.0);
			
			
			final double[] x_values = MathUtils.ppoints(records.size());
			for(int i=0;i< x_values.length;++i) {
				x_values[i] = -Math.log10(x_values[i]);
				}
			Arrays.sort(x_values);
			final double min_x =x_values[0];
			final double max_x =x_values[x_values.length-1];
			
			try(Canvas canvas = new CanvasFactory().
					setWidth(image_width).
					setHeight(image_width).
					setFormat(CanvasFactory.formatFromPath(outputFile).orElse(CanvasFactory.Format.SVG)).open(
							this.outputFile,
							FunctionalMap.make())) {
				canvas.rect(0, 0, image_width, image_width, FunctionalMap.of(
						Canvas.KEY_FILL,Color.LIGHT_GRAY,
						Canvas.KEY_STROKE,null
						));
				
				canvas.line(0, image_width, image_width, 0, FunctionalMap.of(
						Canvas.KEY_STROKE,Color.DARK_GRAY
						));
				
				for(int i=0;i< records.size();++i) {
					final Record rec = records.get(i);
					double x = ((rec.p_value-min_y)/(max_y-min_y))*image_width;
					double y = image_width - ((x_values[i]-min_x)/(max_x-min_x))*image_width;
					canvas.circle(x, y, rec.radius ,FunctionalMap.of(
							Canvas.KEY_STROKE, rec.stroke ,
							Canvas.KEY_FILL,rec.fill,
							Canvas.KEY_TITLE, rec.label,
							Canvas.KEY_HREF, rec.url
							));
					}
				canvas.rect(0, 0, image_width-1, image_width-1, FunctionalMap.of(
						Canvas.KEY_STROKE,Color.DARK_GRAY,
						Canvas.KEY_FILL,null
						));
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new QQPlotter().instanceMainWithExit(args);
		}
	}
