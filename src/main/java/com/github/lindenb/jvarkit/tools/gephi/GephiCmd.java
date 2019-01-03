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
package com.github.lindenb.jvarkit.tools.gephi;

import java.awt.Color;
import java.io.*;
import java.net.URI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Properties;
import java.util.TreeMap;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.gephi.project.api.*;
import org.openide.util.*;
import org.gephi.io.importer.api.*;
import org.gephi.io.importer.spi.*;
import org.gephi.io.processor.plugin.*;
import org.gephi.graph.api.*;
import org.gephi.layout.plugin.*;
import org.gephi.layout.plugin.forceAtlas2.*;
import org.gephi.io.exporter.preview.*;
import org.gephi.io.exporter.plugin.*;
import org.gephi.preview.api.*;
import org.gephi.preview.types.*;
import org.gephi.layout.spi.*;
import org.gephi.io.exporter.api.*;
import com.itextpdf.text.*;
import org.gephi.utils.progress.*;
import org.gephi.layout.plugin.fruchterman.*;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.StringUtil;

import org.gephi.preview.types.*;

import java.awt.Font;
import java.awt.Rectangle;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;

// https://github.com/gephi/gephi-toolkit-demos/blob/master/src/main/java/org/gephi/toolkit/demos/ImportExport.java
// https://github.com/Beziks/webgephi/blob/9897d1a0cb0f3c353c0e316d93a4111786beb41d/webgephi-server/src/main/webapp/WEB-INF/classes/cz/cokrtvac/webgephi/webgephiserver/core/gephi/GephiImporter.java
/** 
BEGIN_DOC

## Compilation

tested with gephi-0.9.2

a variable gephi_home must be defined in local.mk pointing to the local installation of gephi.

```
$ cat local.mk

gephi_home=${HOME}/package/gephi-0.9.2
```


## Input

input is a GEXF file or it reads a GEXF from stdin.



## Example

list the available properties:

```
 java -jar dist/gephicmd.jar -l
```


convert gexf file to SVG

```
 
 $ java -jar dist/gephicmd.jar -e 'layout.duration:10;layout.time.unit:SECONDS;layout.algorithm:fruchtermanReingold;node.label.show:true;node.label.outline.opacity:15;node.label.shorten:true;edge.label.max-char:20;node.opacity:30;edge.color:blue;' -o out.svg in.gexf
Aug 27, 2018 5:51:14 PM org.gephi.io.processor.plugin.DefaultProcessor process
INFO: # Nodes loaded: 497
Aug 27, 2018 5:51:14 PM org.gephi.io.processor.plugin.DefaultProcessor process
INFO: # Edges loaded: 619
[INFO][GephiCmd]running fruchtermanReingold layout for 10 SECONDS(s)
[INFO][GephiCmd]done
[INFO][GephiCmd]exporting to out.svg

```

## Example of property file.

```
layout.duration=10
layout.time.unit=SECONDS
layout.algorithm=fruchtermanReingold

node.label.show=true
node.label.outline.opacity=1
node.label.outline.size=0.5
node.label.proportinalSize=false
node.label.shorten=true
node.label.color=rgb(255,255,0,100)
node.label.font=Nimbus,PLAIN,12
edge.label.max-char=2
node.label.outline.color=blue
node.label.box=true
node.label.box.color=green
node.label.box.opacity=10

edge.radius=0.1
node.opacity=30
edge.color=lightcoral
edge.thickness=0.5

background-color=gainsboro
directed=true
arrow.size=300
```


## Screenshots

https://twitter.com/yokofakun/status/1034107797439504384

![https://twitter.com/yokofakun/status/1034107797439504384](https://pbs.twimg.com/media/DlnjTj1W4AIBB6B.jpg)


https://twitter.com/yokofakun/status/1034397660189523968

![https://twitter.com/yokofakun/status/1034397660189523968](https://pbs.twimg.com/media/DlrqXqvX4AE6r27.jpg)

https://twitter.com/yokofakun/status/1034440031589552128


![https://twitter.com/yokofakun/status/1034440031589552128](https://pbs.twimg.com/media/DlsRaDCW0AAsotv.jpg)


## Note to self:

part of the code was generated using java2xml && xslt see https://gist.github.com/lindenb/890b9957f768c4bad61b271de7747d07




END_DOC
*/
@Program(name="gephicmd",
	description="Cmd-line oriented for gephi.",
	keywords={"graph","gexf","gephi","visualization"}
	)
public class GephiCmd  extends Launcher {
	private static final Logger LOG = Logger.build(GephiCmd.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-l"},description="list available/default properties and exit with success")
	private boolean list_properties = false;
	@Parameter(names={"-f"},description="zero or more java property files. formatted like a java.util.Properties file. https://docs.oracle.com/cd/E23095_01/Platform.93/ATGProgGuide/html/s0204propertiesfileformat01.html ")
	private List<String> customPropertiesFiles = new ArrayList<>();
	@Parameter(names={"-e"},description="override properties. syntax 'key1:value1;key2:value2;...' ")
	private String customProperties = "";


	private final Map<String,ConfigProperty> config_properties = new  TreeMap<>();
		
	
	private abstract class ConfigProperty
		{
		final String name;
		final String description;
		boolean wasSet=false;
		protected ConfigProperty(final String name,final String description) {
			this.name = name;
			this.description = StringUtil.isBlank(description)?name.replace(".", " "):description;
			if(GephiCmd.this.config_properties.containsKey(name)) {
				throw new IllegalArgumentException("duplicate property "+name);
				}
			GephiCmd.this.config_properties.put(name,this);
			}
		abstract void parse(final String v);
		abstract Object getValue();
		abstract String getDataType();
		@Override
		public String toString() {
			return "\""+name+"\"="+getValue();
			}
		}
	
	private class FloatProperty extends ConfigProperty
		{
		float value;
		protected FloatProperty(final String name,final String description,float def) {
			super(name,description);
			this.value = def;
			}
		@Override String getDataType() { return "float";}
		@Override
		void parse(final String v) {
			try {
				this.value = Float.parseFloat(v);
				super.wasSet=true;
				}
			catch(NumberFormatException err) {
				throw new NumberFormatException("Cannot parse float property "+name+" : "+v);
				}
			}
		@Override
		Object getValue() { return value;}
		}
	
	private class DoubleProperty extends ConfigProperty
		{
		double value;
		protected DoubleProperty(final String name,final String description,double def) {
			super(name,description);
			this.value = def;
			}
		@Override String getDataType() { return "double";}
		@Override
		void parse(final String v) {
			try {
				this.value = Double.parseDouble(v);
				super.wasSet = true;
				}
			catch(NumberFormatException err) {
				throw new NumberFormatException("Cannot parse double property "+name+" : "+v);
				}
			}
		@Override
		Object getValue() { return value;}
		}
	
	private class IntProperty extends ConfigProperty
		{
		int value;
		protected IntProperty(final String name,final String description,int def) {
			super(name,description);
			this.value = def;
			}
		@Override String getDataType() { return "int";}
		@Override
		void parse(final String v) {
			try {
				this.value = Integer.parseInt(v);
				super.wasSet = true;
				}
			catch(final NumberFormatException err) {
				throw new NumberFormatException("Cannot parse ineger property "+name+" : "+v);
				}
			}
		@Override
		Object getValue() { return value;}
		}
	private class BoolProperty extends ConfigProperty
		{
		boolean value;
		protected BoolProperty(final String name,final String description,boolean def) {
			super(name,description);
			this.value = def;
			}
		@Override String getDataType() { return "boolean";}
		@Override
		void parse(final String v) {
			if(v.toLowerCase().equals("true") || v.toLowerCase().equals("t")  || v.toLowerCase().equals("1")) {
				this.value=true;
				super.wasSet = true;
			}
			else if(v.toLowerCase().equals("false") || v.toLowerCase().equals("f") || v.toLowerCase().equals("0")) {
				this.value=false;
				super.wasSet = true;
				}
			else {
			throw new IllegalArgumentException("Cannot parse boolean property "+name+" : "+v);
			}
			}
		@Override
		Object getValue() { return this.value;}
		}
	
	private static final String EDGE_COLOR_INFO = 
			"A edge color is a item from org.gephi.preview.types.EdgeColor.Mode or a color parseable with  "+
			ColorUtils.Converter.OPT_DESC;
			
	
	private class StringProperty extends ConfigProperty
		{
		String value;
		protected StringProperty(final String name,final String description,String def) {
			super(name,description);
			this.value = def;
			}
		@Override String getDataType() { return "string";}
		@Override
		void parse(final String v) {
			super.wasSet = true;
			this.value = v;
			}
		@Override
		Object getValue() { return value;}
		
		Color asAwtColor() {
			final ColorUtils cu = new ColorUtils();
			final Color col = cu.parse(this.value);
			if(col==null) throw new IllegalArgumentException("cannot parse color "+this.value+" for "+this.toString());
			return col;
			}
		
		DependantColor asDependantColor() {
			for(final DependantColor.Mode m : DependantColor.Mode.values()) {
				if(m.name().equalsIgnoreCase(this.value)) return new DependantColor(m);
				}
			return new DependantColor(asAwtColor());
			}
		boolean isEmpty() {
			return StringUtil.isBlank(this.value);
			}
		EdgeColor asEdgeColor() {
			for(final EdgeColor.Mode m : EdgeColor.Mode.values()) {
				if(m.name().equalsIgnoreCase(this.value)) return new EdgeColor(m);
				}
			return new EdgeColor(asAwtColor());
			}
		DependantOriginalColor asDependantOriginalColor()  {
			for(final  DependantOriginalColor.Mode m :  DependantOriginalColor.Mode.values()) {
				if(m.name().equalsIgnoreCase(this.value)) return new DependantOriginalColor(m);
				}
			return new DependantOriginalColor(asAwtColor());
		}
		Font asFont() {
			try {
				final String tokens[]= CharSplitter.COMMA.split(this.value);
				String fontFam = tokens[0];
				int fontFace= Font.PLAIN;
				int fontSize = 12;
				if(tokens.length>1) {
					if(tokens[1].equalsIgnoreCase("plain")) {
						fontFace= Font.PLAIN;
						}
					else if(tokens[1].equalsIgnoreCase("italic")) {
						fontFace= Font.ITALIC;
						}
					else if(tokens[1].equalsIgnoreCase("bold")) {
						fontFace= Font.BOLD;
						}
					if(tokens.length>2) {
						fontSize=Integer.parseInt(tokens[2]);
						}
					}
				return new Font(fontFam,fontFace,fontSize);
				}
			catch(final Exception err) {
				 throw new IllegalArgumentException("cannot parse font "+this.value+" for "+this.toString(),err);
				}
			}
		}
	

	
	private enum Algorithm {forceAtlas,forceAtlas2,fruchtermanReingold,yifanHu,noverlap}
	
	private final IntProperty PROP_LAYOUT_DURATION = new IntProperty("layout.duration", "Layout duration", 1);
	private final StringProperty PROP_LAYOUT_TIME_UNIT = new StringProperty("layout.time.unit", "Layout time unit", TimeUnit.MINUTES.name());
	private final StringProperty PROP_EDGE_DIRECTION_DEFAULT = new StringProperty("edge.direction.default", "Edge direction default ",EdgeDirectionDefault.DIRECTED.name());
	private final BoolProperty PROP_CREATE_MISSING_NODES = new BoolProperty("create.missing.nodes", "",false);
	private final BoolProperty PROP_EXPORT_TRANSPARENT = new BoolProperty("export.transparent.background", "",true);
	private final IntProperty PROP_EXPORT_WIDTH = new IntProperty("export.transparent.width", "",1000);
	private final IntProperty PROP_EXPORT_HEIGHT = new IntProperty("export.transparent.height", "",1000);
	private final IntProperty PROP_EXPORT_MARGIN = new IntProperty("export.transparent.margin", "",1000);
	private final StringProperty PROP_EXPORT_PDF_PAGE_SIZE = new StringProperty("export.pdf.size", "PDF page size. One of http://itextsupport.com/apidocs/iText5/5.5.9/com/itextpdf/text/PageSize.html ","A0");
	private final BoolProperty PROP_EXPORT_PDF_LANDSCAPE = new BoolProperty("export.pdf.landscape", "",false);
	private final BoolProperty PROP_EXPORT_SVG_SCALESTROKE = new BoolProperty("export.svg.scalestrokes", "",false);
	private final StringProperty PROP_LAYOUT_ALGORITHM = new StringProperty("layout.algorithm",
			"Layout algorithm one of "+Arrays.stream(Algorithm.values()).
				map(E->E.name()).collect(Collectors.joining(", "))
			,Algorithm.forceAtlas2.name());

	// from https://raw.githubusercontent.com/gephi/gephi/master/modules/PreviewAPI/src/main/java/org/gephi/preview/api/PreviewProperty.java
	private final BoolProperty PROP_DIRECTED = new BoolProperty(PreviewProperty.DIRECTED,null,true);
	private final StringProperty PROP_BACKGROUND_COLOR = new StringProperty(PreviewProperty.BACKGROUND_COLOR,"General Color property of the background color","");
	private final FloatProperty PROP_VISIBILITY_RATIO = new FloatProperty(PreviewProperty.VISIBILITY_RATIO,"the ratio of the visible graph used in preview. For instance if 0.5 only 50% of nodes items are built",1f);
	private final FloatProperty PROP_MARGIN = new FloatProperty(PreviewProperty.MARGIN,"percentage (0-100) describing the margin size.",5f);
	private final FloatProperty PROP_NODE_BORDER_WIDTH = new FloatProperty(PreviewProperty.NODE_BORDER_WIDTH,"node border size.",0.5f);
	private final StringProperty PROP_NODE_BORDER_COLOR = new StringProperty(PreviewProperty.NODE_BORDER_COLOR,null,"");
	private final FloatProperty PROP_NODE_OPACITY = new FloatProperty(PreviewProperty.NODE_OPACITY,"property between 0-100 which defines the opacity. 100 means opaque",100f);
	private final BoolProperty PROP_NODE_PER_NODE_OPACITY = new BoolProperty(PreviewProperty.NODE_PER_NODE_OPACITY,"Indicating whether or not to use the opacity value defined as part of the Node color",true);
	private final BoolProperty PROP_SHOW_EDGES = new BoolProperty(PreviewProperty.SHOW_EDGES,"Edge Boolean property defining whether to show edges.",true);
	private final FloatProperty PROP_EDGE_THICKNESS = new FloatProperty(PreviewProperty.EDGE_THICKNESS,"Edge Float property for the edge's thickness",1f);
	private final BoolProperty PROP_EDGE_CURVED = new BoolProperty(PreviewProperty.EDGE_CURVED,null,true);
	private final StringProperty PROP_EDGE_COLOR = new StringProperty(PreviewProperty.EDGE_COLOR,"defines the edge color." + EDGE_COLOR_INFO ,"");
	private final FloatProperty PROP_EDGE_OPACITY = new FloatProperty(PreviewProperty.EDGE_OPACITY,"property between 0-100 which defines the opacity.",100f);
	private final BoolProperty PROP_EDGE_RESCALE_WEIGHT = new BoolProperty(PreviewProperty.EDGE_RESCALE_WEIGHT,"defining whether edge's weight should be rescaled between fixed bounds.",false);
	private final FloatProperty PROP_EDGE_RESCALE_WEIGHT_MIN = new FloatProperty(PreviewProperty.EDGE_RESCALE_WEIGHT_MIN,"Edge float property defining the minimum weight when edge weight rescaling is enabled.",0.3f);
	private final FloatProperty PROP_EDGE_RESCALE_WEIGHT_MAX = new FloatProperty(PreviewProperty.EDGE_RESCALE_WEIGHT_MAX,"Edge float property defining the maximum weight when edge weight rescaling is enabled.",100f);
	private final FloatProperty PROP_EDGE_RADIUS = new FloatProperty(PreviewProperty.EDGE_RADIUS,"Edge Float property defining an extra distance between the node and the edge.",0);
	private final FloatProperty PROP_ARROW_SIZE = new FloatProperty(PreviewProperty.ARROW_SIZE,"Arrow Float property defining the arrow size.",3f);
	private final BoolProperty PROP_SHOW_NODE_LABELS = new BoolProperty(PreviewProperty.SHOW_NODE_LABELS,"Node Label Boolean property defining whether to show node labels.",true);
	private final StringProperty PROP_NODE_LABEL_FONT = new StringProperty(PreviewProperty.NODE_LABEL_FONT,"Node Label Font property defining node label's font","");
	private final BoolProperty PROP_NODE_LABEL_PROPORTIONAL_SIZE = new BoolProperty(PreviewProperty.NODE_LABEL_PROPORTIONAL_SIZE,"Node Label Boolean property defining whether to use node's size in label size calculation.",true);
	private final StringProperty PROP_NODE_LABEL_COLOR = new StringProperty(PreviewProperty.NODE_LABEL_COLOR,"color label","");
	private final BoolProperty PROP_NODE_LABEL_SHORTEN = new BoolProperty(PreviewProperty.NODE_LABEL_SHORTEN,"Node Label Boolean property defining whether the label is shortened.",false);
	private final IntProperty PROP_NODE_LABEL_MAX_CHAR = new IntProperty(PreviewProperty.NODE_LABEL_MAX_CHAR,"Edge Label Integer property defining the maximum number of characters.",20);
	private final FloatProperty PROP_NODE_LABEL_OUTLINE_SIZE = new FloatProperty(PreviewProperty.NODE_LABEL_OUTLINE_SIZE,"Edge Label Outline Float property defining the outline size.",0f);
	private final FloatProperty PROP_NODE_LABEL_OUTLINE_OPACITY = new FloatProperty(PreviewProperty.NODE_LABEL_OUTLINE_OPACITY,"Edge Label Outline Float property between 0-100 which defines the opacity. 100 means opaque.",100f);
	private final StringProperty PROP_NODE_LABEL_OUTLINE_COLOR = new StringProperty(PreviewProperty.NODE_LABEL_OUTLINE_COLOR,"label outline color","");
	private final BoolProperty PROP_NODE_LABEL_SHOW_BOX = new BoolProperty(PreviewProperty.NODE_LABEL_SHOW_BOX,null,false);
	private final StringProperty PROP_NODE_LABEL_BOX_COLOR = new StringProperty(PreviewProperty.NODE_LABEL_BOX_COLOR,null,"");
	private final FloatProperty PROP_NODE_LABEL_BOX_OPACITY = new FloatProperty(PreviewProperty.NODE_LABEL_BOX_OPACITY,null,100f);
	private final BoolProperty PROP_SHOW_EDGE_LABELS = new BoolProperty(PreviewProperty.SHOW_EDGE_LABELS,"Edge Label Boolean property defining whether to show edge labels.",false);
	private final StringProperty PROP_EDGE_LABEL_FONT = new StringProperty(PreviewProperty.EDGE_LABEL_FONT,null,"");
	private final StringProperty PROP_EDGE_LABEL_COLOR = new StringProperty(PreviewProperty.EDGE_LABEL_COLOR,null,"");
	private final BoolProperty PROP_EDGE_LABEL_SHORTEN = new BoolProperty(PreviewProperty.EDGE_LABEL_SHORTEN,null,false);
	private final IntProperty PROP_EDGE_LABEL_MAX_CHAR = new IntProperty(PreviewProperty.EDGE_LABEL_MAX_CHAR,"Edge Label Integer property defining the maximum number of characters.",20);
	private final FloatProperty PROP_EDGE_LABEL_OUTLINE_SIZE = new FloatProperty(PreviewProperty.EDGE_LABEL_OUTLINE_SIZE,"Edge Label Outline Float property defining the outline size.",0f);
	private final FloatProperty PROP_EDGE_LABEL_OUTLINE_OPACITY = new FloatProperty(PreviewProperty.EDGE_LABEL_OUTLINE_OPACITY,"Edge Label Outline Float property between 0-100 which defines the opacity. 100 means opaque.",100f);
	private final StringProperty PROP_EDGE_LABEL_OUTLINE_COLOR = new StringProperty(PreviewProperty.EDGE_LABEL_OUTLINE_COLOR,"label outline color","");
	private final StringProperty PROP_EDGE_ARROWS = new StringProperty(PreviewProperty.CATEGORY_EDGE_ARROWS,null,"");
	
	
	/** ForceAtlasLayout **/
	private final DoubleProperty PROP_FORCEATLAST_INERTIA = new DoubleProperty("forceatlas.inertia",null,0.0);
	private final DoubleProperty PROP_FORCEATLAST_REPULSIONSTRENGTH = new DoubleProperty("forceatlas.repulsionstrength",null,0.0);
	private final DoubleProperty PROP_FORCEATLAST_ATTRACTIONSTRENGTH = new DoubleProperty("forceatlas.attractionstrength",null,0.0);
	private final DoubleProperty PROP_FORCEATLAST_MAXDISPLACEMENT = new DoubleProperty("forceatlas.maxdisplacement",null,0.0);
	private final BoolProperty PROP_FORCEATLAST_FREEZEBALANCE = new BoolProperty("forceatlas.freezebalance",null,false);
	private final DoubleProperty PROP_FORCEATLAST_FREEZESTRENGTH = new DoubleProperty("forceatlas.freezestrength",null,0.0);
	private final DoubleProperty PROP_FORCEATLAST_FREEZEINERTIA = new DoubleProperty("forceatlas.freezeinertia",null,0.0);
	private final DoubleProperty PROP_FORCEATLAST_GRAVITY = new DoubleProperty("forceatlas.gravity",null,0.0);
	private final DoubleProperty PROP_FORCEATLAST_COOLING = new DoubleProperty("forceatlas.cooling",null,0.0);
	private final BoolProperty PROP_FORCEATLAST_OUTBOUNDATTRACTIONDISTRIBUTION = new BoolProperty("forceatlas.outboundattractiondistribution",null,false);
	private final BoolProperty PROP_FORCEATLAST_ADJUSTSIZES = new BoolProperty("forceatlas.adjustsizes",null,false);
	private final DoubleProperty PROP_FORCEATLAST_SPEED = new DoubleProperty("forceatlas.speed",null,0.0);
	private final BoolProperty PROP_FORCEATLAST_CONVERGED = new BoolProperty("forceatlas.converged",null,false);
	/** YifanHuLayout **/
	private final IntProperty PROP_YIFANHU_QUADTREEMAXLEVEL = new IntProperty("yifanhu.quadtreemaxlevel",null,0);
	private final FloatProperty PROP_YIFANHU_BARNESHUTTHETA = new FloatProperty("yifanhu.barneshuttheta",null,0f);
	private final FloatProperty PROP_YIFANHU_OPTIMALDISTANCE = new FloatProperty("yifanhu.optimaldistance",null,0f);
	private final FloatProperty PROP_YIFANHU_RELATIVESTRENGTH = new FloatProperty("yifanhu.relativestrength",null,0f);
	private final BoolProperty PROP_YIFANHU_ADAPTIVECOOLING = new BoolProperty("yifanhu.adaptivecooling",null,false);
	private final FloatProperty PROP_YIFANHU_STEPRATIO = new FloatProperty("yifanhu.stepratio",null,0f);
	private final FloatProperty PROP_YIFANHU_CONVERGENCETHRESHOLD = new FloatProperty("yifanhu.convergencethreshold",null,0f);
	//private final FloatProperty PROP_YIFANHU_INITIALSTEP = new FloatProperty("yifanhu.initialstep",null,0f);
	//private final BoolProperty PROP_YIFANHU_CONVERGED = new BoolProperty("yifanhu.converged",null,false);
	/** ForceAtlas2 **/
	private final DoubleProperty PROP_FORCEATLAS2_BARNESHUTTHETA = new DoubleProperty("forceatlas2.barneshuttheta",null,0.0);
	private final DoubleProperty PROP_FORCEATLAS2_EDGEWEIGHTINFLUENCE = new DoubleProperty("forceatlas2.edgeweightinfluence",null,0.0);
	private final DoubleProperty PROP_FORCEATLAS2_JITTERTOLERANCE = new DoubleProperty("forceatlas2.jittertolerance",null,0.0);
	private final BoolProperty PROP_FORCEATLAS2_LINLOGMODE = new BoolProperty("forceatlas2.linlogmode",null,false);
	private final DoubleProperty PROP_FORCEATLAS2_SCALINGRATIO = new DoubleProperty("forceatlas2.scalingratio",null,0.0);
	private final BoolProperty PROP_FORCEATLAS2_STRONGGRAVITYMODE = new BoolProperty("forceatlas2.stronggravitymode",null,false);
	private final DoubleProperty PROP_FORCEATLAS2_GRAVITY = new DoubleProperty("forceatlas2.gravity",null,0.0);
	//private final IntProperty PROP_FORCEATLAS2_THREADSCOUNT = new IntProperty("forceatlas2.threadscount",null,0);
	private final BoolProperty PROP_FORCEATLAS2_OUTBOUNDATTRACTIONDISTRIBUTION = new BoolProperty("forceatlas2.outboundattractiondistribution",null,false);
	private final BoolProperty PROP_FORCEATLAS2_ADJUSTSIZES = new BoolProperty("forceatlas2.adjustsizes",null,false);
	private final BoolProperty PROP_FORCEATLAS2_BARNESHUTOPTIMIZE = new BoolProperty("forceatlas2.barneshutoptimize",null,false);
	/** FruchtermanReingold **/
	private final FloatProperty PROP_FRUCHTERMANREINGOLD_AREA = new FloatProperty("fruchtermanreingold.area",null,0f);
	private final DoubleProperty PROP_FRUCHTERMANREINGOLD_GRAVITY = new DoubleProperty("fruchtermanreingold.gravity",null,0.0);
	private final DoubleProperty PROP_FRUCHTERMANREINGOLD_SPEED = new DoubleProperty("fruchtermanreingold.speed",null,0.0);
	//private final BoolProperty PROP_FRUCHTERMANREINGOLD_CONVERGED = new BoolProperty("fruchtermanreingold.converged",null,false);
	/** NoverlapLayout **/
	private final DoubleProperty PROP_NOVERLAP_SPEED = new DoubleProperty("noverlap.speed",null,0.0);
	private final DoubleProperty PROP_NOVERLAP_RATIO = new DoubleProperty("noverlap.ratio",null,0.0);
	private final DoubleProperty PROP_NOVERLAP_MARGIN = new DoubleProperty("noverlap.margin",null,0.0);
	
	private Stream<Node> streamNodes(final Graph graphPart) {
		return StreamSupport.stream( new IterableAdapter<>(graphPart.getNodes().iterator()).spliterator(),false);
	}
	
@Override
public int doWork(final List<String> args) {
	
	 try
		 {
		 if(this.list_properties) {
				final PrintWriter out = openFileOrStdoutAsPrintWriter(this.outputFile);
				this.config_properties.entrySet().stream().forEach(P->{
					out.println("# "+P.getValue().description+" ("+P.getValue().getDataType()+")");
					out.println(P.getKey()+"="+  P.getValue().getValue());
					});
				out.flush();
				out.close();
				return 0;
				}
		 
		 for(final String filename: this.customPropertiesFiles) {
			final Properties props = new Properties();
			final File file = new File(filename);
			IOUtil.assertFileIsReadable(file);
			final FileReader fr=new FileReader(file);
			props.load(fr);
			fr.close();
			for(final String k:props.stringPropertyNames())
				{
				final ConfigProperty prop = this.config_properties.get(k);
				 if(prop==null) {
					 LOG.error("undefined property "+k+" in "+file);
					 return -1;
				 	}
				 prop.parse(props.getProperty(k));
				}
		 	}
		 
		 
		 for(final String kv: CharSplitter.SEMICOLON.split(this.customProperties)) {
			 if(StringUtil.isBlank(kv)) continue;
			 final String tokens[] = CharSplitter.COLON.split(kv,2);
			 if(tokens.length!=2) {
				 LOG.error("bad property "+kv);
				 return -1;
			 }
			 final ConfigProperty prop = this.config_properties.get(tokens[0]);
			 if(prop==null) {
				 LOG.error("undefined property "+tokens[0]+" in "+this.customProperties);
				 return -1;
			 }
			 prop.parse(tokens[1]);
		 }
		 
		 
		
		 
		 
		 final String gexPath = oneFileOrNull(args);
		 
		 final ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);
		 
		 pc.newProject();
		 final Workspace workspace = pc.getCurrentWorkspace();
		 
	     //Get controllers and models
		 ImportController importController = Lookup.getDefault().lookup(ImportController.class);		 
		 FileImporter fileImporter = importController.getFileImporter(".gexf");
		 
		 if(fileImporter==null) {
	         	LOG.error("FileImporter is null");
	         	return -1;
	         }
	       
        final File file ;
        if(gexPath==null)
        	{
        	LOG.info("reading GEXF from stdin");
        	file = File.createTempFile("tmp.", ".gexf", IOUtils.getDefaultTmpDir());
        	file.deleteOnExit();
        	IOUtils.copyTo(stdin(), file);
        	}
        else
            {
            file = new File(gexPath);
            }
        final Container container = importController.importFile(new java.io.FileReader(file), fileImporter);
        if(container==null || container.getLoader()==null) {
        	LOG.error("importController for "+file+" failed");
        	return -1;
        	}
        container.getLoader().setEdgeDefault(EdgeDirectionDefault.valueOf(PROP_EDGE_DIRECTION_DEFAULT.value));   //Force DIRECTED
        container.getLoader().setAllowAutoNode(PROP_CREATE_MISSING_NODES.value);  //Don't create missing nodes
	    importController.process(container, new DefaultProcessor(), workspace);
	       
	    final GraphModel graphModel = Lookup.getDefault().lookup(GraphController.class).getGraphModel();
	       
       //Layout for 1 minute
        final AutoLayout autoLayout = new AutoLayout(
        		PROP_LAYOUT_DURATION.value,
        		TimeUnit.valueOf(PROP_LAYOUT_TIME_UNIT.value)
        		);
        autoLayout.setGraphModel(graphModel);
	    final Layout layout;
	    switch(Algorithm.valueOf(PROP_LAYOUT_ALGORITHM.value))
	    	{
	    	case yifanHu:
			    {
			    final org.gephi.layout.plugin.force.yifanHu.YifanHuLayout  faLayout = 
			    		new org.gephi.layout.plugin.force.yifanHu.YifanHuProportional().
			    		buildLayout();
			    if(PROP_YIFANHU_QUADTREEMAXLEVEL.wasSet)
					{
					faLayout.setQuadTreeMaxLevel(PROP_YIFANHU_QUADTREEMAXLEVEL.value);
					}
				if(PROP_YIFANHU_BARNESHUTTHETA.wasSet)
					{
					faLayout.setBarnesHutTheta(PROP_YIFANHU_BARNESHUTTHETA.value);
					}
				if(PROP_YIFANHU_OPTIMALDISTANCE.wasSet)
					{
					faLayout.setOptimalDistance(PROP_YIFANHU_OPTIMALDISTANCE.value);
					}
				if(PROP_YIFANHU_RELATIVESTRENGTH.wasSet)
					{
					faLayout.setRelativeStrength(PROP_YIFANHU_RELATIVESTRENGTH.value);
					}
				if(PROP_YIFANHU_ADAPTIVECOOLING.wasSet)
					{
					faLayout.setAdaptiveCooling(PROP_YIFANHU_ADAPTIVECOOLING.value);
					}
				if(PROP_YIFANHU_STEPRATIO.wasSet)
					{
					faLayout.setStepRatio(PROP_YIFANHU_STEPRATIO.value);
					}
				if(PROP_YIFANHU_CONVERGENCETHRESHOLD.wasSet)
					{
					faLayout.setConvergenceThreshold(PROP_YIFANHU_CONVERGENCETHRESHOLD.value);
					}
				/*
				if(PROP_YIFANHU_INITIALSTEP.wasSet)
					{
					faLayout.setInitialStep(PROP_YIFANHU_INITIALSTEP.value);
					}
				if(PROP_YIFANHU_CONVERGED.wasSet)
					{
					faLayout.setConverged(PROP_YIFANHU_CONVERGED.value);
					}*/
				layout =    faLayout;
				break;
			    }
	    	case forceAtlas:
	    		{
	    		final org.gephi.layout.plugin.forceAtlas.ForceAtlasLayout  faLayout = 
	    				new org.gephi.layout.plugin.forceAtlas.ForceAtlas().
	    					buildLayout();

	    		if(PROP_FORCEATLAST_INERTIA.wasSet)
		    		{
		    		faLayout.setInertia(PROP_FORCEATLAST_INERTIA.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_REPULSIONSTRENGTH.wasSet)
		    		{
		    		faLayout.setRepulsionStrength(PROP_FORCEATLAST_REPULSIONSTRENGTH.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_ATTRACTIONSTRENGTH.wasSet)
		    		{
		    		faLayout.setAttractionStrength(PROP_FORCEATLAST_ATTRACTIONSTRENGTH.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_MAXDISPLACEMENT.wasSet)
		    		{
		    		faLayout.setMaxDisplacement(PROP_FORCEATLAST_MAXDISPLACEMENT.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_FREEZEBALANCE.wasSet)
		    		{
		    		faLayout.setFreezeBalance(PROP_FORCEATLAST_FREEZEBALANCE.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_FREEZESTRENGTH.wasSet)
		    		{
		    		faLayout.setFreezeStrength(PROP_FORCEATLAST_FREEZESTRENGTH.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_FREEZEINERTIA.wasSet)
		    		{
		    		faLayout.setFreezeInertia(PROP_FORCEATLAST_FREEZEINERTIA.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_GRAVITY.wasSet)
		    		{
		    		faLayout.setGravity(PROP_FORCEATLAST_GRAVITY.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_COOLING.wasSet)
		    		{
		    		faLayout.setCooling(PROP_FORCEATLAST_COOLING.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_OUTBOUNDATTRACTIONDISTRIBUTION.wasSet)
		    		{
		    		faLayout.setOutboundAttractionDistribution(PROP_FORCEATLAST_OUTBOUNDATTRACTIONDISTRIBUTION.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_ADJUSTSIZES.wasSet)
		    		{
		    		faLayout.setAdjustSizes(PROP_FORCEATLAST_ADJUSTSIZES.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_SPEED.wasSet)
		    		{
		    		faLayout.setSpeed(PROP_FORCEATLAST_SPEED.value);
		    		}
	
	
		    	if(PROP_FORCEATLAST_CONVERGED.wasSet)
		    		{
		    		faLayout.setConverged(PROP_FORCEATLAST_CONVERGED.value);
		    		}
	

	    		layout =    faLayout;
	    		break;
	    		}
	    	case forceAtlas2:
	    		{
	    		final ForceAtlas2 faLayout = new ForceAtlas2Builder().buildLayout();
	    		if(PROP_FORCEATLAS2_BARNESHUTTHETA.wasSet)
		    		{
		    		faLayout.setBarnesHutTheta(PROP_FORCEATLAS2_BARNESHUTTHETA.value);
		    		}
		    	if(PROP_FORCEATLAS2_EDGEWEIGHTINFLUENCE.wasSet)
		    		{
		    		faLayout.setEdgeWeightInfluence(PROP_FORCEATLAS2_EDGEWEIGHTINFLUENCE.value);
		    		}
		    	if(PROP_FORCEATLAS2_JITTERTOLERANCE.wasSet)
		    		{
		    		faLayout.setJitterTolerance(PROP_FORCEATLAS2_JITTERTOLERANCE.value);
		    		}
		    	if(PROP_FORCEATLAS2_LINLOGMODE.wasSet)
		    		{
		    		faLayout.setLinLogMode(PROP_FORCEATLAS2_LINLOGMODE.value);
		    		}
		    	if(PROP_FORCEATLAS2_SCALINGRATIO.wasSet)
		    		{
		    		faLayout.setScalingRatio(PROP_FORCEATLAS2_SCALINGRATIO.value);
		    		}
		    	if(PROP_FORCEATLAS2_STRONGGRAVITYMODE.wasSet)
		    		{
		    		faLayout.setStrongGravityMode(PROP_FORCEATLAS2_STRONGGRAVITYMODE.value);
		    		}
		    	if(PROP_FORCEATLAS2_GRAVITY.wasSet)
		    		{
		    		faLayout.setGravity(PROP_FORCEATLAS2_GRAVITY.value);
		    		}
		    	/*if(PROP_FORCEATLAS2_THREADSCOUNT.wasSet)
		    		{
		    		faLayout.setThreadsCount(PROP_FORCEATLAS2_THREADSCOUNT.value);
		    		}*/
		    	if(PROP_FORCEATLAS2_OUTBOUNDATTRACTIONDISTRIBUTION.wasSet)
		    		{
		    		faLayout.setOutboundAttractionDistribution(PROP_FORCEATLAS2_OUTBOUNDATTRACTIONDISTRIBUTION.value);
		    		}
		    	if(PROP_FORCEATLAS2_ADJUSTSIZES.wasSet)
		    		{
		    		faLayout.setAdjustSizes(PROP_FORCEATLAS2_ADJUSTSIZES.value);
		    		}
		    	if(PROP_FORCEATLAS2_BARNESHUTOPTIMIZE.wasSet)
		    		{
		    		faLayout.setBarnesHutOptimize(PROP_FORCEATLAS2_BARNESHUTOPTIMIZE.value);
		    		}

		        layout =    faLayout;
		        break;
	    		}
	    	case fruchtermanReingold:
	    		{
	    		final FruchtermanReingold  faLayout = new FruchtermanReingoldBuilder().buildLayout();
	    		if(PROP_FRUCHTERMANREINGOLD_AREA.wasSet)
		    		{
		    		faLayout.setArea(PROP_FRUCHTERMANREINGOLD_AREA.value);
		    		}
		    	if(PROP_FRUCHTERMANREINGOLD_GRAVITY.wasSet)
		    		{
		    		faLayout.setGravity(PROP_FRUCHTERMANREINGOLD_GRAVITY.value);
		    		}
		    	if(PROP_FRUCHTERMANREINGOLD_SPEED.wasSet)
		    		{
		    		faLayout.setSpeed(PROP_FRUCHTERMANREINGOLD_SPEED.value);
		    		}

	    		layout =    faLayout;
	    		break;
	    		}
	    	case noverlap:
	    		{
	    		final org.gephi.layout.plugin.noverlap.NoverlapLayout  faLayout = 
	    				(org.gephi.layout.plugin.noverlap.NoverlapLayout)
	    				new org.gephi.layout.plugin.noverlap.NoverlapLayoutBuilder().
	    				buildLayout();

	    		if(PROP_NOVERLAP_SPEED.wasSet)
	    			{
	    			faLayout.setSpeed(PROP_NOVERLAP_SPEED.value);
	    			}
	    		if(PROP_NOVERLAP_RATIO.wasSet)
	    			{
	    			faLayout.setRatio(PROP_NOVERLAP_RATIO.value);
	    			}
	    		if(PROP_NOVERLAP_MARGIN.wasSet)
	    			{
	    			faLayout.setMargin(PROP_NOVERLAP_MARGIN.value);
	    			}

    			layout =    faLayout;
	    		break;
	    		}
	    	default:
	    		{
	    		LOG.info("undefined algoritm "+PROP_LAYOUT_ALGORITHM+". Available are: " +
	    				Arrays.stream(Algorithm.values()).map(E->E.name()).collect(Collectors.joining(", "))
	    				);
	    		return -1;
	    		}
	    	}
        
	   
	    
        
	    	
	       autoLayout.addLayout(layout, 1f);
	       LOG.info("running "+PROP_LAYOUT_ALGORITHM.value+" layout for "+ PROP_LAYOUT_DURATION.value+" "+PROP_LAYOUT_TIME_UNIT.value+"(s)");
	       autoLayout.execute();
	       LOG.info("done");
	        
            final Graph graphPart = graphModel.getGraphVisible();
            final double min_node_x = streamNodes(graphPart).mapToDouble(N->N.x()-N.size()).min().orElse(0);
            final double max_node_x = streamNodes(graphPart).mapToDouble(N->N.x()+N.size()).max().orElse(0);
            final double min_node_y = streamNodes(graphPart).mapToDouble(N->N.y()-N.size()).min().orElse(0);
            final double max_node_y = streamNodes(graphPart).mapToDouble(N->N.y()+N.size()).max().orElse(0);
            LOG.warning("width "+min_node_x+":"+max_node_x+" "+(max_node_x-min_node_x));
            LOG.warning("height "+min_node_y+":"+max_node_y+" "+(max_node_y-min_node_y));
            
            
            PreviewController previewController = Lookup.getDefault().lookup(PreviewController.class);
            final PreviewModel previewModel = previewController.getModel();
            
            
            if(PROP_EDGE_COLOR.wasSet && !PROP_EDGE_COLOR.isEmpty())
            	{
            	 previewModel.getProperties().putValue(PROP_EDGE_COLOR.name,
            			 PROP_EDGE_COLOR.asEdgeColor()
            			 );	
            	}
            
            
         // handle properties for colors
            for(final StringProperty colorProp : new StringProperty[] {
            		PROP_NODE_LABEL_OUTLINE_COLOR,
            		PROP_NODE_LABEL_BOX_COLOR,
            	})
	            {
            	if(!colorProp.wasSet) continue;
	            if(colorProp.isEmpty()) continue;
	            LOG.debug("setting "+colorProp+" "+colorProp.asDependantColor());
	            previewModel.getProperties().putValue(colorProp.name,colorProp.asDependantColor());
	            }
            
            
            // handle properties for colors
            for(final StringProperty colorProp : new StringProperty[] {
            		PROP_BACKGROUND_COLOR,
            		PROP_NODE_BORDER_COLOR,            		
            		PROP_EDGE_LABEL_COLOR,
            		PROP_EDGE_LABEL_OUTLINE_COLOR
            	})
	            {
            	if(!colorProp.wasSet) continue;
	            if(colorProp.isEmpty()) continue;
	            LOG.debug("setting "+colorProp);
	            previewModel.getProperties().putValue(colorProp.name,colorProp.asAwtColor());	
	            }
            
            for(final StringProperty colorProp : new StringProperty[] {
            		PROP_NODE_LABEL_COLOR,
            	})
	            {
            	if(!colorProp.wasSet) continue;
	            if(colorProp.isEmpty()) continue;
	            LOG.debug("setting "+colorProp);
	            previewModel.getProperties().putValue(colorProp.name,colorProp.asDependantOriginalColor());	
	            }
            
         // handle properties for font
            for(final StringProperty fontProp : new StringProperty[] {
            		 PROP_NODE_LABEL_FONT,
            		 PROP_EDGE_LABEL_FONT
            	})
	            {
            	if(!fontProp.wasSet) continue;
	            if(fontProp.isEmpty()) continue;
	            LOG.debug("setting "+fontProp);
	            previewModel.getProperties().putValue(fontProp.name,fontProp.asFont());	
	            }
           
            
            //handle other previewModel properties
            for(final ConfigProperty other : new ConfigProperty[] {
            		PROP_EDGE_DIRECTION_DEFAULT,
            		PROP_SHOW_EDGES,
            		PROP_EDGE_THICKNESS,
            		PROP_EDGE_CURVED,
            		PROP_EDGE_OPACITY,
            		PROP_EDGE_RESCALE_WEIGHT,
            		PROP_EDGE_RESCALE_WEIGHT_MIN,
            		PROP_EDGE_RESCALE_WEIGHT_MAX,
            		PROP_EDGE_RADIUS,
            		PROP_SHOW_EDGE_LABELS,
            		PROP_EDGE_LABEL_SHORTEN,
            		PROP_EDGE_LABEL_MAX_CHAR,
            		PROP_EDGE_LABEL_OUTLINE_SIZE,
            		PROP_EDGE_LABEL_OUTLINE_OPACITY,
            		//
            		PROP_NODE_BORDER_WIDTH,
            		PROP_NODE_OPACITY,
            		PROP_NODE_PER_NODE_OPACITY,
            		PROP_SHOW_NODE_LABELS,
            		PROP_NODE_LABEL_PROPORTIONAL_SIZE,
            		PROP_NODE_LABEL_SHORTEN,
            		PROP_NODE_LABEL_MAX_CHAR,
            		PROP_NODE_LABEL_OUTLINE_SIZE,
            		PROP_NODE_LABEL_OUTLINE_OPACITY,
            		PROP_NODE_LABEL_SHOW_BOX,
            		PROP_NODE_LABEL_BOX_OPACITY,
            		//
            		PROP_ARROW_SIZE
            	})
	            {
            	if(!other.wasSet) continue;
            	LOG.debug("setting "+other);
	            previewModel.getProperties().putValue(other.name,other.getValue());	
	            }
            
            previewController.refreshPreview();
            
            final ExportController ec = Lookup.getDefault().lookup(ExportController.class);
            
            if(this.outputFile==null) 
            	{
            	LOG.info("exporting to stdout/PNG");
            	}
            else
            	{
            	LOG.info("exporting to "+this.outputFile);
            	}
            
            //Export full graph
            if(this.outputFile==null || this.outputFile.getName().toLowerCase().endsWith(".png"))
            	{
            	  final PNGExporter pngExporter = (PNGExporter) ec.getExporter("png");
                  
            	  int png_width = PROP_EXPORT_WIDTH.value; 
            	  if(png_width > (int)(max_node_x-min_node_x))
            	  	{
            		  png_width = (int)(max_node_x-min_node_x);
            		 LOG.warn("adjusting png width to "+png_width);
            	  	}
            	  
            	  int png_height = PROP_EXPORT_HEIGHT.value;
            	  if(png_height > (int)(max_node_y-min_node_y))
	          	  	{
            		png_height = (int)(max_node_y-min_node_y);
	          		LOG.warn("adjusting png height to "+png_height);
	          	  	}
            	  
                  pngExporter.setHeight(png_height);
                  pngExporter.setWidth(png_width);
                  pngExporter.setMargin(PROP_EXPORT_MARGIN.value);
                  pngExporter.setTransparentBackground(PROP_EXPORT_TRANSPARENT.value);
                  pngExporter.setWorkspace(workspace);
                  if(this.outputFile==null) {
                	 final File tmpFile = File.createTempFile("gephi.", ".png",IOUtils.getDefaultTmpDir());
                	 tmpFile.deleteOnExit();
                	 ec.exportFile(tmpFile, pngExporter);
                	 IOUtils.copyTo(tmpFile, stdout());
                	 tmpFile.delete();
                  	}
                  else
                  	{
                	ec.exportFile(this.outputFile, pngExporter);
                  	}
                  
            	}
            else if( this.outputFile.getName().toLowerCase().endsWith(".pdf"))
	            {
	            final PDFExporter pdfExporter = (PDFExporter) ec.getExporter("pdf");
	            final com.itextpdf.text.Rectangle pageSize =  PageSize.getRectangle(PROP_EXPORT_PDF_PAGE_SIZE.value);
	            if(pageSize==null ) {
	            	LOG.error("bad pdf page size "+PROP_EXPORT_PDF_PAGE_SIZE);
	            	return -1;
	            }
	            pdfExporter.setPageSize(pageSize);
	            pdfExporter.setLandscape(PROP_EXPORT_PDF_LANDSCAPE.value);
	            pdfExporter.setWorkspace(workspace);
	            
	            ec.exportFile(this.outputFile, pdfExporter);
	            /*
	            ByteArrayOutputStream baos = new ByteArrayOutputStream();
	            ec.exportStream(baos, pdfExporter);
	            byte[] pdf = baos.toByteArray();
	            
	            final FileOutputStream fos = new FileOutputStream(this.outputFile);
	            fos.write(pdf);
	            fos.flush();
	            fos.close();
	            */
	            }
            else  if( this.outputFile.getName().toLowerCase().endsWith(".svg"))
	        	{
	        	final SVGExporter exporter = (SVGExporter) ec.getExporter("svg");
	        	exporter.setWorkspace(workspace);
	        	exporter.setScaleStrokes(PROP_EXPORT_SVG_SCALESTROKE.value);
	        	ec.exportFile(this.outputFile, exporter);
	        	}
            else  if( this.outputFile.getName().toLowerCase().endsWith(".gexf"))
            	{
            	final ExporterGEXF exporter = (ExporterGEXF) ec.getExporter("gexf");
            	exporter.setExportPosition(true); 
            	exporter.setWorkspace(workspace);
            	ec.exportFile(this.outputFile, exporter);
            	}
            else
            	{
            	LOG.error("bad output format for "+this.outputFile+". suffix must be .svg, .png, .pdf or .gexf .");
            	return -1;
            	}
            LOG.info("Completed.");
            return 0;
		 }
	 catch(final Throwable err)
	 	{
		 LOG.error(err);
		return -1;
	 	}
	}
public static void main(String[] args) {
	new GephiCmd().instanceMainWithExit(args);
	}
}
