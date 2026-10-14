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
package com.github.lindenb.jvarkit.tools.genome2svg;

import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.springframework.context.support.FileSystemXmlApplicationContext;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.BamCoverageTrack;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.BasesTrack;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.GffTrack;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.KnownGeneTrack;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.Track;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.VcfTrack;
import com.github.lindenb.jvarkit.tools.genome2svg.beans.WiggleTrack;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

@Program(name="genome2svg",
description="static genome browser as SVG",
keywords={"vcf","svg","xml","visualization"},
modificationDate="20240120",
creationDate="20261013",
generate_doc = false
)
public class GenomeToSvg extends Launcher {
	private static final Logger LOG=Logger.of(GenomeToSvg.class);
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputPath=null;
	@Parameter(names={"-r","--interval","--region"},description=IntervalParser.OPT_DESC,required = true)
	private String intervalStr = null;
	@Parameter(names={"-w","--width"},description="Image width")
	private double image_width = 1000;
	
	private BasesTrack _force_compile1 =null;
	private KnownGeneTrack _force_compile2 =null;
	private BamCoverageTrack _force_compile3 =null;
	private GffTrack _force_compile4 =null;
	private VcfTrack _force_compile5 =null;
	private WiggleTrack _force_compile6=null;

	abstract class ConfigNode {
		ConfigNode parent;
		ConfigNode(ConfigNode parent) {
			this.parent = parent;
			}
		ConfigObjectNode asObject() { throw new IllegalArgumentException();}
		ConfigArrayNode asArray() { throw new IllegalArgumentException();}
		}
	
	class ConfigObjectNode extends ConfigNode  implements Iterable<Map.Entry<String, ConfigNode>>  {
		private final Map<String,ConfigNode> hash=new LinkedHashMap<>();
		ConfigObjectNode(ConfigNode parent) {
			super(parent);
			}
		public Set<String> keySet() { return hash.keySet();}
		public boolean has(String key) { return hash.containsKey(key);}
		public ConfigNode get(String key) { return hash.get(key);}
		@Override
		public Iterator<Entry<String, ConfigNode>> iterator() {
			return hash.entrySet().iterator();
			}
		@Override ConfigObjectNode asObject() { return this;}
		}
	
	class ConfigArrayNode extends ConfigNode implements Iterable<ConfigNode> {
		private final List<ConfigNode> array = new ArrayList<>();
		ConfigArrayNode(ConfigNode parent) {
			super(parent);
			}
		@Override
		public Iterator<ConfigNode> iterator() {
			return array.iterator();
			}
		ConfigArrayNode asArray() { return this;}
		}
	
	class StringNode extends ConfigNode {
		private final String value;
		StringNode(ConfigNode parent,String value) {
			super(null);
			this.value=value;
			}
		}
	
	private ConfigNode fromJSON(ConfigNode parent,JsonElement root) {
		if(root.isJsonNull()) return null;
		if(root.isJsonObject()) {
			final ConfigObjectNode o = new ConfigObjectNode(parent);
			root.getAsJsonObject().entrySet().stream().forEach(KV->{
				if(!KV.getValue().isJsonNull()) {
					o.hash.put(KV.getKey(), fromJSON(o,KV.getValue()));
					}
				});
			return o;
			}
		else if(root.isJsonArray()) {
			final ConfigArrayNode o = new ConfigArrayNode(parent);
			root.getAsJsonArray().forEach(E->{
				if(!E.isJsonNull()) {
					o.array.add(fromJSON(o, E));
					}
				});
			return o;
			}
		else
			{
			return new StringNode(parent, root.getAsJsonPrimitive().getAsString());
			}
		}
	
	
	
	
	private static JsonObject mix(final JsonObject root, final JsonObject o2) {
		for(Map.Entry<String,JsonElement> kv2 : o2.entrySet()) {
			if(root.has(kv2.getKey())) {
				final JsonElement x1 = root.get(kv2.getKey());
				final JsonElement x2 = kv2.getValue();
				if((x1 instanceof JsonObject) && (x2 instanceof JsonObject)) {
					root.add(kv2.getKey(), mix(x1.getAsJsonObject(),x2.getAsJsonObject()));
					}
				else
					{
					root.add(kv2.getKey(), kv2.getValue());
					}
				}
			else
				{
				root.add(kv2.getKey(), kv2.getValue());
				}
			}
		return root;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final JsonObject root = new JsonObject();
			final JsonParser jsonParser=new JsonParser();
			for(Path jsonPath: IOUtils.unrollPaths(args)) {
				try (Reader r= Files.newBufferedReader(jsonPath) {
					root = mix(root,jsonParser.parse(r).getAsJsonObject());
					}
				}
			
			final SVGContext svgContext = new SVGContext();
			svgContext.image_width = this.image_width;
			
			configuration =  new FileSystemXmlApplicationContext(args.toArray(new String[args.size()]));
			
			/*
		    DefaultConversionService service = new DefaultConversionService();
		    service.addConverter(Converter);
			ConfigurableEnvironment environment = new StandardEnvironment();
			environment.setConversionService(service);
			configuration.setEnvironment(environment);
			*/
			svgContext.loc = new IntervalParser().
					apply(this.intervalStr).
					get();

		
			
			final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			final DocumentBuilder db = dbf.newDocumentBuilder();
			final Document svgDom = db.newDocument();
			svgContext.svgDom = svgDom;
			final Element svgRoot = svgContext.element("svg");
			svgRoot.setAttribute("width", String.valueOf(1+this.image_width + svgContext.left_margin));
			svgDom.appendChild(svgRoot);
			
			svgContext.styleNode = svgContext.element("style");
			svgRoot.appendChild(svgContext.styleNode);

			
			svgContext.defsNode = svgContext.element("defs");
			svgRoot.appendChild(svgContext.defsNode);
			svgContext.tracksNode = svgContext.element("g");
			svgContext.tracksNode.setAttribute("transform", "translate("+svgContext.left_margin+",0)");
			svgRoot.appendChild(svgContext.tracksNode);
			final Element frame =  svgContext.element("rect");
			svgContext.tracksNode.appendChild(frame);
			
			frame.setAttribute("style", "fill:white;stroke:black;");
			frame.setAttribute("x", "0");
			frame.setAttribute("y", "0");
			frame.setAttribute("width", String.valueOf(svgContext.image_width+ svgContext.left_margin-1));
			

			for(Object ot: configuration.getBean("tracks",List.class)) {
				final Track track = Track.class.cast(ot);
				if(!track.isVisible()) continue;
				double prev_y = svgContext.y;
				track.paint(svgContext);
				if(svgContext.y>prev_y) {
					final Element frame2 =  svgContext.element("rect");
					frame2.setAttribute("style","stroke:darkgray;fill:none;");
					frame2.setAttribute("x", String.valueOf(-svgContext.left_margin));
					frame2.setAttribute("y",String.valueOf(prev_y));
					frame2.setAttribute("width", String.valueOf(svgContext.image_width+ svgContext.left_margin-1));
					frame2.setAttribute("height", String.valueOf(svgContext.y-prev_y));
					svgContext.tracksNode.appendChild(frame2);
					}
				}
			
			
			svgRoot.setAttribute("height", String.valueOf(svgContext.y));
			frame.setAttribute("height", String.valueOf(svgContext.y));
			
			TransformerFactory.newInstance().newTransformer().
				transform(new DOMSource(svgDom),
				this.outputPath==null?
				new StreamResult(stdout()):
				new StreamResult(outputPath.toFile())
				);
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			configuration = null;
			}
		}
	
	public static void main(String[] args) {
		new GenomeToSvg().instanceMain(args);
	}
}
