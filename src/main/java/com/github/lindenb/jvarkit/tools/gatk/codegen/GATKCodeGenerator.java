package com.github.lindenb.jvarkit.tools.gatk.codegen;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Collection;
import java.util.List;

import org.apache.velocity.Template;
import org.apache.velocity.VelocityContext;
import org.apache.velocity.app.VelocityEngine;
import org.apache.velocity.runtime.RuntimeConstants;
import org.apache.velocity.runtime.resource.loader.FileResourceLoader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="",description="Gatk code generator")
public class GATKCodeGenerator extends Launcher {
private static final Logger LOG = Logger.build(GATKCodeGenerator.class).make();
@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private File outputFile = null;

@Parameter(names="-T",description="Velocity template",required=true)
private File velocityTemplate = null;

public static class Utils
	{
		
	}
public GATKCodeGenerator() {
	
}
@Override
	public int doWork(List<String> args) {
	if(args.size()!=1) {
		LOG.error("expected one and only one argument.");
		return -1;
	}
	if(this.velocityTemplate==null) {
		LOG.error("undefined option -T");
		return -1;
	}
	String command = args.get(0);
	Reader jsonReader = null;
	PrintWriter pw = null;
	JsonElement jsonDoc=null;
	try
		{
		final JsonParser parser=new JsonParser();
		
		try {
			File jsonF = new File(command);
			if(!jsonF.exists() && jsonF.isFile()) {
				jsonReader = new FileReader(jsonF);
			}
		} catch (IOException e) {
			jsonReader=null;
		}
		
		if(jsonReader==null) {
			if(command.endsWith(".php")) {
				command+=".json";
			}
			
			if(command.startsWith("org_")) {
				command = "https://software.broadinstitute.org/gatk/gatkdocs/"+command;
			}

			jsonReader = new InputStreamReader(IOUtils.openURIForReading(command));
		}
		
		jsonDoc = parser.parse(jsonReader);
		jsonReader.close();
		
		LOG.info("URL is :" + command);
		
		IOUtil.assertFileIsReadable(this.velocityTemplate);
		pw = openFileOrStdoutAsPrintWriter(this.outputFile);
		VelocityContext context = new VelocityContext();
		
		context.put("json", jsonDoc);
		context.put("utils", new Utils());
		final VelocityEngine ve = new VelocityEngine();
        ve.setProperty(RuntimeConstants.RESOURCE_LOADER, "file");
        ve.setProperty("file.resource.loader.class",FileResourceLoader.class.getName());
        ve.setProperty("file.resource.loader.path",this.velocityTemplate.getParent());
        ve.init();
        final Template template = ve.getTemplate(this.velocityTemplate.getName());
        template.merge( context, pw);
        pw.flush();
		
		return RETURN_OK;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(pw);
		CloserUtil.close(jsonReader);
		}
	}
public static void main(String[] args) {
	new GATKCodeGenerator().instanceMainWithExit(args);
}
}
