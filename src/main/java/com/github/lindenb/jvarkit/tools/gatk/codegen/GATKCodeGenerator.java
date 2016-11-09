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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

public class GATKCodeGenerator extends AbstractGATKCodeGenerator {
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(GATKCodeGenerator.class);
	public static class Utils
	{
		
	}
public GATKCodeGenerator() {
	
}
@Override
public Collection<Throwable> call() throws Exception {
	final List<String> args = super.getInputFiles();
	if(args.size()!=1) {
		return wrapException("expected one and only one argument.");
	}
	if(super.velocityTemplate==null) {
		return wrapException("undefined option -" + OPTION_VELOCITYTEMPLATE);
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
		
		IOUtil.assertFileIsReadable(super.velocityTemplate);
		pw = openFileOrStdoutAsPrintWriter();
		VelocityContext context = new VelocityContext();
		
		context.put("json", jsonDoc);
		context.put("utils", new Utils());
		final VelocityEngine ve = new VelocityEngine();
        ve.setProperty(RuntimeConstants.RESOURCE_LOADER, "file");
        ve.setProperty("file.resource.loader.class",FileResourceLoader.class.getName());
        ve.setProperty("file.resource.loader.path",super.velocityTemplate.getParent());
        ve.init();
        final Template template = ve.getTemplate(super.velocityTemplate.getName());
        template.merge( context, pw);
        pw.flush();
		
		return RETURN_OK;
		}
	catch(Exception err)
		{
		return wrapException(err);
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
