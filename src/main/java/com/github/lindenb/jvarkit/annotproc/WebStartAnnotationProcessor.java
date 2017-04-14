package com.github.lindenb.jvarkit.annotproc;

import java.util.Collections;
import java.util.Set;

import javax.annotation.processing.AbstractProcessor;
import javax.annotation.processing.Filer;
import javax.annotation.processing.ProcessingEnvironment;
import javax.annotation.processing.RoundEnvironment;
import javax.annotation.processing.SupportedSourceVersion;
import javax.lang.model.SourceVersion;
import javax.lang.model.element.TypeElement;
import javax.tools.FileObject;
import javax.tools.StandardLocation;

import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@SupportedSourceVersion(SourceVersion.RELEASE_8)
public class WebStartAnnotationProcessor extends AbstractProcessor{
	private static final Logger LOG = Logger.build(WebStartAnnotationProcessor.class).make();
	
	private void show() {
	}
	@Override
	public synchronized void init(ProcessingEnvironment processingEnv) {
	
	   super.init(processingEnv);
	}
	@Override
	public Set<String> getSupportedAnnotationTypes() {
		return Collections.singleton(Program.class.getName());
		}
	
	@Override
	public boolean process(final Set<? extends TypeElement> annotations,
		final RoundEnvironment roundEnv
		) {
		//LOG.info("################### PROCESSING ON THE WAY");
		show();
		if(roundEnv.processingOver() && "A".equals("B"))
			{
			final Filer filer = super.processingEnv.getFiler();
			LOG.info(System.getProperty("jvarkit.libs.jars"));
			try 
				{
				FileObject fo=filer.createResource(StandardLocation.CLASS_OUTPUT, "", "META-INF/TEST" );
				fo.openWriter().append("Hello word").close();
				LOG.info(fo.getName());
				}
			catch(Exception err)
				{
				LOG.warn(err);
				}
			}
		return false;
	}
	
	

}
