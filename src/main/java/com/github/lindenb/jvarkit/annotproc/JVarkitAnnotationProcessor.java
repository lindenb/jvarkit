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
package com.github.lindenb.jvarkit.annotproc;

import java.io.File;
import java.io.FileReader;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import javax.annotation.processing.AbstractProcessor;
import javax.annotation.processing.Filer;
import javax.annotation.processing.RoundEnvironment;
import javax.lang.model.SourceVersion;
import javax.lang.model.element.Element;
import javax.lang.model.element.ElementKind;
import javax.lang.model.element.TypeElement;
import javax.lang.model.type.DeclaredType;
import javax.lang.model.type.TypeKind;
import javax.lang.model.type.TypeMirror;
import javax.tools.FileObject;
import javax.tools.StandardLocation;



import com.github.lindenb.jvarkit.util.jcommander.Program;

public class JVarkitAnnotationProcessor extends AbstractProcessor{
	private static final Logger LOG = Logger.getLogger(JVarkitAnnotationProcessor.class.getSimpleName());
	
	
	@Override
	public SourceVersion getSupportedSourceVersion() {
		return SourceVersion.latest();
		}
	
	@Override
	public Set<String> getSupportedAnnotationTypes() {
		return Arrays.asList(
				IncludeSourceInJar.class,
				Program.class
				).stream().
			map(C->C.getName()).
			collect(Collectors.toSet())
			;
		}
	
	private void copySource(final javax.lang.model.element.Element element)
			{
			final boolean include_super_class = false;
			if(element==null || element.getKind()!=ElementKind.CLASS) return;
			final String thisDir = System.getProperty("jvarkit.this.dir","src/main/java");	
			if(thisDir==null || thisDir.isEmpty()) {
				LOG.warning("[PROC] jvarkit.basedir is not defined");
				return ;
				}
			try 
				{
				do
					{
					String className= element.toString();
					if(className==null ||  className.isEmpty() || className.equals("java.lang.Object") ) return;
		
					final int dollar  = className.indexOf('$');
					if(dollar!=-1) className = className.substring(0,dollar);
					final File javaFile = new File(thisDir+File.separator+className.replace('.',File.separatorChar) +".java");
		
					if(!javaFile.exists()) {
						LOG.warning("[PROC] File not found: "+javaFile);
						break;
						}
					final Filer filer = super.processingEnv.getFiler();					
					final String packageName;
					final String fileName;
					int dot= className.lastIndexOf('.');
					if(dot==-1)
						{
						packageName = "";
						fileName = className;
						}
					else
						{
						packageName = className.substring(0,dot);
						fileName = className.substring(dot+1);
						}
					final FileObject fo=filer.createResource(StandardLocation.CLASS_OUTPUT,
							packageName,
							fileName+".java");
					
					if( new File(fo.getName()).exists()) {
						LOG.info("[PROC] ## skip "+ javaFile+ " because it exists");
						break;
						}
					
					LOG.info("[PROC] ## Copying "+ javaFile+ " -> "+ fo.getName());
					try(FileReader reader = new FileReader(javaFile)) {
						try(Writer writer=fo.openWriter()) {
							final char  buffer[]=new char[1024];
							int nRead;
							while((nRead=reader.read(buffer))!=-1)
								{
								writer.write(buffer,0,nRead);
								}
							writer.flush();
							}
						}
					} while(false);
				}
			catch(final Exception err)
				{
				LOG.warning(err.getMessage());
				}
			/* include super class */
			if(include_super_class) {
				final TypeMirror superclass = ((TypeElement) element).getSuperclass();
				if(superclass==null || superclass.getKind() == TypeKind.NONE) {
					return ;
					}
				final DeclaredType kind = (DeclaredType) superclass;
				copySource(kind.asElement());
				}
			}
	
		
		@Override
		public boolean process(final Set<? extends TypeElement> annotations,
			final RoundEnvironment roundEnv
			) {
			LOG.info("[PROC] process");
			final Set<Element> set = new HashSet<>();
			set.addAll(roundEnv.getElementsAnnotatedWith(IncludeSourceInJar.class));
			set.addAll(roundEnv.getElementsAnnotatedWith(Program.class));
			
			set.stream().
					filter(E->E.getKind()==ElementKind.CLASS).
					filter(E-> E.getAnnotation(Program.class)!=null || E.getAnnotation(IncludeSourceInJar.class)!=null).
					forEach(E->{					
						copySource(E);
						});
			return true;
			}
		}
