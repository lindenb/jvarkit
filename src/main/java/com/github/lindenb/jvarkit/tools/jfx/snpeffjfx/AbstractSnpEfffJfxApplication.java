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
package com.github.lindenb.jvarkit.tools.jfx.snpeffjfx;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.List;

import com.github.lindenb.jvarkit.tools.jfx.AbstractJfxApplication;

import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.stage.Stage;

public abstract class AbstractSnpEfffJfxApplication extends AbstractJfxApplication {
	private AbstractSnpEffRunner runningThread=null; 

	protected AbstractSnpEfffJfxApplication(){
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		stage.setTitle(getClass().getSimpleName());
		super.start(stage);
		}


	protected abstract class AbstractSnpEffRunner
		implements Runnable
		{
		protected final Class<?> snpEffClass;
		private final File fileout;
		private PrintStream stdout=null;
		protected final String args[];
		private final AbstractSnpEfffJfxApplication owner;
		AbstractSnpEffRunner(final String className,final List<String> args,final File fileout) throws JFXException
			{
			this.owner = AbstractSnpEfffJfxApplication.this;
			this.args = args.toArray(new String[args.size()]);
			this.fileout = fileout;
			try {
				this.snpEffClass = Class.forName(className);
			} catch (ClassNotFoundException err) {
				throw new JFXException(err);
			}			
			
			}
		
		protected abstract Object createInstance() throws JFXException;
		protected abstract boolean isExitSuccess(Object o) throws JFXException;
		
		@Override
		public void run() {
			try {
				final Object snpEffInstance = createInstance();
				
				final Method run = snpEffClass.getMethod("run");
				this.stdout = new PrintStream(this.fileout);
				System.setOut(this.stdout);
				System.setErr(AbstractSnpEfffJfxApplication.this.printToConsole);
				final Object is_ok =run.invoke(snpEffInstance);
				AbstractSnpEfffJfxApplication.this.printToConsole.println("Done ("+is_ok+")");
				this.stdout.flush();
				this.stdout.close();
				this.stdout=null;
				if(isExitSuccess(is_ok) )
					{
					try {
						Platform.runLater(new Runnable() {
							@Override
							public void run() {
								//if(AbstractSnpEffRunner.this != owner.runningThread) return;
								final Alert alert = new Alert(AlertType.CONFIRMATION);
								alert.setContentText("Completed:"+Arrays.toString(AbstractSnpEffRunner.this.args));
								
								alert.showAndWait();
								owner.runningThread=null;
							}
						});
					} catch (Exception e) {
						e.printStackTrace(AbstractJfxApplication.realStderr);
						}
					}
				else
					{
					try {
						Platform.runLater(new Runnable() {
							@Override
							public void run() {
								//if(AbstractSnpEffRunner.this != owner.runningThread) return;
								final Alert alert = new Alert(AlertType.ERROR);
								alert.setContentText("Failure:"+Arrays.toString(AbstractSnpEffRunner.this.args));
	
								owner.runningThread=null;
							}
						});
					} catch (Exception e) {
						e.printStackTrace(AbstractJfxApplication.realStderr);
						}
					}
			} catch (Exception e) {
				e.printStackTrace(AbstractSnpEfffJfxApplication.realStderr);
				}
			finally
				{
				if(this.stdout!=null) 
					{
					this.stdout.flush();
					this.stdout.close();
					}
				}
			}
		}
}
