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
package com.github.lindenb.jvarkit.jfx.components;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import htsjdk.samtools.util.StringUtil;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.scene.control.Alert;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TextArea;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.stage.Window;
import javafx.stage.FileChooser.ExtensionFilter;

public class JfxUtils {
	
	public static class DialogBuilder
		{
		private Object err;
		private String title=null;
		private String headerText=null;
		private AlertType alertType = AlertType.WARNING;
		
		public DialogBuilder cause(final Object o) { this.err=o; return this;}
		public DialogBuilder title(final String s) { this.title=s; return this;}
		public DialogBuilder header(final String s) { this.headerText=s; return this;}
		public DialogBuilder type(final AlertType alertType) { this.alertType=alertType; return this;}
		
		private String getHeaderText() {
			return StringUtil.isBlank(this.headerText)?
				   alertType.name():
				   this.headerText;
			}
		private String getTitle() {
			return StringUtil.isBlank(this.title)?
				   getHeaderText():
				   this.title;
			}

		
		public void show(final Window owner) {
			final Alert alert = new Alert(this.alertType);
			alert.setTitle(this.getTitle());
		
			if(this.err!=null && this.err instanceof Throwable)
				{
				alert.setHeaderText(getHeaderText());
				final Throwable error=(Throwable) err;
				alert.setContentText(String.valueOf(error.getMessage()));
				error.printStackTrace();
		
				final StringWriter sw = new StringWriter();
				final PrintWriter pw = new PrintWriter(sw);
				error.printStackTrace(pw);
				final String exceptionText = sw.toString();
		
				final Label label = new Label("The exception stacktrace was:");
				final TextArea textArea = new TextArea(exceptionText);
				textArea.setEditable(false);
				textArea.setWrapText(false);
		
		
				textArea.setMaxWidth(Double.MAX_VALUE);
				textArea.setMaxHeight(Double.MAX_VALUE);
				GridPane.setVgrow(textArea, Priority.ALWAYS);
				GridPane.setHgrow(textArea, Priority.ALWAYS);
		
				final BorderPane expContent = new BorderPane(new ScrollPane(textArea));
				expContent.setTop(label);
				expContent.setMinHeight(500);
				expContent.setMinWidth(800);
				expContent.setPrefWidth(1000);
				expContent.setPrefHeight(900);
		
				// Set expandable Exception into the dialog pane.
				alert.getDialogPane().setExpandableContent(expContent);
				}
			else
				{
				final String msg=String.valueOf(err);
				alert.setHeaderText(msg);
				alert.setContentText(msg);
				alert.getDialogPane().setExpandableContent(new Label(msg));
				}
			alert.showAndWait();
			}
		}
	
	
	public static DialogBuilder dialog() {
		return new DialogBuilder();
		}	
	
	public static List<ExtensionFilter> getIndexedVcfExtensionFilters() {
	    return  Arrays.asList(
	    		new ExtensionFilter("Indexed VCF File", "*.vcf", "*.vcf.gz"),
	    		new ExtensionFilter("Tabix-Indexed VCF Files", "*.vcf.gz"),
	    		new ExtensionFilter("Tribble-Indexed VCF File", "*.vcf")
	    		);
		}
	
	public static <T,R> TableColumn<T,R> makeTableColumn(final String label,final Function<T,R> supplier)
		{
	    final TableColumn<T,R>  col = new TableColumn<>(label);
	    col.setCellValueFactory(param-> new ReadOnlyObjectWrapper<R>(supplier.apply(param.getValue())));
	    return col;
		}

}
