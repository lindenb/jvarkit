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
package com.github.lindenb.jvarkit.tools.jfx.gatkjfx;

import javafx.scene.Scene;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.TextField;
import javafx.scene.layout.GridPane;
import javafx.stage.Stage;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;


public class VariantFiltrationJfx extends AbstractGatkJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	
	@FXML private FileChooserPane mask;
	@FXML private TextField maskName;
 	@FXML private CheckBox filterNotInMask;
 	@FXML private TextField maskExtension;
 	
 	@FXML private CheckBox invalidatePreviousFilters;
 	@FXML private CheckBox invertFilterExpression;
 	@FXML private CheckBox invertGenotypeFilterExpression;
 	@FXML private CheckBox missingValuesInExpressionsShouldEvaluateAsFailing;
 	@FXML private CheckBox setFilteredGtToNocall;
 	@FXML private GridPane expressionGridPane;
 	@FXML private GridPane gtExpressionGridPane;
	private List<Expression> expression=new ArrayList<>();
	private List<Expression> gtexpressions=new ArrayList<>();
 	
 	private static class Expression
 		{
 		final TextField expression;
 		final TextField name;
 		final boolean gtExpr;
 		Expression(boolean gtExpr) {
 			this.gtExpr=gtExpr;
 			this.expression = new  TextField();
 			this.expression.setPrefColumnCount(40);
 			this.expression.setPromptText("JEXL expression");
 			this.name = new  TextField();
 			this.name.setPrefColumnCount(10);
 			this.name.setPromptText("Filter Name");
 			}
 		void fill(final List<String> args) throws JFXException
 			{
 			final String ex=this.expression.getText().trim();
 			final String n=this.name.getText().trim();
 			if(ex.isEmpty() && n.isEmpty()) return;
 			if(ex.isEmpty() && !n.isEmpty())throw new JFXException("No expression for filter ("+n+")");
 			if(!ex.isEmpty() && n.isEmpty())throw new JFXException("No name for expression ("+ex+")");
 			args.add(gtExpr?"--genotypeFilterName":"--filterName");
 			args.add(n);
 			args.add(gtExpr?"--genotypeFilterExpression":"--filterExpression");
 			args.add(ex);
 			}
 		}
 	
	public VariantFiltrationJfx() {
	}
	
	@Override
	protected String getAnalysisType() {
		return "VariantFiltration";
		}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Scene scene = new Scene(fxmlLoad("VariantFiltrationJfx.fxml"));
        
		stage.setScene(scene);
		
		
		 this.expressionGridPane.add(new Label("Filter Expression"),0,0);
		 this.expressionGridPane.add(new Label("Filter Name"),1,0);
		for(int y=1;y<=10;++y){
		final Expression ex=new Expression(false);
		this.expression.add(ex);
		
		this.expressionGridPane.add(ex.expression,0,y);
		this.expressionGridPane.add(ex.name,1,y);
	    
		}
		
		 this.gtExpressionGridPane.add(new Label("GTFilter Expression"),0,0);
		 this.gtExpressionGridPane.add(new Label("GTFilter Name"),1,0);
		for(int y=1;y<=10;++y){
		final Expression ex=new Expression(true);
		this.gtexpressions.add(ex);
		
		this.gtExpressionGridPane.add(ex.expression,0,y);
		this.gtExpressionGridPane.add(ex.name,1,y);
	    
		}

		
		
        super.start(stage);
    	}
	
	
	public static void main(String[] args)
		{
		launch(args);
		}
	
	@Override
	protected  List<String> buildArgs() throws JFXException {
		final List<String> args= super.buildArgs();
		new OptionBuilder(inputvcf,"--variant").fill(args);
		new OptionBuilder(outputvcf,"-o").fill(args);
		
		for(final Expression ex:this.expression)
			{
			ex.fill(args);
			}
		for(final Expression ex:this.gtexpressions)
			{
			ex.fill(args);
			}

		if( mask.getSelectedFile()!=null){
			final String makName= maskName.getText().trim();
			if(makName.isEmpty()) throw new JFXException("Mask name is empty");
			new OptionBuilder(mask,"--mask").fill(args);
			new OptionBuilder(maskName,"--maskName").fill(args);
			new OptionBuilder(filterNotInMask,"--filterNotInMask").fill(args);
			new OptionBuilder(maskExtension,"--maskExtension").itemClass(Integer.class).fill(args);
		}
		
		new OptionBuilder(this.invalidatePreviousFilters,"--invalidatePreviousFilters").fill(args);		
		new OptionBuilder(this.invertFilterExpression,"--invertFilterExpression").fill(args);		
		new OptionBuilder(this.invertGenotypeFilterExpression,"--invertGenotypeFilterExpression").fill(args);		

		new OptionBuilder(this.missingValuesInExpressionsShouldEvaluateAsFailing,"--missingValuesInExpressionsShouldEvaluateAsFailing").fill(args);		
		new OptionBuilder(this.setFilteredGtToNocall,"--setFilteredGtToNocall").fill(args);		

		
		return args;
	}
	
	
	

	
	}
