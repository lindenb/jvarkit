package com.github.lindenb.jvarkit.tools.jfx.picardjfx;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextArea;
import javafx.scene.layout.*;
import javafx.stage.Stage;
import picard.vcf.filter.FilterVcf;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.ResourceBundle;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;
import javafx.scene.*;
import javafx.geometry.Rectangle2D;
import javafx.stage.Screen;


public class FilterVcfJfx extends AbstractPicardJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	@FXML
	private Spinner<Double> min_ab;
	
	
	public FilterVcfJfx() {
		super(FilterVcf.class);
	}
	
	@Override
	public void start(Stage stage) throws Exception {
		final Parent root;
		
		try
			{
			java.net.URL in= getClass().getResource("FilterVcfJfx.fxml");
			if(in==null) throw new java.io.IOException("cannot get resource");
			FXMLLoader loader = new FXMLLoader(in);
			loader.setController(this);
			root = loader.load();
    		
			}
    	catch(Exception err)
    		{
    		err.printStackTrace();
    		throw err;
    		}
    		
		
       
        Scene scene = new Scene(root);
        stage.setScene(scene);
       
        
        super.start(stage);
    	}
	
	
	public static void main(String[] args)
		{
		launch(args);
		}
	
	@Override
	protected  List<String> buildArgs() {
		final List<String> args= new ArrayList<>();
		new OptionBuilder(inputvcf,"I=").fill(args);
		new OptionBuilder(outputvcf,"O=").fill(args);
		new OptionBuilder(min_ab,"MIN_AB=").fill(args);
		
		return args;
	}
	
	
	

	
	}
