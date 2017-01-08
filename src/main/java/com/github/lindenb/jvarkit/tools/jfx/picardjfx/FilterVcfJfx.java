package com.github.lindenb.jvarkit.tools.jfx.picardjfx;

import javafx.scene.Scene;
import javafx.scene.control.Spinner;
import javafx.stage.Stage;
import picard.vcf.filter.FilterVcf;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.jfx.components.FileChooserPane;

import javafx.fxml.*;
import javafx.scene.*;


public class FilterVcfJfx extends AbstractPicardJfxApplication
	{
	@FXML
	private FileChooserPane inputvcf;
	@FXML
	private FileChooserPane outputvcf;
	@FXML
	private Spinner<Double> min_ab;
	@FXML
	private Spinner<Integer> min_dp;
	@FXML
	private Spinner<Integer> min_gq;
	@FXML
	private Spinner<Double> max_fs;
	@FXML
	private Spinner<Double> min_qd;
	@FXML
	private FileChooserPane javascript;
	
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
    		
			if( inputvcf==null) throw new java.io.IOException("inputvcf is null");
			
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
	protected  List<String> buildArgs() throws JFXException {
		final List<String> args= new ArrayList<>();
		new OptionBuilder(inputvcf,"I=").fill(args);
		new OptionBuilder(outputvcf,"O=").fill(args);
		new OptionBuilder(min_ab,"MIN_AB=").fill(args);
		new OptionBuilder(min_dp,"MIN_DP=").fill(args);
		new OptionBuilder(min_gq,"MIN_GQ=").fill(args);
		new OptionBuilder(max_fs,"MAX_FS=").fill(args);
		new OptionBuilder(min_qd,"MIN_QD=").fill(args);
		new OptionBuilder(javascript,"JS=").fill(args);
		return args;
	}
	
	
	

	
	}
