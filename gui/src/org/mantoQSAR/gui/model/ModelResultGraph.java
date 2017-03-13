/* This file is part of mantoQSAR.

mantoQSAR - Quantitative structure-activity relationship descriptor 
			calculation and modeling for biomolecules.
			
Copyright (C) 2016  Jörg Kittelmann


mantoQSAR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, 
or any later version.

mantoQSAR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mantoQSAR. If not, see <http://www.gnu.org/licenses/>.
*/

package org.mantoQSAR.gui.model;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JPanel;
import javax.swing.UIManager;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.TickUnitSource;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.labels.CustomXYToolTipGenerator;
import org.jfree.chart.labels.StandardXYItemLabelGenerator;
import org.jfree.chart.labels.XYItemLabelGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.mantoQSAR.core.ObservationSet;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.model.ModelGroup;
import org.mantoQSAR.core.model.ModelList;
import org.mantoQSAR.core.model.ModelStatic;
import org.mantoQSAR.core.util.ColorStatic;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;



public class ModelResultGraph extends JPanel implements PropertyChangeListener {

    private final Logger logger;

    private ChartPanel chartPanel;
     ModelController modelControl;
     ModelGroup modelGroup; 
     ModelList modelList; 
 
    
    double[] scale = new double[2];
    double scaleMax;

    private boolean logScale = false;
    private boolean modelDetail = false;

    public ModelResultGraph() {

        modelControl = ModelController.getInstance(); 
        modelGroup = modelControl.getModelGroup(); 
        modelList = modelGroup.getModelList();
        
        
        logger = LoggerFactory.getLogger(ModelResultGraph.class);

        scale[0] = 0.0;
        scale[1] = 50.0;
        scaleMax = 100.0;
 
        this.setBackground(Color.WHITE);

        initChart();
        initListeners(); 

    }

   
    
   private void initChart() {


        JFreeChart chart = ChartFactory.createScatterPlot(
                null, // chart title
                "x", // domain (x) axis label
                "y", // range (y) axis label
                null, // initial series
                PlotOrientation.VERTICAL, // orientation
                true, // include legend
                true, // tooltips?
                false // URLs?
        );
        
        chart.setBackgroundPaint(Color.white);
       
        chartPanel = new ChartPanel(chart, false) {

            @Override
            public final Dimension getPreferredSize() {
                Dimension d = super.getPreferredSize();
                Dimension prefSize = null;
                Component c = getParent();
                if (c == null) {
                    prefSize = new Dimension(
                            (int) d.getWidth(), (int) d.getHeight());
                } else if (c.getWidth() > d.getWidth()
                        && c.getHeight() > d.getHeight()) {
                    prefSize = c.getSize();
                } else {
                    prefSize = d;
                }
                int w = (int) prefSize.getWidth() - 100;
                int h = (int) prefSize.getHeight() - 20;

                int s = (w > h ? h : w);

                return new Dimension(s, s);
            }
        };

        chartPanel.setBackground(Color.WHITE);
        chartPanel.setPreferredSize(new Dimension(640, 480));

        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL;

        this.add(chartPanel, c);
    }
    

    private void createChart(){

   
        System.out.println(ColorStatic.PURPLE + "call to draw chart" + ColorStatic.RESET);
        JFreeChart chart = ChartFactory.createScatterPlot(
                "main model", 
                "experimental value ", 
                "predicted value ", 
                null, // initial series
                PlotOrientation.VERTICAL, // orientation
                true, // include legend
                true, // tooltips?
                false // URLs?
        );

        // set chart background and legend border
        chart.setBackgroundPaint(Color.white);
        chart.getLegend().setFrame(BlockBorder.NONE);
        
        // set fonts to theme font
        Font themeFont = UIManager.getDefaults().getFont("TextField.font");

        TextTitle title = chart.getTitle();
        title.setFont(themeFont);

        // set a few custom plot features
        XYPlot plot = (XYPlot) chart.getPlot();
        // plot.setBackgroundPaint(new Color(0xffffe0));
        plot.setBackgroundAlpha(new Float(0.0));

        plot.setDomainGridlinesVisible(true);
        plot.setDomainGridlinePaint(Color.lightGray);
        plot.setRangeGridlinePaint(Color.lightGray);

        // render shapes and lines
        List<Color> colorList = new ArrayList<>();
        colorList.add(Color.GRAY);
        plot.setRenderer(newDataRenderer(colorList));

        
        if (logScale == true) {

            LogarithmicAxis logAxisRange = new LogarithmicAxis("predicted value ");
            // logAxisRange.setBase(10);
            logAxisRange.setLog10TickLabelsFlag(false);
            logAxisRange.setMinorTickMarksVisible(true);
            logAxisRange.setAllowNegativesFlag(true);
            logAxisRange.setRange(this.scale[0], this.scale[1]);

            plot.setRangeAxis(logAxisRange);

            LogarithmicAxis logAxisDomain = new LogarithmicAxis("experimental value ");
            logAxisDomain.setLog10TickLabelsFlag(false);
            logAxisDomain.setMinorTickMarksVisible(true);
            logAxisDomain.setAllowNegativesFlag(true);
            logAxisDomain.setRange(this.scale[0], this.scale[1]);

            plot.setDomainAxis(logAxisDomain);

        } else {

            // set the plot's axes to display integers
            TickUnitSource ticks = NumberAxis.createIntegerTickUnits();
            NumberAxis domain = (NumberAxis) plot.getDomainAxis();
            domain.setLabelFont(themeFont);
            domain.setStandardTickUnits(ticks);

            NumberAxis range = (NumberAxis) plot.getRangeAxis();
            range.setLabelFont(themeFont);
            range.setStandardTickUnits(ticks);
        }
        
        this.chartPanel.setChart(chart);
        
        this.drawCenterLine();
        
            if (this.modelDetail == true) {
                this.drawModelDetail();
            }

            this.drawModelMain();
            this.drawPredictionMain();
            this.revalidate();
    }

    private XYLineAndShapeRenderer newDataRenderer(List<Color> colorList) {

        XYLineAndShapeRenderer renderer
                = new XYLineAndShapeRenderer(false, true);
        renderer.setBaseShapesVisible(true);
        renderer.setBaseShapesFilled(false);

        Ellipse2D circle = new Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0);

        for (int i = 0; i < colorList.size(); i++) {
            renderer.setSeriesPaint(i, colorList.get(i));
            renderer.setSeriesShape(i, circle);
        }

        Stroke stroke = new BasicStroke(
                1.5f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL);
        renderer.setBaseOutlineStroke(stroke);

        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(2);
        XYItemLabelGenerator generator
                = new StandardXYItemLabelGenerator(
                        StandardXYItemLabelGenerator.DEFAULT_ITEM_LABEL_FORMAT,
                        format, format);
        renderer.setBaseItemLabelGenerator(generator);
        renderer.setBaseItemLabelsVisible(false);

        return renderer;

    }
    private XYLineAndShapeRenderer newCenterLineRenderer() {

        XYLineAndShapeRenderer renderer
                = new XYLineAndShapeRenderer(true, false);
        renderer.setBaseShapesVisible(false);
        renderer.setBaseShapesFilled(false);

        renderer.setSeriesPaint(0, Color.LIGHT_GRAY);

        Stroke stroke = new BasicStroke(
                3f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL);
        renderer.setBaseOutlineStroke(stroke);

        return renderer;

    }


    public void drawCenterLine() {

        XYPlot plot = (XYPlot) this.chartPanel.getChart().getPlot();
  
        XYSeries dataAdded = new XYSeries("center line");
        dataAdded.add(this.scale[0], this.scale[0]);
        dataAdded.add(this.scale[1], this.scale[1]);

        XYSeriesCollection dataCollection = new XYSeriesCollection();
        dataCollection.addSeries(dataAdded);

        plot.setDataset(0, (XYDataset) dataCollection);
        plot.setRenderer(0, this.newCenterLineRenderer());
    }

    
    public void drawLabel(List<String> label, String type){
        
        XYPlot plot = (XYPlot) this.chartPanel.getChart().getPlot();
        
        
        int pos = 1; 
        
        if ("modelmain".equalsIgnoreCase(type)) {pos = 2;}
        if ("modeldetail".equalsIgnoreCase(type)) {pos = 3;}
        if ("predictiondetail".equalsIgnoreCase(type)) {pos = 4;}
        if ("predictionmain".equalsIgnoreCase(type)) {pos = 5;}

        CustomXYToolTipGenerator generator = new CustomXYToolTipGenerator();
        generator.addToolTipSeries(label);
        plot.getRenderer().setSeriesToolTipGenerator(pos, generator);
    }
    
    
    public void drawSeries(List<List<Double>> X, List<List<Double>> Y, String type) {

        XYPlot plot = (XYPlot) this.chartPanel.getChart().getPlot();

        if (X.size() != Y.size()) {
            System.out.println(ColorStatic.RED + "Draw vectors X/Y are not of same length. \n"
                    + " Xvector is " + X.size() + " and Y is" + Y.size() + ColorStatic.RESET);
            return;
        }

        List<Double> x;
        List<Double> y;

        XYSeriesCollection dataCollection = new XYSeriesCollection();

        String title = "unknown";
        Color color = Color.GRAY;
        int pos = 1; 

        if ("modelmain".equalsIgnoreCase(type)) {
            title = "model data";
            color = Color.DARK_GRAY;
            pos = 2;
        }
        if ("modeldetail".equalsIgnoreCase(type)) {
            title = "single model data";
            color = Color.LIGHT_GRAY;
            pos = 3;
        }

        if ("predictiondetail".equalsIgnoreCase(type)) {
            title = "prediction data";
            color = Color.YELLOW;
            pos = 4;
        }

        if ("predictionmain".equalsIgnoreCase(type)) {
            title = "single prediction data";
            color = Color.RED;
            pos = 5;
        }

        XYSeries dataAdded = new XYSeries(title);

        int n = 0;
        for (int i = 0; i < X.size(); i++) {
            x = X.get(i);
            y = Y.get(i);

            System.out.println(ColorStatic.BLUE + "Drawing x " + x.toString() + ColorStatic.RESET);
            System.out.println(ColorStatic.BLUE + "Drawing y " + y.toString() + ColorStatic.RESET);

            if (x.isEmpty() || y.isEmpty()) {
                continue;
            }

            if (x.size() != y.size()) {
                logger.error("Draw vector x/y are not of same length.");
                continue;
            }

            for (int j = 0; j < x.size(); j++) {

                Double xs = x.get(j);
                Double ys = y.get(j);

                if (ys != null) {
                    dataAdded.add(xs, ys);

                    if (xs > scale[1]) {
                        if (xs < scaleMax) {
                            scale[1] = xs;
                        } else {
                            scale[1] = scaleMax;
                        }
                    }

                    if (ys > scale[1]) {
                        if (ys < scaleMax) {
                            scale[1] = ys;
                        } else {
                            scale[1] = scaleMax;
                        }
                    }
                }
            }
        }

        List<Color> colorList = new ArrayList<>();

        dataCollection.addSeries(dataAdded);
        colorList.add(color);

        plot.setDataset(pos, (XYDataset) dataCollection);
        plot.setRenderer(pos, this.newDataRenderer(colorList));

        plot.getRangeAxis().setRange(this.scale[0], this.scale[1]);
        plot.getDomainAxis().setRange(this.scale[0], this.scale[1]);
        
        this.repaint();
    }
        
    /*
     draw ModelGroup content including single model property/response and 
     prediction property/responses pairs
     */
    public void drawModelDetail() {

        ModelGroup modelGroup = this.modelControl.getModelGroup();

        List<List<Double>> modelPropertyDetail = modelGroup.getModelPropertyDetail(modelGroup.lowerModelRSquare, modelGroup.lowerPredictRSquare);
        List<List<Double>> modelResponseDetail = modelGroup.getModelResponseDetail(modelGroup.lowerModelRSquare, modelGroup.lowerPredictRSquare);

        if (modelPropertyDetail != null
                && modelResponseDetail != null) {

            drawSeries(modelPropertyDetail, modelResponseDetail, "modelDetail");
        } else {
            System.out.println(ColorStatic.PURPLE + "No models to draw.");
        }
    }

    public void drawModelMain() {
        
        ModelGroup modelGroup = this.modelControl.getModelGroup();

        List<Double> modelProperty = modelGroup.getModelProperty();
        List<Double> modelResponse = modelGroup.getModelResponse();

        if (modelProperty != null && modelResponse != null) {
            
            List<List<Double>> Xmodel = new ArrayList<>();
            Xmodel.add(modelProperty);
            List<List<Double>> Ymodel = new ArrayList<>();
            Ymodel.add(modelResponse);

            List<String> label = new ArrayList<>();
            List<ObservationSet> osl = Screen.getInstance().getModelObservationSetList(); 
            
            for (ObservationSet os : osl) {
                label.add(os.getName() + " " + os.getCondition().getpH());
            }

            this.drawSeries(Xmodel, Ymodel, "modelMain");
            this.drawLabel(label, "modelMain");
        }
        
    }   
        
    public void drawPredictionMain() {

        ModelGroup modelGroup = this.modelControl.getModelGroup();
        List<Double> predictProperty = modelGroup.getPredictProperty();
        List<Double> predictResponse = modelGroup.getPredictResponse();

        if (predictProperty != null
                && predictResponse != null) {

            List<List<Double>> Xpred = new ArrayList<>();
            Xpred.add(predictProperty);
            List<List<Double>> Ypred = new ArrayList<>();
            Ypred.add(predictResponse);

            List<String> label = new ArrayList<>();
            List<ObservationSet> osl = Screen.getInstance().getPredictObservationSetList(); 
            
            for (ObservationSet os : osl) {
                label.add(os.getName() + " " + os.getCondition().getpH());
            }
            
            this.drawSeries(Xpred, Ypred, "predictionMain");
            this.drawLabel(label, "predictionMain");
            
            System.out.println("Prediction printed to graph.");
        } else {
            System.out.println("Prediction not printed to graph.");
        }
    }

    public boolean isLogScale() {
        return logScale;
    }

    public void setLogScale(boolean logS) {
        if (this.logScale != logS) {
            this.logScale = logS;
            this.createChart();
        }
    }

    public boolean isModelDetail() {
        return modelDetail;
    }

    public void setModelDetail(boolean modelD) {
        if (this.modelDetail != modelD) {
            this.modelDetail = modelD;
            createChart();
        }
    }

     protected final void initListeners(){
        modelControl.addPropertyChangeListener(new PropertyChangeListener() {

            @Override
            public void propertyChange(PropertyChangeEvent evt) {
            

        if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CALCULATED)) {
            createChart();
        }
        if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CHANGED)) {
            createChart();
        }
            }
        });
    }
    
    
    @Override
    public void propertyChange(PropertyChangeEvent evt) {


        if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CALCULATED)) {
            this.createChart();
        }
        if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CHANGED)) {
            this.createChart();

        }
    }
}
