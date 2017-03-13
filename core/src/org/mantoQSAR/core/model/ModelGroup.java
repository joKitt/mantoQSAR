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
package org.mantoQSAR.core.model;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.mantoQSAR.core.math.MathUtil;
import org.mantoQSAR.core.math.Matrix;
import org.mantoQSAR.core.util.ColorStatic;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ModelGroup {

    PropertyChangeSupport changes;

    ModelList modelList;
    Logger logger;

    Matrix modelMatrix = null;
    Matrix predictMatrix = null;

    List<Double> modelProperty = null;
    List<Double> predictProperty;

    List<List<Double>> predictResponse;
    List<List<Double>> modelResponse;

    // set values for development purpose
    public Double randSelect = 1.0;
    public Double randPredict = 0.15;

    public double lowerModelRSquare = 0.5;
    public double lowerPredictRSquare = 0.2;

    Boolean normalize = true;
    Boolean logScale = true;
    Double[] normalizationFactor;

    List<Double> coefVal;
    List<Integer> coefCount;

    public ModelGroup(int numM) {

        this.logger = LoggerFactory.getLogger(ModelGroup.class);
        this.modelList = new ModelList();
        this.predictProperty = new ArrayList<>();
        this.modelProperty = new ArrayList<>();
        this.predictResponse = new ArrayList<>();
        this.modelResponse = new ArrayList<>();

        for (int i = 0; i < numM; i++) {
            ReplacementModel repMod = new ReplacementModel();
            this.modelList.addModel(repMod);
        }

        changes = new PropertyChangeSupport(this);
    }

    private void setCoefSetting() {

        int n = this.getModelDescriptor().getColumnCount();
        coefVal = new ArrayList<>(n);
        coefCount = new ArrayList<>(n);

        for (int i = 0; i < n; i++) {
            coefVal.add(0.0);
            coefCount.add(0);
        }
    }

    /*
     returns  list of integers representing observation distribution; 
     null and and 0: samples are excluded according to fraction by randSelect;  
     1 - used for model generation;
     2 - used for internal model test;
     */
    private List<Integer> randGroupObservation(int length) {

        Random rand = new Random();
        List<Integer> listSelect = new ArrayList<>();
        for (int i = 0; i < length; i++) {
            if (rand.nextDouble() <= this.randSelect) {

                if (rand.nextDouble() <= this.randPredict) {
                    listSelect.add(2);
                } else {
                    listSelect.add(1);
                }
            } else {
                listSelect.add(-1);
            }
        }

        System.out.println("randGroupObservation " + listSelect.toString());
        return listSelect;
    }

    public List<Double> sampleList(List<Double> data, List<Integer> selector) {

        List<Double> newData = new ArrayList<>();

        for (int i = 0; i < data.size(); i++) {
            if (selector.get(i) != null) {
                newData.add(data.get(i));
            }
        }
        return newData;
    }

    public List<Double> completeList(List<Double> sublist, List<Integer> selector) {
        List<Double> newData = new ArrayList<>();

        int mi = 0;
        for (Integer selector1 : selector) {
            if (selector1 != null) {
                newData.add(sublist.get(mi));
                mi++;
            } else {
                newData.add(null);
            }
        }
        return newData;
    }

    public boolean processInputsSingle()
            throws InterruptedException, ExecutionException {

        for (int i = 0; i < modelList.getModelList().size(); i++) {

            ReplacementModel mc = (ReplacementModel) modelList.getModel(i);

            mc.setDescriptor(getModelDescriptor());
            mc.setProperty(setScale(modelProperty));
            mc.setSelector(randGroupObservation(modelMatrix.getRowCount()));
            mc.calcModel();
            modelList.setModel(i, mc);

            /*  System.out.println(ColorStatic.ANSI_BLUE + " model " + i + " response:" + ColorStatic.RESET);
             System.out.println(ColorStatic.ANSI_BLUE + " model property: " + mc.getModelProperty() + ColorStatic.RESET);
             System.out.println(ColorStatic.ANSI_BLUE + " model response: " + mc.getModelResponse() + ColorStatic.RESET);
             System.out.println(ColorStatic.ANSI_BLUE + " model prediction: " + mc.getModelResponse(predictMatrix) + ColorStatic.RESET);
             */
            changes.firePropertyChange(ModelStatic.EVENT_MODEL_CALCULATED, null, null);
        }
        return true;
    }

    public List<ModelClass> processInputsSingle(List<ModelClass> inputs)
            throws InterruptedException, ExecutionException {

        List<ModelClass> mList = new ArrayList<>();

        for (final ModelClass modelClass : inputs) {
            ReplacementModel mc = (ReplacementModel) modelClass;

            mc.setDescriptor(getModelDescriptor());
            mc.setProperty(setScale(modelProperty));
            mc.setSelector(randGroupObservation(modelMatrix.getRowCount()));
            mc.calcModel();
            mList.add(mc);

            changes.firePropertyChange(ModelStatic.EVENT_MODEL_CALCULATED, null, null);
        }
        return mList;
    }

    public List<ModelClass> processInputsParallel(List<ModelClass> inputs)
            throws InterruptedException, ExecutionException {

        int threads = Runtime.getRuntime().availableProcessors();
        ExecutorService service = Executors.newFixedThreadPool(threads);

        List<Future<ModelClass>> futures = new ArrayList<>();

        for (final ModelClass modelClass : inputs) {
            Callable<ModelClass> callable = new Callable<ModelClass>() {

                @Override
                public ModelClass call() throws Exception {

                    ReplacementModel mc = (ReplacementModel) modelClass;
                    mc.setDescriptor(getModelDescriptor());
                    mc.setProperty(setScale(modelProperty));
                    mc.setSelector(randGroupObservation(modelMatrix.getRowCount()));
                    mc.calcModel();
                    return mc;
                }
            };
            futures.add(service.submit(callable));
        }
        service.shutdown();

        List<ModelClass> outputs = new ArrayList<>();
        for (Future<ModelClass> future : futures) {
            outputs.add(future.get());
        }
        return outputs;
    }

    public boolean calcModel() {

        Thread t2 = new Thread(new Runnable() {

            @Override
            public void run() {

                try {

                    long startTime = System.currentTimeMillis();
                    boolean suc = processInputsSingle();
                    changes.firePropertyChange(ModelStatic.EVENT_MODEL_CALCULATED, null, null);

                    System.out.println(ColorStatic.BLUE
                            + "Calculating regression model with " + modelList.getModelList().size()
                            + " sub models." + ColorStatic.RESET);

                    long stopTime = System.currentTimeMillis();
                    long elapsedTime = stopTime - startTime;
                    System.out.println(ColorStatic.BLUE + "Calculation time was " + elapsedTime / 1000 + " sec." + ColorStatic.RESET);
                    logger.info("Models calculated in {} sec", elapsedTime / 1000);

                    System.out.println("ModelGroup thinks calculation finished.");
                    changes.firePropertyChange(ModelStatic.EVENT_MODEL_CALCULATED, null, null);

                } catch (InterruptedException | ExecutionException ex) {
                    logger.error(ex.getMessage());
                    System.out.println(ColorStatic.RED + ex.getMessage() + ColorStatic.RESET);
                }
            }
        });
        t2.start();

        return true;
    }

    public Object getModel(int i) {
        return this.modelList.getModel(i);

    }

    public void setModel(int i, ModelClass model) {
        this.modelList.setModel(i, model);
    }

    public void addModel(ModelClass model) {
        this.modelList.addModel(model);
        int n = this.modelList.getModelList().size();
    }

    public Matrix getModelDescriptor() {

        Matrix m = this.modelMatrix.denormalize(normalizationFactor);
        return m;

    }

    public void setModelDescriptor(Matrix mM) {

        if (normalize == true) {
            this.normalizationFactor = mM.getNormalizationFactor();
            this.modelMatrix = mM.normalize(normalizationFactor);
        } else {
            this.modelMatrix = mM;
        }
    }

    public Matrix getPredictDescriptor() {

        Matrix m = this.predictMatrix.denormalize(normalizationFactor);
        return m;
    }

    public void setPredictDescriptor(Matrix pM) {

        if (normalize == true) {
            Matrix m = pM.normalize(normalizationFactor);
            this.predictMatrix = m;
        } else {
            this.predictMatrix = pM;
        }
    }

    private List<Double> setScale(List<Double> valueList) {
        List<Double> logValueList = new ArrayList<>(valueList.size());

        if (this.logScale == true) {
            for (int i = 0; i < valueList.size(); i++) {
                logValueList.add((Math.log(valueList.get(i)) / Math.log(2)));
            }
            return logValueList;
        } else {
            return valueList;
        }
    }

    private List<Double> reverseScale(List<Double> valueList) {

        List<Double> normValueList = new ArrayList<>(valueList.size());

        if (this.logScale == true) {
            for (int i = 0; i < valueList.size(); i++) {
                normValueList.add(Math.pow(2, valueList.get(i)));
            }
            return normValueList;

        } else {
            return valueList;
        }
    }

    public List<Double> getModelProperty() {

        return this.modelProperty;
    }

    public List<List<Double>> getModelPropertyDetail(Double modelPar, Double predPar) {

        List<List<Double>> modelPropDetail = new ArrayList<>();
        for (ModelClass modelList1 : this.modelList.getModelList()) {

            ReplacementModel rm = (ReplacementModel) modelList1;
            if (rm.modelRSquare() > modelPar && rm.predictiveRSquare() > predPar) {

                modelPropDetail.add(this.reverseScale(this.getModelProperty()));
            }
        }
        return modelPropDetail;
    }

    public List<List<Double>> getModelPropertyDetail() {

        List<List<Double>> modelPropDetail = new ArrayList<>();

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            modelPropDetail.add(this.reverseScale(this.getModelProperty()));
        }
        return modelPropDetail;
    }

    public List<List<Double>> getPredictPropertyDetail(Double modelPar, Double predPar) {
        List<List<Double>> predictPropDetail = new ArrayList<>();
        for (ModelClass modelList1 : this.modelList.getModelList()) {

            ReplacementModel rm = (ReplacementModel) modelList1;
            if (rm.modelRSquare() > modelPar && rm.predictiveRSquare() > predPar) {
                predictPropDetail.add(this.reverseScale(this.getPredictProperty()));
            }
        }
        return predictPropDetail;
    }

    public List<List<Double>> getPredictPropertyDetail() {
        List<List<Double>> predictPropDetail = new ArrayList<>();
        for (ModelClass modelList1 : this.modelList.getModelList()) {
            predictPropDetail.add(this.reverseScale(this.getPredictProperty()));
        }
        return predictPropDetail;
    }

    public List<List<Double>> getModelResponseDetail(Double modelPar, Double predPar) {

        List<List<Double>> modelRespDetail = new ArrayList<>();

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;
            if (rm.modelRSquare() > modelPar && rm.predictiveRSquare() > predPar) {

                modelRespDetail.add(this.reverseScale(rm.getModelResponse(this.getModelDescriptor())));
            }
        }
        return modelRespDetail;

    }

    public List<List<Double>> getModelResponseDetail() {

        List<List<Double>> modelRespDetail = new ArrayList<>();

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;
            modelRespDetail.add(this.reverseScale(rm.getModelResponse(this.getModelDescriptor())));
        }
        return modelRespDetail;

    }

    public List<List<Double>> getPredictResponseDetail(Double modelPar, Double predPar) {
        List<List<Double>> predRespDetail = new ArrayList<>();

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;
            if (rm.modelRSquare() > modelPar && rm.predictiveRSquare() > predPar) {
                predRespDetail.add(this.reverseScale(rm.getModelResponse(this.getPredictDescriptor())));
            }
        }
        return predRespDetail;
    }

    public List<List<Double>> getPredictResponseDetail() {
        List<List<Double>> predRespDetail = new ArrayList<>();

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;
            predRespDetail.add(this.reverseScale(rm.getModelResponse(this.getPredictDescriptor())));
        }
        return predRespDetail;
    }

    private List<Double> getMedian(List<List<Double>> array) {
        List<Double> respMedian = new ArrayList<>();

        int n = array.get(0).size();

        for (int i = 0; i < n; i++) {
            List<Double> obs = new ArrayList<>();
            for (List<Double> array1 : array) {
                if (array1.get(i) != null && !Double.isNaN(array1.get(i)) && !Double.isInfinite(array1.get(i))) {
                    obs.add(array1.get(i));
                }
            }
            Collections.sort(obs);

            int middle = obs.size() / 2;
            if (obs.size() % 2 == 1) {
                respMedian.add(obs.get(middle));
            } else {
                respMedian.add((obs.get(middle - 1) + obs.get(middle)) / 2.0);
            }
        }
        return respMedian;
    }

    private List<Double> getMean(List<List<Double>> array) {

        if (array == null) {
            return null;
        }

        List<Double> respMean = new ArrayList<>();

        List<Double> pRM = new ArrayList<>();
        List<Integer> pRc = new ArrayList<>();

        for (List<Double> array1 : array) {
            for (int j = 0; j < array1.size(); j++) {

                if (array1.get(j) == null) {
                    continue;
                }
                if (Double.isNaN(array1.get(j))) {
                    continue;
                }

                try {
                    pRM.set(j, pRM.get(j) + array1.get(j));
                    pRc.set(j, pRc.get(j) + 1);
                } catch (IndexOutOfBoundsException e) {
                    pRM.add(array1.get(j));
                    pRc.add(1);
                }
            }
        }

        for (int i = 0; i < pRM.size(); i++) {
            respMean.add(pRM.get(i) / pRc.get(i));
        }

        return respMean;
    }

    public List<Double> getModelDescriptorValue() {

        int n = this.getModelDescriptor().getColumnCount();
        coefVal = new ArrayList<>(n);

        for (int i = 0; i < n; i++) {
            coefVal.add(0.0);
        }

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;

            if (rm.modelRSquare() > this.lowerModelRSquare && rm.predictiveRSquare() > this.lowerPredictRSquare) {
                try {
                    List<Double> coef = rm.getCoefficients();
                    List<Integer> vec = rm.getDescriptorsSelected();

                    for (int i = 0; i < vec.size(); i++) {
                        int ni = vec.get(i);
                        coefVal.set(ni, coefVal.get(ni) + coef.get(i + 1));
                    }
                } catch (Exception e) {
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }

            }
        }
        System.out.println(ColorStatic.PURPLE + coefVal.toString() + ColorStatic.RESET);
        return coefVal;
    }

    public List<Integer> getModelDescriptorPos() {

        int n = this.getModelDescriptor().getColumnCount();
        coefCount = new ArrayList<>(n);

        for (int i = 0; i < n; i++) {
            coefCount.add(0);
        }

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;

            if (rm.modelRSquare() > this.lowerModelRSquare && rm.predictiveRSquare() > this.lowerPredictRSquare) {
                try {
                    List<Integer> vec = rm.getDescriptorsSelected();
                    for (Integer vec1 : vec) {
                        int ni = vec1;
                        coefCount.set(ni, coefCount.get(ni) + 1);
                    }
                } catch (Exception e) {
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }
            }
        }

        System.out.println(ColorStatic.PURPLE + coefCount.toString() + ColorStatic.RESET);
        return coefCount;
    }

    public List<Double> getModelResponse() {

        List<List<Double>> modelRespDetail = new ArrayList<>();
        int i = 0;
        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;

            List<Double> n = rm.getCoefficients();
            if (n == null) {
                return null;
            }

            System.out.println("vecFinal from ModelGroup getModelResponse() " + n.get(0));

            if (rm.modelRSquare() != null) {
                if (rm.modelRSquare() > this.lowerModelRSquare && rm.predictiveRSquare() > this.lowerPredictRSquare) {
                    modelRespDetail.add(this.reverseScale(rm.getModelResponse(this.getModelDescriptor())));

                    System.out.println(ColorStatic.BLUE + "Added model with r2 "
                            + rm.modelRSquare() + " and "
                            + rm.predictiveRSquare() + ColorStatic.RESET);
                } else {
                    System.out.println("Ignored model with r2 "
                            + rm.modelRSquare() + " and "
                            + rm.predictiveRSquare());
                }
            }
        }

        if (modelRespDetail.isEmpty()) {
            System.out.println("No model found.");
            return null;
        }

        List<Double> medianResponse = getMedian(modelRespDetail);

        return medianResponse;
    }

    public void setModelProperty(List<Double> modProp) {
        this.modelProperty = modProp;
    }

    public List<Double> getPredictProperty() {
        return this.predictProperty;
    }

    public void setPredictProperty(List<Double> pProp) {

        this.predictProperty = pProp;
    }

    /*
     returns a List of prediction results, this result is a mean of single or sub model predictions 
     */
    public List<Double> getPredictResponse() {
        List<List<Double>> predRespDetail = new ArrayList<>();

        for (ModelClass modelList1 : this.modelList.getModelList()) {
            ReplacementModel rm = (ReplacementModel) modelList1;

            Double mrsquare = rm.modelRSquare();
            if (mrsquare != null) {

                if (mrsquare > this.lowerModelRSquare && rm.predictiveRSquare() > this.lowerPredictRSquare) {
                    predRespDetail.add(this.reverseScale(rm.getModelResponse(this.getPredictDescriptor())));
                    System.out.println(ColorStatic.GREEN
                            + "Added model with r2 " + rm.modelRSquare()
                            + "; predictive r2 " + rm.predictiveRSquare() + ColorStatic.RESET);
                }
            }
        }
        if (predRespDetail.isEmpty()) {
            System.out.println("No models found fitting specifications of model "
                    + "rsquare > " + this.lowerModelRSquare
                    + "; predictive rsquare > " + this.lowerPredictRSquare
                    + ".\n Adapt model fitting specifications or dataset.");
            return null;
        }

        List<Double> medianResponse = getMedian(predRespDetail);
        return medianResponse;
    }

    public void addPredictResponse(List<Double> singleResponse) {
        List<List<Double>> oldPredictResponse = null;

        if (this.predictResponse != null) {
            oldPredictResponse = new ArrayList<>(this.predictResponse);
        }
        this.predictResponse.add(singleResponse);
    }

    public void addModelResponse(List<Double> singleResponse) {
        List<List<Double>> oldModelResponse = null;

        if (this.modelResponse != null) {
            oldModelResponse = new ArrayList<>(this.modelResponse);
        }
        this.modelResponse.add(singleResponse);
    }

    public void setPredictResponse(int i, List<Double> singleResponse) {
        this.predictResponse.set(i, singleResponse);
    }

    public void setPredictResponse(List<List<Double>> predictResponse) {
        this.predictResponse = predictResponse;
    }

    public Double getRSquare() {

        Double rsquare = null;

        try {
            List<Double> yPred = this.getModelResponse();
            List<Double> y = this.getModelProperty();
            Double sumA = 0.0;
            Double sumB = 0.0;
            Double meanY = MathUtil.mean(y);

            for (int i = 0; i < yPred.size(); i++) {

                sumA = sumA + Math.pow((yPred.get(i) - y.get(i)), 2);
                sumB = sumB + Math.pow((y.get(i) - meanY), 2);

            }
            rsquare = 1 - (sumA / sumB);
        } catch (Exception e) {

            logger.error(e.getMessage());
        }

        return rsquare;
    }

    public Double getPredRSquare() {

        Double predrsquare = 0.00;
        try {
            List<Double> yPred = this.getPredictResponse();
            List<Double> y = this.getPredictProperty();
            Double sumA = 0.0;
            Double sumB = 0.0;
            Double meanY = MathUtil.mean(this.getModelProperty());

            for (int i = 0; i < yPred.size(); i++) {
                sumA = sumA + Math.pow((yPred.get(i) - y.get(i)), 2);
                sumB = sumB + Math.pow((y.get(i) - meanY), 2);
            }
            predrsquare = 1 - (sumA / sumB);
        } catch (Exception e) {
            logger.error(e.getMessage());
        }

        return predrsquare;
    }

    public void displayModelData() {

        List<Double> modelProp = this.getModelProperty();
        List<Double> modelPred = this.getModelResponse();

        System.out.println(ColorStatic.PURPLE + "model data");
        for (int i = 0; i < modelProp.size(); i++) {

            System.out.println(ColorStatic.PURPLE + modelProp.get(i)
                    + "  " + modelPred.get(i) + "; ..." + ColorStatic.RESET);
        }

        List<Double> predProp = this.getPredictProperty();
        List<Double> predPred = this.getPredictResponse();

        System.out.println(ColorStatic.PURPLE + "predict data");
        for (int i = 0; i < predProp.size(); i++) {

            System.out.println(ColorStatic.PURPLE + predProp.get(i)
                    + "  " + predPred.get(i) + "; ..." + ColorStatic.RESET);
        }

        List<Integer> descrPos = this.getModelDescriptorPos();
        List<Double> descrVal = this.getModelDescriptorValue();

        System.out.println(ColorStatic.PURPLE + "descriptor data");
        for (int i = 0; i < descrPos.size(); i++) {

            System.out.println(ColorStatic.BLACK + descrPos.get(i)
                    + "  " + descrVal.get(i) + "; ..." + ColorStatic.RESET);
        }
    }

    public ModelList getModelList() {
        return modelList;
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        changes.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        changes.removePropertyChangeListener(listener);
    }

}
