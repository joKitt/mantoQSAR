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

import org.mantoQSAR.core.model.ModelStatic;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.model.*;
import org.mantoQSAR.gui.molecule.MoleculeController;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public final class ModelController implements ActionListener{
    
    MoleculeController moleculeController;     
    Screen screen; 
    final private PropertyChangeSupport changes = new PropertyChangeSupport(this);

    private static ModelController instance = null; 
    ReplacementModel model; 
    ModelGroup modelGroup;
    ModelList modelList; 
    
    Logger logger; 
    
   protected ModelController(){
       
      this.initModelController();
       
   } 
    
   public static ModelController getInstance(){
       
       if(instance == null){
            instance = new ModelController(); 
        }
        return instance; 
   } 
   
   public void initModelController(){
       
       logger = LoggerFactory.getLogger(ModelController.class);
       logger.trace("instantiated");
       logger.trace("moleculeController instantiated");
       moleculeController = MoleculeController.getInstance();

       // set value for development purpose
       modelGroup = new ModelGroup(250);
       modelList = modelGroup.getModelList();
       screen = moleculeController.getScreen();      
       initEventListener();
   }

    @Override
    public void actionPerformed(ActionEvent e) {
         logger.debug("Action call: " + e.getActionCommand());

        if(e.getActionCommand().equalsIgnoreCase(ModelStatic.ACTION_TEST_MODEL)){
        } 
        if(e.getActionCommand().equalsIgnoreCase(ModelStatic.ACTION_TEST_MATRIX)){
            model.testMatrix();
        }   
        
        if(e.getActionCommand().equalsIgnoreCase(ModelStatic.ACTION_CALC_MODEL)){

                    modelGroup.setModelDescriptor(screen.getModelDescriptorMatrix());
                    modelGroup.setModelProperty(screen.getModelObservationValues());
                    modelGroup.setPredictDescriptor(screen.getPredictDescriptorMatrix());
                    modelGroup.setPredictProperty(screen.getPredictObservationValues());

                    modelGroup.calcModel();

        }
        
        if(e.getActionCommand().equalsIgnoreCase(ModelStatic.ACTION_CALC_STEPWISE)){

            ReplacementModel m = (ReplacementModel) modelGroup.getModel(0);
            m.calcStepwise(8);
            System.out.println("Stepwise model calculation finished");
        }
        
        if(e.getActionCommand().equalsIgnoreCase(ModelStatic.ACTION_CALC_PREDICTION)){

        } 
   
        
    }

    public ModelGroup getModelGroup() {
        return modelGroup;
    }

    public void setModelGroup(ModelGroup modelGroup) {
        this.modelGroup = modelGroup;
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        changes.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        changes.removePropertyChangeListener(listener);
    }

    private void initEventListener() {

        modelGroup.addPropertyChangeListener(new PropertyChangeListener() {

            @Override
            public void propertyChange(PropertyChangeEvent evt) {
            
                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CALCULATED)) {

                    changes.firePropertyChange(ModelStatic.EVENT_MODEL_CALCULATED, null, null);
                }

                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CHANGED)) {
                    changes.firePropertyChange(ModelStatic.EVENT_MODEL_CHANGED, null, null);

                }

            }
        });


        modelList.addPropertyChangeListener(new PropertyChangeListener() {

            @Override
            public void propertyChange(PropertyChangeEvent evt) {
             
                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CALCULATED)) {

                    changes.firePropertyChange(ModelStatic.EVENT_MODEL_CALCULATED, null, null);
                }

                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CHANGED)) {
                    changes.firePropertyChange(ModelStatic.EVENT_MODEL_CHANGED, null, null);
                }
            }
        });
    }
}
