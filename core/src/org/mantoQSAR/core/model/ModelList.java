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
import java.util.List;


public class ModelList {
    
    List<ModelClass> modelList;
    final private PropertyChangeSupport changes = new PropertyChangeSupport(this);
    
    ModelList(){
        modelList = new ArrayList<>(); 
    }

    public List<ModelClass> getModelList() {
        return modelList;
    }

    public void setModelList(List<ModelClass> modelList) {
        this.modelList = modelList;
    }
    
    public void addModel(ModelClass model){
        modelList.add(model);
        changes.firePropertyChange(ModelStatic.EVENT_MODEL_CHANGED, null, null);
    }
    
    public void setModel(int pos, ModelClass model){
        modelList.set(pos, model);
        changes.firePropertyChange(ModelStatic.EVENT_MODEL_CHANGED, null, null);
    }
    
    public ModelClass getModel(int pos){
        return modelList.get(pos);
    }
    
    
    public void addPropertyChangeListener(PropertyChangeListener listener) {
        changes.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        changes.removePropertyChangeListener(listener);
    }

}
