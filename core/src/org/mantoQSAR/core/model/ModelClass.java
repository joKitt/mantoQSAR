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

import java.beans.PropertyChangeSupport;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang3.NotImplementedException;
import org.mantoQSAR.core.math.Matrix;
import org.mantoQSAR.core.util.ColorStatic;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class ModelClass {

    private final PropertyChangeSupport changes = new PropertyChangeSupport(this);

    Logger logger;
    List<Double> property = null;
    Matrix matrix = null;

    public List<Integer> entryPosition = null;

    ModelClass() {
        this.property = new ArrayList<>();
        logger = LoggerFactory.getLogger(ModelClass.class);

    }

    public void calcModel() {
        throw new NotImplementedException("This function is overwritten by children class");
    }

    public void setProperty(List<Double> property) {
        List<Double> oldP = null;
        if (this.property != null) {
            oldP = new ArrayList<>(this.property);
        }
        this.property = property;
    }

    public List<Double> getModelProperty() {
        if (this.entryPosition == null) {
            logger.error("entryPosition not defined prior to getModelProperty()");
            return null;
        }
        if (this.entryPosition.size() != this.property.size()) {
            logger.error("entryPosition has false dimensions "
                    + this.entryPosition.size() + "; "
                    + this.property.size() + " expected");
            return null;
        }

        List<Double> modelP = new ArrayList<>();

        for (int i = 0; i < this.entryPosition.size(); i++) {
            if (this.entryPosition.get(i) == 1) {
                modelP.add(this.property.get(i));
            }
        }
        return modelP;
    }

    public Matrix getModelDescriptor() {

        if (this.entryPosition == null) {
            return null;
        }

        if (this.entryPosition.size() != matrix.getRowCount()) {
            System.out.println(ColorStatic.RED + "entryPosition List with wrong size"
                    + this.entryPosition.size() + " / " + this.matrix.getRowCount() + ColorStatic.RESET);
            return null;
        }

        int mn = 0;
        for (Integer entryPosition1 : this.entryPosition) {
            if (entryPosition1 == 1) {
                mn++;
            }
        }

        Matrix rMatrix = new Matrix(mn, matrix.getColumnCount());

        int mn2 = 0;

        for (int i = 0; i < this.entryPosition.size(); i++) {
            if (this.entryPosition.get(i) == 1) {
                for (int j = 0; j < this.matrix.getColumnCount(); j++) {
                    rMatrix.data[mn2][j] = this.matrix.data[i][j];
                }
                mn2++;
            }
        }
        return rMatrix;
    }

    public Matrix getPredictDescriptor() {

        int mn = 0;
        for (Integer entryPosition1 : this.entryPosition) {
            if (entryPosition1 == 2) {
                mn++;
            }
        }

        Matrix rMatrix = new Matrix(mn, matrix.getColumnCount());

        int mn2 = 0;

        for (int i = 0; i < this.entryPosition.size(); i++) {
            if (this.entryPosition.get(i) == 2) {
                for (int j = 0; j < matrix.getColumnCount(); j++) {
                    rMatrix.data[mn2][j] = this.matrix.data[i][j];
                }
                mn2++;
            }
        }
        return rMatrix;
    }

    public void setSelector(List<Integer> selector) {
        this.entryPosition = selector;
    }

    public void setDescriptor(Matrix matrix) {
        this.matrix = matrix;
    }

    public List<Double> getPredictProperty() {

        if (this.entryPosition != null
                && this.entryPosition.size() == matrix.getRowCount()) {

            List<Double> predP = new ArrayList<>();

            for (int i = 0; i < this.entryPosition.size(); i++) {
                if (this.entryPosition.get(i) == 2) {
                    predP.add(this.property.get(i));
                }
            }
            return predP;
        }
        return null;
    }
}
