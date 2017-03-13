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


public class ModelStatic {
    public static final String ACTION_CALC_MODEL = "action_calc_model";
    public static final String ACTION_TEST_MODEL = "action_test_model";
    public static final String ACTION_TEST_MATRIX = "action_test_matrix";
    public static final String ACTION_CALC_STEPWISE = "action_calc_stepwise";
    public static final String ACTION_CALC_PREDICTION = "action_calc_prediction";
    public static final String ACTION_RAND_PREDICTION = "action_rand_prediction";
    
    public static final String CHANGE_PREDICT_PROPERTY = "change_predictProperty";
    public static final String CHANGE_PREDICT_RESPONSE = "change_predictResponse";
    public static final String CHANGE_PREDICT_DESCRIPTOR = "change_predictDescriptor";
    
    public static final String CHANGE_MODEL_PROPERTY = "change_modelProperty";
    public static final String CHANGE_MODEL_RESPONSE = "change_modelResponse";
    public static final String CHANGE_MODEL_DESCRIPTOR = "change_modelDescriptor";
    
    public static final String EVENT_MODEL_CALCULATED = "event_modelCalculated";
    public static final String EVENT_MODEL_CHANGED = "event_modelChanged";
    
}
