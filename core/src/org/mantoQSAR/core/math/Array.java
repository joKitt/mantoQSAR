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


package org.mantoQSAR.core.math;

public abstract class Array {
    
    public static int minIndex(Integer[] array){
        int minIndex = 0;
        
        for (int i = 0; i < array.length; i++){
            if ((array[i] < array[minIndex])){
            minIndex = i;
        }
        }
     return minIndex;    
    }
    
    public static int minIndex(Double[] array){
        int minIndex = 0;
        
        for (int i = 0; i < array.length; i++){
            if ((array[i] < array[minIndex])){
            minIndex = i;
        }
        }
     return minIndex;    
    }
    
    
    public static int maxIndex(Integer[] array){
        int maxIndex = 0;
        
        for (int i = 0; i < array.length; i++){
            if ((array[i] > array[maxIndex])){
            maxIndex = i;
        }
        }
     return maxIndex;    
    }
    
    public static int maxIndex(Double[] array){
        int maxIndex = 0;
        
        for (int i = 0; i < array.length; i++){
            if ((array[i] > array[maxIndex])){
            maxIndex = i;
        }
        }
     return maxIndex;    
    }
    
    public static double rms(Double[] nums){
        double ms = 0;
        for (Double num : nums) {
            ms += num * num;
        }
        ms /= nums.length;
        return Math.sqrt(ms);
    }
}
