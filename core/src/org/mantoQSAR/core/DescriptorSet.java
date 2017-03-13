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

package org.mantoQSAR.core;


public class DescriptorSet {
    
    String name = "empty"; 
    Surface surface; 
    Projection projection;
    DescriptorProperty descriptor; 
    
public void DescriptorSet(){
    
}
    @Override
    public String toString(){
        return this.name + "\n" + descriptor.toString();
}

    public String getName() {
        return name;
    }

    public Surface getSurface() {
        return surface;
    }

    public Projection getProjection() {
        return projection;
    }

    public DescriptorProperty getDescriptor() {
        return descriptor;
    }

    public class Surface{
     float resolution;  
     String typeofsurface; 
     Double sizeofsphere; 
     String property;
     int mapfunc; 

        public Surface() {
            this.resolution = 0.1f;
            this.typeofsurface = "unknown";
            this.sizeofsphere = 1.1;
            this.property = "esp";
            this.mapfunc = 0; 
            
        }

        public String getProperty() {
            return property;
        }

        public void setProperty(String property) {
            this.property = property;
        }

        public int getMapFunction() {
            return mapfunc;
        }

        public void setMapFunction(int mapfunc) {
            this.mapfunc = mapfunc;
        }

        public float getResolution() {
            return resolution;
        }

        public void setResolution(float resolution) {
            this.resolution = resolution;
        }

        public String getTypeOfSurface() {
            return typeofsurface;
        }

        public void setTypeOfSurface(String typeOfSurface) {
            this.typeofsurface = typeOfSurface;
        }

        public Double getSizeOfSphere() {
            return this.sizeofsphere;
        }

        public void setSizeOfSphere(Double sizeOfSphere) {
            this.sizeofsphere = sizeOfSphere;
        }
        
       
        
}

    public class Projection{
       public String type; 
       public Double density; 
       public Double size;
       public Double distance; 
       public int mapfunc; 
       public int orientation; 
       public String select; 
       public String selectFunc; 
       public int selectIO; 
       
        public Projection(){
            this.selectIO = 1;

        }

        public int getSelectIO() {
            return selectIO;
        }
        
        public Double getDensity() {
            return density;
        }

        public void setDensity(Double density) {
            this.density = density;
        }

        public void setSelectIO(int selectIO) {
            this.selectIO = selectIO;
        }
        public void setSelectIO(boolean selectIO) {
            if(selectIO == true){
                this.selectIO = 1; 
            } else{
                this.selectIO = 0;
            }
        }

        public String getType() {
            return type;
        }

        public void setType(String type) {
            this.type = type;
        }

        public Double getSize() {
            return size;
        }

        public void setSize(Double size) {
            this.size = size;
        }

        public Double getDistance() {
            return distance;
        }

        public void setDistance(Double distance) {
            this.distance = distance;
        }

        public int getMapFunction() {
            return mapfunc;
        }

        public void setMapFunction(int mapFunction) {
            this.mapfunc = mapFunction;
        }

        public String getSelectID() {
            return select;
        }

        public void setSelectID(String selectID) {
            this.select = selectID;
        }

        public String getSelectFunction() {
            return selectFunc;
        }

        public void setSelectFunction(String selectFunction) {
            this.selectFunc = selectFunction;
        }

        public int getMapfunc() {
            return mapfunc;
        }

        public void setMapfunc(int mapfunc) {
            this.mapfunc = mapfunc;
        }

        public int getOrientation() {
            return orientation;
        }

        public void setOrientation(int orientation) {
            this.orientation = orientation;
        }

        public String getSelect() {
            return select;
        }

        public void setSelect(String select) {
            this.select = select;
        }

        public String getSelectFunc() {
            return selectFunc;
        }

        public void setSelectFunc(String selectFunc) {
            this.selectFunc = selectFunc;
        }
        
        
    }

    public class DescriptorProperty{
    public String name; 
    public String property; 
    public Double binScale; 
    
    @Override
    public String toString(){
        String out = "descriptor property " + name + " \n " + 
                     "mapping " + property + "\n" +
                     "scale " + binScale + "\n";
    return out; 
}
}
}

