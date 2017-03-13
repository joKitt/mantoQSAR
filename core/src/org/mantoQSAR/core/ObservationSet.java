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


public class ObservationSet {
    
   private String name; 
   private String species; 
   private Boolean rcsb; 
   private String pdbID; 
   private String file; 
   private Condition condition; 
   private Double response = null; 
   private int count; 
   private String note; 
   private boolean active; // active observations are considered in model generation and prediction
   private boolean predict; // active observations are to be predicted; 
    
    public ObservationSet(){
        this.name = "empty";
        this.species = "empty";
        this.rcsb = false; 
        this.file = "empty";
        this.condition = new Condition(); 
        this.count = 0; 
        this.active = true; 
        this.predict = false; 
    }

    public void setActive(boolean active) {
        this.active = active;
    }

    public void setPredict(boolean predict) {
        this.predict = predict; 
    }

    public class Condition{
        
        public Double pH; 
        public Double ionicStrength; 
        public Double concentration; 
        
        public Condition(){
        }

        public Double getpH() {
            return pH;
        }

        public void setpH(Double pH) {
            this.pH = pH;
        }

        public Double getIonicStrength() {
            return ionicStrength;
        }

        public void setIonicStrength(Double ionicStrength) {
            this.ionicStrength = ionicStrength;
        }

        public Double getConcentration() {
            return concentration;
        }

        public void setConcentration(Double concentration) {
            this.concentration = concentration;
        }
        
        
        
    }

    public boolean isActive() {
        if(active == true & this.response == null){
           return false; 
        }else{
           return active; 
        }
    }

    public boolean isPredict() {
        return predict;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public Boolean getRcsb() {
        return rcsb;
    }

    public void setRcsb(Boolean rcsb) {
        this.rcsb = rcsb;
    }

    public String getPdbID() {
        return pdbID;
    }

    public void setPdbID(String pdbID) {
        this.pdbID = pdbID;
    }

    public String getFile() {
        return file;
    }

    public void setFile(String file) {
        this.file = file;
    }

    public Condition getCondition() {
        return condition;
    }

    public void setCondition(Condition condition) {
        this.condition = condition;
    }

    public Double getResponse() {
        return response;
    }

    public void setResponse(Double response) {
        this.response = response;
    }

    public String getNote() {
        return note;
    }

    public void setNote(String note) {
        this.note = note;
    }
}
