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

package org.mantoQSAR.gui.resources;

import java.awt.Color;


public class ColorRange {
    
    private double maxPos = 0; 
    private double min = 0;
    private double maxNeg = 0; 
    private boolean sigmoidal = false; 
    private int scheme = 1; // sets color scheme to use
    
    public ColorRange(){
        
        
    }

    public Color getColor(Double value) {

        value = this.normalize(value);
        float hue = 1.0f; 
        float b = 0.9f;  // brightness
        float s = 0.9f; // saturation; 
        
        switch (scheme) {
            case 1:// color range green yellow red
                hue = (float)(this.normalize(value)+1.0)* (60/360);
                
                break;

            case 2: // color range black - red - yellow
                
                if(value < 0.50 && value >= 0.0){
                    hue = 10; 
                    b = (float) (value*2.0);

                }
                if(value >= 0.50){
                    hue = (float)(value-0.5)*2*(45/360); 
                    b = 1; 
                }
                
                 if(value <= 0.0){
                    hue = (210/360); 
                    b = (float) (value*-1.0);
                }
                break; 
            
            case 3: 
                if(value > 0.0){
                hue = (float) 0.65;
                s = (float) Math.abs(value);
                }else{
                hue = (float) (value*0.1);
                s = (float) Math.abs(value);
                }
                break;
                
            case 4:     
               hue = (float) 0.65; 
                break; 
            default:
                
                System.err.println("No viable color scheme found.");
                break;
        }
        return Color.getHSBColor(hue, s, b); 
    }
    
    private double normalize(double value) {
        
      double val = 0;    
      
      if(this.sigmoidal != true){
          
          if (value > 0) {
              val = value / this.maxPos;

              if (val > this.maxPos) {
                  val = this.maxPos;
              }
          } else {
              val = value / this.maxNeg;

              if (val > this.maxNeg) {
                  val = this.maxNeg;
              }

          }
          
      }else{
          // sigmoidal map of -infinite to +infinite to range 0 to 1
          double v0 = value / Math.sqrt(1 + value * value); 
         val = v0; 

      }
      
    return val; 
}
    
    public double getMaxPos() {
        return this.maxPos;
    }

    public void setMaxPos(double max) {
        this.maxPos = max;
    }
    
    public void setMaxMin(double max){
        this.maxNeg = max; 
    }

    public double getMaxMin(){
        return this.maxNeg; 
    }
    
    public double getMin() {
        return min;
    }

    public void setMin(double min) {
        this.min = min;
    }

    public boolean isSigmoidal() {
        return sigmoidal;
    }

    public void setSigmoidal(boolean sigmoidal) {
        this.sigmoidal = sigmoidal;
    }

    public int getScheme() {
        return scheme;
    }

    public void setScheme(int scheme) {
        this.scheme = scheme;
    }
    
    
    
    
    
}
