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

import org.biojava.nbio.structure.jama.Maths;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class Matrix {
    
    public double[][] data; 
    private final int m; 
    private final int n; 
    
    private double[] s;
    private double[][] U, V; 
    
    Logger logger; 
    
    public Matrix(int m, int n) {
        
        logger = LoggerFactory.getLogger(Matrix.class);
        this.m = m;
        this.n = n;
        this.data = new double[m][n];
        
        
    }

    public Matrix(double[][] data) {
        
        logger = LoggerFactory.getLogger(Matrix.class);
        
        this.m = data.length;
        this.n = data[0].length;
        this.data = new double[m][n];
        for (int i = 0; i < this.m; i++)
            for (int j = 0; j < this.n; j++)
                    this.data[i][j] = data[i][j];

    }
    
    private void calcParameter(){
        
         int nu = Math.min(m,n);
      s = new double [Math.min(m+1,n)];
      U = new double [m][nu];
      V = new double [n][n];
      double[] e = new double [n];
      double[] work = new double [m];
      boolean wantu = true;
      boolean wantv = true; 
      
      int nct = Math.min(m-1,n);
      int nrt = Math.max(0,Math.min(n-2,m));
      
      for (int k = 0; k < Math.max(nct,nrt); k++) {
          if (k < nct) {
           s[k] = 0;
            for (int i = k; i < m; i++) {
               s[k] = Maths.hypot(s[k],this.data[i][k]);
            }
            if (s[k] != 0.0) {
               if (this.data[k][k] < 0.0) {
                  s[k] = -s[k];
               }
               for (int i = k; i < m; i++) {
                  this.data[i][k] /= s[k];
               }
               this.data[k][k] += 1.0;
            }
            s[k] = -s[k];  
              
          }
          for (int j = k+1; j < n; j++) {
              if ((k < nct) & (s[k] != 0.0))  {

               double t = 0;
               for (int i = k; i < m; i++) {
                  t += this.data[i][k]*this.data[i][j];
               }
               t = -t/this.data[k][k];
               for (int i = k; i < m; i++) {
                  this.data[i][j] += t*this.data[i][k];
               }    
              }

            e[j] = this.data[k][j];
          }
          
         if (wantu & (k < nct)) {

            for (int i = k; i < m; i++) {
               U[i][k] = this.data[i][k];
            }
         } 
         if (k < nrt) {
            e[k] = 0;
            for (int i = k+1; i < n; i++) {
               e[k] = Maths.hypot(e[k],e[i]);
            }
            if (e[k] != 0.0) {
               if (e[k+1] < 0.0) {
                  e[k] = -e[k];
               }
               for (int i = k+1; i < n; i++) {
                  e[i] /= e[k];
               }
               e[k+1] += 1.0;
            } 
            e[k] = -e[k];
            if ((k+1 < m) & (e[k] != 0.0)) {
               for (int i = k+1; i < m; i++) {
                  work[i] = 0.0;
               }
               for (int j = k+1; j < n; j++) {
                  for (int i = k+1; i < m; i++) {
                     work[i] += e[j]*this.data[i][j];
                  }
               }
               for (int j = k+1; j < n; j++) {
                  double t = -e[j]/e[k+1];
                  for (int i = k+1; i < m; i++) {
                     this.data[i][j] += t*work[i];
                  }
               }
            }
            if (wantv) {
               for (int i = k+1; i < n; i++) {
                  V[i][k] = e[i];
               } 
            }
         }
      }

      int p = Math.min(n,m+1);
      if (nct < n) {
         s[nct] = this.data[nct][nct];
      }
      if (m < p) {
         s[p-1] = 0.0;
      }
      if (nrt+1 < p) {
         e[nrt] = this.data[nrt][p-1];
      }
      e[p-1] = 0.0;
      
      if (wantu) {
         for (int j = nct; j < nu; j++) {
            for (int i = 0; i < m; i++) {
               U[i][j] = 0.0;
            }
            U[j][j] = 1.0;
         } 
          for (int k = nct-1; k >= 0; k--) {
            if (s[k] != 0.0) {
               for (int j = k+1; j < nu; j++) {
                  double t = 0;
                  for (int i = k; i < m; i++) {
                     t += U[i][k]*U[i][j];
                  }
                  t = -t/U[k][k];
                  for (int i = k; i < m; i++) {
                     U[i][j] += t*U[i][k];
                  }
               }
               for (int i = k; i < m; i++ ) {
                  U[i][k] = -U[i][k];
               }
               U[k][k] = 1.0 + U[k][k];
               for (int i = 0; i < k-1; i++) {
                  U[i][k] = 0.0;
               }
            } else {
               for (int i = 0; i < m; i++) {
                  U[i][k] = 0.0;
               }
               U[k][k] = 1.0;
            }
         }  
      }

      if (wantv) {
         for (int k = n-1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
               for (int j = k+1; j < nu; j++) {
                  double t = 0;
                  for (int i = k+1; i < n; i++) {
                     t += V[i][k]*V[i][j];
                  }
                  t = -t/V[k+1][k];
                  for (int i = k+1; i < n; i++) {
                     V[i][j] += t*V[i][k];
                  }
               }
            }
            for (int i = 0; i < n; i++) {
               V[i][k] = 0.0;
            }
            V[k][k] = 1.0;
         }
      }

      int pp = p-1;
      int iter = 0;
      double eps = Math.pow(2.0,-52.0);
      double tiny = Math.pow(2.0,-966.0);
      
      while (p > 0) {
         int k,kase;

        for (k = p-2; k >= -1; k--) {
            if (k == -1) {
               break;
            }
            if (Math.abs(e[k]) <=
                  tiny + eps*(Math.abs(s[k]) + Math.abs(s[k+1]))) {
               e[k] = 0.0;
               break;
            }
         }
        if (k == p-2) {
            kase = 4;
         } else {
            int ks;
            for (ks = p-1; ks >= k; ks--) {
               if (ks == k) {
                  break;
               }
               double t = (ks != p ? Math.abs(e[ks]) : 0.) + 
                          (ks != k+1 ? Math.abs(e[ks-1]) : 0.);
               if (Math.abs(s[ks]) <= tiny + eps*t)  {
                  s[ks] = 0.0;
                  break;
               }
            }
            if (ks == k) {
               kase = 3;
            } else if (ks == p-1) {
               kase = 1;
            } else {
               kase = 2;
               k = ks;
            }
         }
         k++;  

        switch (kase) {

            case 1: {
               double f = e[p-2];
               e[p-2] = 0.0;
               for (int j = p-2; j >= k; j--) {
                  double t = Maths.hypot(s[j],f);
                  double cs = s[j]/t;
                  double sn = f/t;
                  s[j] = t;
                  if (j != k) {
                     f = -sn*e[j-1];
                     e[j-1] = cs*e[j-1];
                  }
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][p-1];
                        V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
                        V[i][j] = t;
                     }
                  }
               }
            }
            break;

            case 2: {
               double f = e[k-1];
               e[k-1] = 0.0;
               for (int j = k; j < p; j++) {
                  double t = Maths.hypot(s[j],f);
                  double cs = s[j]/t;
                  double sn = f/t;
                  s[j] = t;
                  f = -sn*e[j];
                  e[j] = cs*e[j];
                  if (wantu) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][k-1];
                        U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
                        U[i][j] = t;
                     }
                  }
               }
            }
            break;

            case 3: {

               double scale = Math.max(Math.max(Math.max(Math.max(
                       Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])), 
                       Math.abs(s[k])),Math.abs(e[k]));
               double sp = s[p-1]/scale;
               double spm1 = s[p-2]/scale;
               double epm1 = e[p-2]/scale;
               double sk = s[k]/scale;
               double ek = e[k]/scale;
               double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
               double c = (sp*epm1)*(sp*epm1);
               double shift = 0.0;
               if ((b != 0.0) | (c != 0.0)) {
                  shift = Math.sqrt(b*b + c);
                  if (b < 0.0) {
                     shift = -shift;
                  }
                  shift = c/(b + shift);
               }
               double f = (sk + sp)*(sk - sp) + shift;
               double g = sk*ek;
   
               for (int j = k; j < p-1; j++) {
                  double t = Maths.hypot(f,g);
                  double cs = f/t;
                  double sn = g/t;
                  if (j != k) {
                     e[j-1] = t;
                  }
                  f = cs*s[j] + sn*e[j];
                  e[j] = cs*e[j] - sn*s[j];
                  g = sn*s[j+1];
                  s[j+1] = cs*s[j+1];
                  if (wantv) {
                     for (int i = 0; i < n; i++) {
                        t = cs*V[i][j] + sn*V[i][j+1];
                        V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
                        V[i][j] = t;
                     }
                  }
                  t = Maths.hypot(f,g);
                  cs = f/t;
                  sn = g/t;
                  s[j] = t;
                  f = cs*e[j] + sn*s[j+1];
                  s[j+1] = -sn*e[j] + cs*s[j+1];
                  g = sn*e[j+1];
                  e[j+1] = cs*e[j+1];
                  if (wantu && (j < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = cs*U[i][j] + sn*U[i][j+1];
                        U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
                        U[i][j] = t;
                     }
                  }
               }
               e[p-2] = f;
               iter = iter + 1;
            }
            break;

            case 4: {

               if (s[k] <= 0.0) {
                  s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                  if (wantv) {
                     for (int i = 0; i <= pp; i++) {
                        V[i][k] = -V[i][k];
                     }
                  }
               }
   
               while (k < pp) {
                  if (s[k] >= s[k+1]) {
                     break;
                  }
                  double t = s[k];
                  s[k] = s[k+1];
                  s[k+1] = t;
                  if (wantv && (k < n-1)) {
                     for (int i = 0; i < n; i++) {
                        t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
                     }
                  }
                  if (wantu && (k < m-1)) {
                     for (int i = 0; i < m; i++) {
                        t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
                     }
                  }
                  k++;
               }
               iter = 0;
               p--;
            }
            break;
         }
      }
        
        
    }
   
    public static Matrix random(int m, int n) {
        Matrix A = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A.data[i][j] = Math.random();
        return A;
    }

    
    public int getColumnCount(){
        return this.n;
    }
    
    public int getRowCount(){
        return this.m; 
    }
    
    public static Matrix identity(int n) {
        
        Matrix I = new Matrix(n, n);
        for (int i = 0; i < n; i++)
            I.data[i][i] = 1;
        return I;
    }

    private void swap(int i, int j) {
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }

    public Matrix transpose() {
        Matrix A = new Matrix(n,m);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }

    public double[][] getData(){
        return this.data;
    }
    
    public double getData(int m, int n){
        return this.data[m][n];
    }
    
    // return C = A + B
    public Matrix plus(Matrix B) {
        Matrix A = this;
        if (B.m != A.m || B.n != A.n) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }

    // return C = A + b
    public Matrix plus(Double b) {
        Matrix A = this;
        
        Matrix C = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C.data[i][j] = A.data[i][j] + b;
        return C;
    }
    
    // return C = A - B
    public Matrix minus(Matrix B) {
        Matrix A = this;
        if (B.m != A.m || B.n != A.n) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        return C;
    }
    
    // return C = A - b
    public Matrix minus(Double b) {
        Matrix A = this;
        
        Matrix C = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                C.data[i][j] = A.data[i][j] - b;
        return C;
    }

    // does A = B exactly?
    public boolean eq(Matrix B) {
        Matrix A = this;
        if (B.m != A.m || B.n != A.n) throw new RuntimeException("Illegal matrix dimensions.");
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                if (A.data[i][j] != B.data[i][j]) return false;
        return true;
    }

    // return C = A * B
    public Matrix times(Matrix B) {
        Matrix A = this;
        if (A.n != B.m) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.m, B.n);
        for (int i = 0; i < C.m; i++)
            for (int j = 0; j < C.n; j++)
                for (int k = 0; k < A.n; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }

    // return C = A .* b
    public Matrix times(double b){
        Matrix A = this;
        Matrix C = new Matrix(A.m, A.n);
        for (int i = 0; i < C.m; i++){
            for (int j = 0; j < C.n; j++){
                C.data[i][j] = A.data[i][j]*b;
            }
        }
          return C;
    }
    
    public double[] getDiagonal(){
        
        Matrix A = this;
        int c = m;
        
        if(n<m){
            c = n;
        }
        double resp[] = new double[c]; 
        
        for(int i = 0; i < c; i++){
            resp[i] = A.data[i][i];
        }
        return resp; 
    }
    
    public Matrix getColumns(int[] vec){
        try{
        double[][] tmp = new double[this.m][vec.length];
        Matrix A = this;
        
        for(int i = 0; i < vec.length; i++){
            if(vec[i] > this.n){ 
                throw new RuntimeException("Illegal vector extends matrix size.");
            }
        }
        for(int i = 0; i < this.m; i++){
            for(int j = 0; j < vec.length; j++){
                tmp[i][j] = A.data[i][vec[j]];
            }
        }
        return new Matrix(tmp);
        }catch(Exception e){
            System.out.println(e.getLocalizedMessage());
            return null; 
            
        }
    }

    // return x = A^-1 b, assuming A is square and has full rank
    public Matrix solve(Matrix rhs) {
        if (m != n || rhs.m != n || rhs.n != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this.data);
        Matrix b = new Matrix(rhs.data);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < n; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < n; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < n; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < n; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < n; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        // back substitution
        Matrix x = new Matrix(n, 1);
        for (int j = n - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < n; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        return x;
   
    }

    // print matrix to standard output
    public void show() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                System.out.printf("%9.4f ", data[i][j]);
            }
            System.out.println();
        }
    }
    
    public void showSize(){
        System.out.println("size " + this.m + " rows; " + this.n  + " columns;");
    }
    
    public Matrix invert(){
        
        double x[][] = new double[this.n][this.n];
        double b[][] = new double[this.n][this.n];
        int index[] = new int[this.n];
        
        for (int i=0; i<this.n; ++i){
            b[i][i] = 1;
        }
        
        // Transform the matrix into an upper triangle
        Matrix A = gaussian(index);
        double a[][] = A.getData();
        
        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<this.n-1; ++i)
            for (int j=i+1; j<this.n; ++j)
                for (int k=0; k<this.n; ++k)
                    b[index[j]][k]
                    	    -= a[index[j]][i]*b[index[i]][k];

        for (int i=0; i<this.n; ++i) 
        {
            x[this.n-1][i] = b[index[this.n-1]][i]/a[index[this.n-1]][this.n-1];
            for (int j=this.n-2; j>=0; --j) 
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<this.n; ++k) 
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
 
        return new Matrix(x); 
    }

    private Matrix gaussian(int[] index) {
   
       Matrix A = new Matrix(this.data);
       double a[][] = A.data; 
        
       double c[] = new double[this.n];
        
       // Initialize the index
        for (int i=0; i<this.n; ++i){
              index[i] = i;
        }
          
       // Find the rescaling factors, one from each row
        for (int i=0; i<this.n; ++i) 
        {
            double c1 = 0;
            for (int j=0; j<this.n; ++j) 
            {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }
        
        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<this.n-1; ++j) 
        {
            double pi1 = 0;
            for (int i=j; i<this.n; ++i) 
            {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }
        
        // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<this.n; ++i) 	
            {
                double pj = a[index[i]][j]/a[index[j]][j];
 
            // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;
 
            // Modify other elements accordingly
                for (int l=j+1; l<this.n; ++l){
                    a[index[i]][l] -= pj*a[index[j]][l];
                }
            }
        }
        
        return new Matrix(a); 
        
    }
    
     public double norm1 () {
      double f = 0;
      for (int j = 0; j < this.n; j++) {
         double s = 0;
         for (int i = 0; i < this.m; i++) {
            s += Math.abs(this.data[i][j]);
         }
         f = Math.max(f,s);
      }
      return f;
   }

     

    public double norm2 () {
        
      calcParameter();  
      return s[0];
   } 
    
     public double normInf () {
      double f = 0;
      for (int i = 0; i < m; i++) {
         double s = 0;
         for (int j = 0; j < n; j++) {
            s += Math.abs(this.data[i][j]);
         }
         f = Math.max(f,s);
      }
      return f;
   }
     

   public double normF () {
      double f = 0;
      for (int i = 0; i < this.m; i++) {
         for (int j = 0; j < this.n; j++) {
            f = Maths.hypot(f,this.data[i][j]);
         }
      }
      return f;
   }  

   
   public Matrix normalize(Double[] factor){
       
       Matrix nData = new Matrix(this.data); 
       
       for(int i = 0; i < nData.getColumnCount(); i++){
           
           for(int j = 0; j < nData.getRowCount(); j++){
              if(factor[i] != null && !factor[i].equals(0.0) && !factor[i].isNaN()){
                  try{
               nData.data[j][i] =  nData.data[j][i]/factor[i];
               }catch(Exception e){
                   logger.error(e.getMessage() + "for entry " + j + ":" + i);
               }
              }
           }
       }
       return nData;
   }
   
   public Matrix denormalize(Double[] factor){
       
       Matrix nData = new Matrix(this.data); 
       
       for(int i = 0; i < nData.getColumnCount(); i++){
           
           for(int j = 0; j < nData.getRowCount(); j++){
              if(factor[i] != null && !factor[i].equals(0.0) && !factor[i].isNaN()){
                  try{
               nData.data[j][i] =  nData.data[j][i]*factor[i];
               }catch(Exception e){
                   logger.error(e.getMessage() + "for entry " + j + ":" + i);
               }
              }
           }
       }
       return nData;
   }
   
   public Double[] getNormalizationFactor(){
       
       Double[] factor = new Double[this.getColumnCount()];
       
       for(int i = 0; i < this.getColumnCount(); i++){
           factor[i] = 0.0; 
           for(int j = 0; j < this.getRowCount(); j++){
               if(Math.abs(this.data[j][i]) > factor[i] ){
                   factor[i] = Math.abs(this.data[j][i]);
               }
           }
       }
       return factor; 
   }
   
   
   public int rank () {
       
      calcParameter(); 
      double eps = Math.pow(2.0,-52.0);
      double tol = Math.max(m,n)*s[0]*eps;
      int r = 0;
      for (int i = 0; i < s.length; i++) {
         if (s[i] > tol) {
            r++;
         }
      }
      return r;
   }

    public void display() {
        for (int row = 0; row < data.length; row++) {
            for (int column = 0; column < data[row].length; column++) {
                System.out.print(data[row][column] + " ");
            }
            System.out.println();
        }
    }


   
}
