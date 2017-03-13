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
import com.google.gson.FieldNamingPolicy;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import flexjson.JSONDeserializer;
import flexjson.JSONSerializer;
import java.awt.AWTException;
import java.awt.Rectangle;
import java.awt.Robot;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Reader;
import javax.imageio.ImageIO;
import javax.swing.JPanel;
import javax.swing.table.TableModel;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.jdesktop.swingx.JXTable;
import org.mantoQSAR.core.descriptor.ProjectDescriptor;
import org.mantoQSAR.core.util.ColorStatic;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;





public class ImportExport {

    Logger logger; 
    
    public ImportExport(){

        this.logger = LoggerFactory.getLogger(ImportExport.class);
    }

    public void saveTableToExcel(JXTable table, File file) throws FileNotFoundException, IOException {

    Workbook wb = new XSSFWorkbook(); //Excell workbook
    Sheet sheet = wb.createSheet(); //WorkSheet
    Row row = sheet.createRow(2); //Row created at line 3
    TableModel model = table.getModel(); //Table model


    Row headerRow = sheet.createRow(0); //Create row at line 0
    for(int headings = 0; headings < model.getColumnCount(); headings++){ //For each column
        headerRow.createCell(headings).setCellValue(model.getColumnName(headings));//Write column name
    }

    for(int rows = 0; rows < model.getRowCount(); rows++){ //For each table row
        for(int cols = 0; cols < table.getColumnCount(); cols++){ //For each table column
            row.createCell(cols).setCellValue(model.getValueAt(rows, cols).toString()); //Write value
        }

        row = sheet.createRow((rows + 3)); 
    }
    wb.write(new FileOutputStream(file.toString()));//Save the file     
}
    

    public void saveTableToCSV(JXTable table, File file)  throws FileNotFoundException, IOException{

        FileWriter fw;

            TableModel model = table.getModel();
            fw = new FileWriter(file);
            for (int i = 0; i < model.getColumnCount(); i++) {
                fw.write(model.getColumnName(i) + "\t");
            }   fw.write("\n");
            for (int i = 0; i < model.getRowCount(); i++) {
                for (int j = 0; j < model.getColumnCount(); j++) {
                    fw.write(model.getValueAt(i, j).toString() + "\t");
                }
                fw.write("\n");
            }   fw.close();
       
            try {
                fw.close();
            } catch (IOException ex) {
                System.out.println(ColorStatic.RED + ex.getMessage() + ColorStatic.RESET);
            }
    }

    public void savePanelToImage(JPanel panel, File file){
        
        try {
            BufferedImage image = new Robot().createScreenCapture(new Rectangle(panel.getLocationOnScreen().x, 
                    panel.getLocationOnScreen().y, panel.getWidth(), panel.getHeight()));

                ImageIO.write(image, "png", file);
           
        } catch (IOException | AWTException ex) {
           logger.error(ex.getMessage());
        }
        
        
    }
    
    public void exportToJson(Object object, File f){
        
       Gson gson = new GsonBuilder()
             .disableHtmlEscaping()
             .setFieldNamingPolicy(FieldNamingPolicy.UPPER_CAMEL_CASE)
             .serializeSpecialFloatingPointValues()
             .setPrettyPrinting()
             .serializeNulls()
             .create();
        
       
        try {
            String sFile = f.getCanonicalPath();
            f.delete();
            
            f = new File(sFile);
        } catch (IOException ex) {
            logger.error(ex.getMessage());
        }

	String json = gson.toJson(object);
        
	try {
           try (
                   FileWriter writer = new FileWriter(f)) {
               writer.write(json);

           }

	} catch (IOException e) {
            logger.error(e.getMessage());
	}   
    }

    public void exportDescriptorListToJson(ProjectDescriptor projDescr, File f){
                
       
        try {
            String sFile = f.getCanonicalPath();
            f.delete();
            
            f = new File(sFile);
        } catch (IOException ex) {
            logger.error(ex.getMessage());
        }
       
	JSONSerializer js = new JSONSerializer(); 
        
        String json = js.exclude("molecule", "structure", "description")
                       .prettyPrint(true)
                       .deepSerialize(projDescr); 
	
        
        try {
           try (
                   FileWriter writer = new FileWriter(f)) {
               writer.write(json);
  
           }
	} catch (IOException e) {
            logger.error(e.getMessage());
	}   
    }
    

    
    public ProjectDescriptor importDescriptorList(File f){ 
     
        ProjectDescriptor descriptorList = new ProjectDescriptor();
      Reader reader;
             
                try {

                    reader = new BufferedReader(new FileReader(f));

                 descriptorList = new JSONDeserializer<ProjectDescriptor>().deserialize(reader, ProjectDescriptor.class);

                } catch (FileNotFoundException ex) {
                    logger.error(ex.getMessage());
                    System.out.println(ColorStatic.RED + "Import of descriptors failed." + ColorStatic.RESET);
                }         
                return descriptorList; 
    }

    
    public Structure importStructure(File file) {

        Structure s = null;
        
        System.out.println(ColorStatic.GREEN + "Import structure from " + file.toString() + ColorStatic.RESET);

        if (file.exists() && !file.isDirectory()) {

            try {
                PDBFileReader pdbreader = new PDBFileReader();

                s = pdbreader.getStructure(file.getAbsolutePath());

            } catch (IOException e) {

                logger.error("General I/O exception: " + e.getMessage());
            }

        } else {
            System.out.println(ColorStatic.RED + "No structure file found. Please check path and file name. " + ColorStatic.RESET);

        }

        return s;
    }
}
