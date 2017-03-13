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

package org.mantoQSAR.gui.molecule;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.List;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import org.jdesktop.swingx.JXCollapsiblePane;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.descriptor.Descriptor;
import org.mantoQSAR.core.descriptor.OrientationDescriptorGroup;
import org.mantoQSAR.core.descriptor.PlaneDescriptorGroup;
import org.mantoQSAR.core.util.ColorStatic;
import org.openscience.jmol.app.jmolpanel.console.AppConsole;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class MoleculeViewController extends JPanel implements PropertyChangeListener{
    
    private MoleculeController moleculeControl;
    private Screen screen; 
    private Logger logger; 
    
    JPanel controlPanel; 
    JPanel controlDetailPanel; 
    JPanel controlTopPanel; 
    JPanel comboBoxPanel;
    JComboBox<String> comboBox;
    JmolPanel viewPanel; 
    
    public MoleculeViewController(){

        logger = LoggerFactory.getLogger(MoleculeViewController.class);
        initMoleculeControl();
        initControlPanel();
        this.add(controlPanel, BorderLayout.CENTER);
    }
    
    private void initMoleculeControl() {
        // get screen
        moleculeControl = MoleculeController.getInstance();
        moleculeControl.addPropertyChangeListener(this);
        
        screen = moleculeControl.getScreen();
        screen.addPropertyChangeListener(this);
        viewPanel = moleculeControl.getViewPanel();
    }
    
  
    private void initControlPanel() {

        controlPanel = new JPanel();
        controlPanel.setBackground(Color.WHITE);
        controlPanel.setLayout(new BoxLayout(controlPanel, BoxLayout.PAGE_AXIS));

        JPanel groupPanelView = new JPanel();
        groupPanelView.setBackground(Color.WHITE);
        groupPanelView.setLayout(new BoxLayout(groupPanelView, BoxLayout.PAGE_AXIS));
        groupPanelView.setBorder(new TitledBorder("Edit View"));

        final JXCollapsiblePane cPanel = new JXCollapsiblePane();
        cPanel.setBackground(Color.WHITE);
        cPanel.setSize(new Dimension(0, 400));

        this.controlTopPanel = new JPanel();
        this.controlTopPanel.setBackground(Color.WHITE);

        this.controlTopPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));

        this.comboBoxPanel = new JPanel();
        this.comboBoxPanel.setBackground(Color.WHITE);
        this.controlTopPanel.add(comboBoxPanel);

        this.controlDetailPanel = new JPanel();
        this.controlDetailPanel.setBackground(Color.WHITE);

        final JLabel labelToggle = new JLabel("show");
        cPanel.setCollapsed(true);
        labelToggle.setForeground(Color.BLUE);

        labelToggle.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (cPanel.isCollapsed()) {
                    cPanel.setCollapsed(false);
                    labelToggle.setText("hide");
                } else {

                    cPanel.setCollapsed(true);
                    labelToggle.setText("show");
                }
            }
        });

        controlTopPanel.add(labelToggle);

        JPanel panelCon = new JPanel();
        panelCon.setBackground(Color.WHITE);

        panelCon.setLayout(new BorderLayout());
        panelCon.setPreferredSize(new Dimension(500, 300));
        AppConsole console = new AppConsole(this.viewPanel.getViewer(), panelCon, "History State Clear");
        
        
        panelCon.setBackground(Color.WHITE);
        this.viewPanel.getViewer().setJmolCallbackListener(console);
        
        
        try {
        Component[] cL = console.jcd.getComponents();
        for (Component cL1 : cL) {           
                cL1.setBackground(Color.WHITE);   
                }
        }catch(Exception e){
                logger.error(e.getMessage());
            }
        
        cPanel.add(panelCon, BorderLayout.NORTH);
        cPanel.setCollapsed(true);

        groupPanelView.add(controlTopPanel);
        groupPanelView.add(cPanel);

        controlPanel.add(groupPanelView);

    }

    // adds a drop down list in edit view panel for selection of orientation to be displayed
    private void drawOrientationSelector() {
       
       
        comboBoxPanel.removeAll();
        if(screen.getDescriptorGroup(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor())
           instanceof PlaneDescriptorGroup){
          
        OrientationDescriptorGroup oDG = (OrientationDescriptorGroup) screen.getDescriptorGroup(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor()); 
        
        comboBox = new JComboBox<>();
        
        List<Double[]> vector = oDG.getVector(); 
        int nOrientation = vector.size();
        System.out.println(ColorStatic.GREEN + "orientations to display " + nOrientation + ColorStatic.RESET);

        int n = this.screen.getDescriptorPositionInGroup(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor());
        
        for(int i = 0; i < nOrientation; i++){
            Descriptor d = oDG.getDescriptor().get(n);
            ArrayList<ArrayList<Descriptor>> descriptorDetail = oDG.getDescriptorDetail();
            
            Double y = vector.get(i)[1]*90;
            Double x = 0.0; 
            
            if(vector.get(i)[2] > 0){
                if(vector.get(i)[0] > 0){
                    x = (vector.get(i)[0])* 90+90; 
                }else{
                    x = (vector.get(i)[0]* -90)-90; 
                }
                    }else{
                if(vector.get(i)[0] > 0){
                    x = vector.get(i)[0]* 90; 
                }else{
                    x = vector.get(i)[0]* -90;
                }
                
            }
            
            String s = (i+1) + ": "  + String.format("%.2g%n",  x) + " " 
                                     + String.format("%.2g%n",  y) + "\t "
                                     + String.format("%.4g%n",  descriptorDetail.get(i).get(n).getValue());

            comboBox.addItem(s);  
            System.out.println(ColorStatic.PURPLE + s + ColorStatic.RESET);
            
        }
        comboBox.setBackground(Color.WHITE);
        comboBox.setEditable(false);
        
        Integer select = this.moleculeControl.getOrientationID(); 
        if(select != null){
          comboBox.setSelectedIndex(select);  
        }
        
        comboBox.addActionListener(new ActionListener(){
            
            @Override
            public void actionPerformed(ActionEvent e){
                @SuppressWarnings("unchecked")
                JComboBox<String> combo = (JComboBox)e.getSource();
                int currentIndex = combo.getSelectedIndex();
                moleculeControl.setOrientationID(currentIndex);
                System.out.println(ColorStatic.PURPLE + "set active orientation index " + currentIndex + ColorStatic.RESET);
            }
        });
        
        comboBoxPanel.add(comboBox);   
        }else{

            comboBoxPanel.removeAll();
        }

        this.repaint();
    }
    
    public JPanel getViewControlPanel() {
        return controlPanel;
    }
    
    @Override
    public void propertyChange(PropertyChangeEvent evt) {
        
        if (evt.getPropertyName().equalsIgnoreCase(MoleculeGuiStatic.CHANGE_ACTIVE_DESCRIPTOR)) {
            drawOrientationSelector();
        }
    }
}
