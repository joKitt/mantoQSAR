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
package org.mantoQSAR.gui.model;

import java.awt.Color;
import java.awt.GridLayout;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import org.mantoQSAR.gui.resources.CollapsiblePane;

public class ModelStatisticPane extends JPanel {

    CollapsiblePane cPane;
    private String rSquare = "0.0";
    private String predRSquare = "0.0";

    public ModelStatisticPane() {

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        cPane = new CollapsiblePane();
        cPane.setCollapsed(false);
        cPane.setTitle("model statistics");
        setPanelContent();

        this.add(cPane);
    }

    private void setPanelContent() {

        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {

                JPanel mainPanel = cPane.getMainPanel();

                mainPanel.removeAll();
                mainPanel.setLayout(new GridLayout(0, 2));

                // rsquare text
                JLabel rsquareDescr = new JLabel("rsquare value ");
                rsquareDescr.setForeground(Color.DARK_GRAY);

                // rsquare value
                JLabel rsquareVal = new JLabel(rSquare);
                rsquareVal.setForeground(Color.DARK_GRAY);

                // pred rsquare text
                JLabel predrsquareDescr = new JLabel("pred. rsquare value ");
                predrsquareDescr.setForeground(Color.DARK_GRAY);

                // pred rsquare value
                JLabel predrsquareVal = new JLabel(predRSquare);
                predrsquareVal.setForeground(Color.DARK_GRAY);

                mainPanel.add(rsquareDescr);
                mainPanel.add(rsquareVal);
                mainPanel.add(predrsquareDescr);
                mainPanel.add(predrsquareVal);

                cPane.revalidate();
                cPane.repaint();
            }
        });
    }

    public void setRSquare(Double rSquare) {

        if (rSquare != null) {
            this.rSquare = Double.toString(rSquare);
        } else {
            this.rSquare = "null";
        }
        System.out.println("r square " + rSquare);
        setPanelContent();
    }

    public void setPredRSquare(Double predRSquare) {

        if (rSquare != null) {
            this.predRSquare = Double.toString(predRSquare);
        } else {
            this.predRSquare = "null";
        }

        System.out.println("pred r square " + predRSquare);
        setPanelContent();
    }
}
