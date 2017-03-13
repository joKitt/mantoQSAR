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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import static java.awt.GridBagConstraints.HORIZONTAL;
import static java.awt.GridBagConstraints.REMAINDER;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import org.jdesktop.swingx.JXCollapsiblePane;

public class CollapsiblePane extends JPanel {

    public TitledBorder border;
    public JXCollapsiblePane cPane;
    private JPanel topPane;
    private JPanel mainPane;
    private JLabel labelToggle;

    public String title = " ";

    public CollapsiblePane() {

        cPane = new JXCollapsiblePane();
        cPane.setBackground(Color.WHITE);

        cPane.setLayout(new BorderLayout());

        topPane = new JPanel(new GridBagLayout());
        topPane.setBackground(Color.WHITE);

        mainPane = new JPanel(new GridBagLayout());
        mainPane.setBackground(Color.WHITE);

        border = new TitledBorder(BorderFactory.createLineBorder(Color.LIGHT_GRAY), title);

        this.setBorder(border);
        this.setLayout(new BorderLayout());
        this.setBackground(Color.WHITE);

        labelToggle = new JLabel("show");
        labelToggle.setForeground(Color.BLUE);

        labelToggle.addMouseListener(new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent e) {
                if (cPane.isCollapsed()) {
                    cPane.setCollapsed(true);
                    labelToggle.setText("hide");
                } else {
                    cPane.setCollapsed(false);
                    labelToggle.setText("show");
                }
            }
        });

        topPane.add(labelToggle, newHideShowConstraints());
        cPane.setCollapsed(true);
        cPane.add(mainPane, BorderLayout.PAGE_START);

        this.add(topPane, BorderLayout.PAGE_START);
        this.add(cPane, BorderLayout.CENTER);

    }

    public void setTitle(String title) {
        border.setTitle(title);

    }

    public void setCollapsed(boolean b) {
        cPane.setCollapsed(b);

    }

    private GridBagConstraints newConstraints() {
        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(2, 2, 2, 2);
        return c;
    }

    private GridBagConstraints newBtnConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.NORTHEAST;

        c.weightx = 1.0;
        c.fill = HORIZONTAL;
        c.gridwidth = REMAINDER;
        return c;
    }

    private GridBagConstraints newHideShowConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.FIRST_LINE_END;
        c.weightx = 0.0;
        c.gridwidth = GridBagConstraints.REMAINDER;
        return c;
    }

    public JPanel getTopPane() {
        return topPane;
    }

    public void setTopPane(JPanel topPane) {
        this.topPane = topPane;
    }

    public JPanel getMainPanel() {
        return mainPane;
    }

    public void setMainPanel(JPanel mainPane) {
        this.mainPane = mainPane;
    }

}
