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
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.List;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.table.AbstractTableModel;
import org.jdesktop.swingx.JXTable;
import org.jdesktop.swingx.decorator.Highlighter;
import org.jdesktop.swingx.decorator.HighlighterFactory;
import org.mantoQSAR.core.ObservationSet;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.descriptor.MoleculeStatic;
import org.mantoQSAR.core.math.Matrix;
import org.mantoQSAR.core.util.ColorStatic;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DescriptorTablePanel extends JPanel implements PropertyChangeListener {

    private Matrix descriptorData;
    Logger logger;
    DescriptorTableModel descriptorTableModel;
    Screen screen;
    MoleculeController moleculeControl;
    List<ObservationSet> observationSetList;

    public DescriptorTablePanel() {

        this.logger = LoggerFactory.getLogger(DescriptorTablePanel.class);
        this.moleculeControl = MoleculeController.getInstance();
        this.screen = moleculeControl.getScreen();
        this.screen.addPropertyChangeListener(this);
        this.observationSetList = screen.getObservationSetList();
        this.descriptorData = screen.getDescriptorMatrix();

        descriptorTableModel = new DescriptorTableModel();
        moleculeControl.descriptorTable = createTable();

    }

    private JXTable createTable() {

        final JXTable table = new JXTable(descriptorTableModel);

        table.setAutoResizeMode(JXTable.AUTO_RESIZE_OFF);
        table.setAutoCreateColumnsFromModel(true);
        table.setPreferredScrollableViewportSize(table.getPreferredSize());

        configureHighlighters(table);

        table.addMouseListener(new MouseAdapter() {

            @Override
            public void mouseClicked(MouseEvent evt) {

                if (evt.getClickCount() == 2 && !evt.isConsumed()) {
                    evt.consume();

                    int row = table.rowAtPoint(evt.getPoint());
                    int column = table.columnAtPoint(evt.getPoint());

                    moleculeControl.setActivePosition(row);

                    if (column > 3) {
                        moleculeControl.setActiveDescriptor(column - 4);
                    }
                }
            }
        });

        JScrollPane scrollPane = new JScrollPane(table);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

        this.setLayout(new BorderLayout());
        this.add(scrollPane, BorderLayout.CENTER);

        return table;
    }

    private static void configureHighlighters(JXTable table) {

        Highlighter highlighter = HighlighterFactory.createAlternateStriping(Color.WHITE, new Color(235, 235, 255));
        table.addHighlighter(highlighter);
    }

    private JXTable setColumnWidth(JXTable table) {

        for (int i = 0; i < table.getColumnModel().getColumnCount(); i++) {
            table.getColumnModel().getColumn(i).setMinWidth(100);
            table.getColumnModel().getColumn(i).setPreferredWidth(100);
        }
        return table;
    }

    @Override
    public void propertyChange(PropertyChangeEvent evt) {

        logger.debug("PropertyChange call " + evt.getPropertyName());

        if (evt.getPropertyName().equalsIgnoreCase(MoleculeStatic.CHANGE_CALC_DESCRIPTOR)) {

            this.descriptorData = screen.getDescriptorMatrix();
            this.observationSetList = screen.getObservationSetList();
            this.descriptorTableModel.fireTableStructureChanged();

            moleculeControl.descriptorTable = this.setColumnWidth(moleculeControl.descriptorTable);
            moleculeControl.descriptorTable.revalidate();
        }

    }

    private class DescriptorTableModel extends AbstractTableModel {

        @Override
        public int getRowCount() {
            return descriptorData.getRowCount();
        }

        @Override
        public int getColumnCount() {
            return descriptorData.getColumnCount() + 4;
        }

        @Override
        public String getColumnName(int column) {

            String result;

            switch (column) {
                case 0:
                    result = "index";
                    break;
                case 1:
                    result = "molecule name";
                    break;
                case 2:
                    result = "pH value";
                    break;
                case 3:
                    result = "ionic strength";
                    break;
                default:
                    List<String> descrName = screen.getDescriptorName();
                    if (descrName != null) {
                        result = descrName.get(column - 4);
                    } else {
                        result = Integer.toString(column);
                    }
                    break;
            }

            return result;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {

            Object result = null;

            switch (columnIndex) {

                case 0:
                    result = rowIndex;
                    break;

                case 1:

                    if (observationSetList != null && !observationSetList.isEmpty() && observationSetList.size() > rowIndex) {
                        result = observationSetList.get(rowIndex).getName();
                    } else {
                        result = "empty";
                    }

                    break;

                case 2:
                    if (observationSetList != null && !observationSetList.isEmpty() && observationSetList.size() > rowIndex) {
                        result = observationSetList.get(rowIndex).getCondition().pH.toString();
                    } else {
                        result = "empty";
                    }

                    break;

                case 3:
                    if (observationSetList != null && !observationSetList.isEmpty() && observationSetList.size() > rowIndex) {
                        result = observationSetList.get(rowIndex).getCondition().getIonicStrength().toString();
                    } else {
                        result = "empty";
                    }

                    break;

                default:
                    try {
                        result = descriptorData.data[rowIndex][columnIndex - 4];
                    } catch (Exception e) {
                        System.out.println(ColorStatic.PURPLE + e.getMessage() + ColorStatic.RESET);
                    }
                    ;

            }
            return result;
        }

    }

}
