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
import java.awt.FlowLayout;
import java.awt.Point;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.ObservationSet;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.descriptor.MoleculeStatic;
import org.mantoQSAR.gui.Statics;

public class MoleculeContentGroupPanel extends JPanel implements PropertyChangeListener {

    private Screen screen;
    private MoleculeController moleculeController;
    private JPanel containerPanel;
    private final JPanel contentPane;
    private final JScrollPane scrollPane;

    public MoleculeContentGroupPanel() {

        moleculeController = MoleculeController.getInstance();
        screen = moleculeController.getScreen();
        this.setLayout(new BorderLayout());
        contentPane = new JPanel();

        HeaderPanel headerPane = new HeaderPanel();
        JPanel centerPanel = new JPanel();
        centerPanel.setLayout(new BorderLayout());
        centerPanel.add(headerPane, BorderLayout.PAGE_START);
        centerPanel.add(contentPane, BorderLayout.CENTER);

        this.scrollPane = new JScrollPane(centerPanel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        this.add(scrollPane, BorderLayout.CENTER);

        initContentPanel();
        initListeners();
    }

    protected final void initListeners() {
        moleculeController.addPropertyChangeListener(this);
        screen.addPropertyChangeListener(this);

    }

    private void initContentPanel() {

        containerPanel = new JPanel();
        containerPanel.setLayout(new BoxLayout(containerPanel, BoxLayout.Y_AXIS));
        contentPane.setLayout(new BorderLayout());
        contentPane.setBackground(Color.WHITE);
        contentPane.add(containerPanel, BorderLayout.NORTH);

    }

    public void drawObservationList() {

        Point point = scrollPane.getViewport().getViewPosition();

        containerPanel.removeAll();

        for (int i = 0; i < screen.getMoleculeListSize(); i++) {
            Molecule m = screen.getMolecule(i);
            ObservationSet os = screen.getObservationSetList().get(i);

            // set an id for position in list 
            MoleculeContentPanel mContentPanel = new MoleculeContentPanel(os, m, i);

            if (i == moleculeController.getActivePosition()) {
                mContentPanel.setActive(true);
            } else {
                mContentPanel.setActive(false);
            }
            containerPanel.add(mContentPanel);
        }

        this.revalidate();
        this.repaint();
        this.scrollPane.getViewport().setViewPosition(point);
    }

    @Override
    public void propertyChange(PropertyChangeEvent evt) {

        if (evt.getPropertyName().equalsIgnoreCase(MoleculeStatic.CHANGE_LOAD_PROPERTY_SETTING)) {
            drawObservationList();
        }

        if (evt.getPropertyName().equalsIgnoreCase(MoleculeGuiStatic.CHANGE_ACTIVE_OBSERVATION)) {
            drawObservationList();
        }

    }

    private class HeaderPanel extends JPanel {

        private final JButton runCalcBtn;
        private final MoleculeController moleculeControl;

        public HeaderPanel() {

            moleculeControl = MoleculeController.getInstance();

            runCalcBtn = new JButton("Run Calculation");
            runCalcBtn.addActionListener(moleculeControl);
            runCalcBtn.setActionCommand(Statics.ACTION_CALC_DESCRIPTOR);
            runCalcBtn.setContentAreaFilled(false);

            JButton showSurfaceBtn = new JButton("Show surface");
            showSurfaceBtn.addActionListener(moleculeControl);
            showSurfaceBtn.setActionCommand(Statics.ACTION_SHOW_SURFACE);
            showSurfaceBtn.setContentAreaFilled(false);
            showSurfaceBtn.setForeground(Color.LIGHT_GRAY);

            setLayout(new FlowLayout(FlowLayout.RIGHT));
            setBackground(Color.WHITE);

            this.add(runCalcBtn);
            this.add(showSurfaceBtn);
        }

    }

}
