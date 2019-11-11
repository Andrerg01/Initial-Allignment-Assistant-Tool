#Importing sys so we can importy Optics.py
import sys
sys.path.insert(0, '~/Gaussian Beam Propagation')
#Importing all the definitions form Optics.py
from Optics import *

#Imports numpy
import numpy as np
#Imports improtant functions for juptyer display
from IPython.core.display import display, HTML, Markdown, clear_output
#Improts matplotlib
import matplotlib.pyplot as plt
#Improts widgets for buttons and text boxes for jupyter display
import ipywidgets as widgets
#Imports parser for ini files
import configparser
#Import os to see files and directories
import os

#Determines current working directory
#global workingDir
workingDir = os.getcwd()

#List of all the systems in the 'Systems' folder
#global systemsFiles
systemsFiles = [f for f in os.listdir(workingDir + '/Systems') if os.path.isfile(os.path.join(workingDir + '/Systems', f))]

#Removes the ".ini" from the end of the files
#global systems
systems = [system[:-4] for system in systemsFiles]
#Alphabetizes the systems
systems.sort()

#global systemIndex
systemIndex = 0

def beamFromFile(filepath):
    config = configparser.ConfigParser()
    config.read(filepath)
    return Beam(radiusOfCurvature = eval(config['Beam']['radiusOfCurvature']), width = eval(config['Beam']['width']), direction = eval(config['Beam']['direction']), position = eval(config['Beam']['position']), wavelength = eval(config['Beam']['wavelength']), indexOfRefraction = eval(config['Beam']['indexOfRefraction']), verbose = eval(config['Beam']['verbose']))

def elementsFromFile(filepath):
    config = configparser.ConfigParser()
    config.read(filepath)
    count = 1
    elements = []
    elementNow = 'Element'+'{:02d}'.format(count)
    while elementNow in config:
        if config[elementNow]['Type'] == 'Mirror':
            elements.append(Mirror(ID = config[elementNow]['ID'], radiusOfCurvature = eval(config[elementNow]['radiusOfCurvature']), positionOfCM = eval(config[elementNow]['positionOfCM']), parameter_d = eval(config[elementNow]['parameter_d']), yaw = eval(config[elementNow]['yaw']), pitch = eval(config[elementNow]['pitch']), diameter = eval(config[elementNow]['diameter'])))
        elif config[elementNow]['Type'] == 'Lens':
            elements.append(Lens(ID = config[elementNow]['ID'], radiusOfCurvature = eval(config[elementNow]['radiusOfCurvature']), positionOfCM = eval(config[elementNow]['positionOfCM']), parameter_d = eval(config[elementNow]['parameter_d']), yaw = eval(config[elementNow]['yaw']), pitch = eval(config[elementNow]['pitch']), diameter = eval(config[elementNow]['diameter']), indexOfRefraction = eval(config[elementNow]['indexOfRefraction'])))
        
        count = count + 1
        elementNow = 'Element'+'{:02d}'.format(count)
    return elements

def calculatePlotParameters(beam0, elements0, elementView, elementControl, extraYaws, extraPitches):
    #Makes a copy of the beam
    beam = beam0.copy()
    beamDefault = beam0.copy()
    #Makes a copy of all the elements
    elements = [element.copy() for element in elements0]
    elementsDefault = [element.copy() for element in elements0]

    #Shifts the pitch and yaw of all the elements accordingly
    for i in range(len(elements)):
        elements[i].pitch = elements[i].pitch + extraPitches[i]
        elements[i].yaw = elements[i].yaw + extraYaws[i]
    #Finds the ID os the element that is being controlled (not useful yet, but no harm in it).
    for iCtrl in range(len(elements)):
        if elements[iCtrl].ID == elementControl:
            break
    #Finds the ID os the element that is being viewed.
    for iView in range(len(elements)):
        if elements[iView].ID == elementView:
            break
            
    #Calculates the states of the beam at each element.
    states = beam.calculateStates(elements)
    statesDefault = beam.calculateStates(elementsDefault)
    #Calculates the flags for non-intersection and clipping
    flags = beam.calculateFlags(elements)
    
    #Returns the dictionary with all necessary information for the plot.
    return {
            "Title": "View: " + str(elementView) + " | Control: " + str(elementControl),
            "Range": np.array([-1, 1])*elements[iView].diameter*1.1/2.0,
            "BeamCenter": rotatePitchYaw(states[elementView].position - elements[iView].positionOfCM, -elements[iView].pitch, -elements[iView].yaw)[1:],
            "BeamRadius": states[elementView].width,
            "BeamCenterDefault": rotatePitchYaw(statesDefault[elementView].position - elementsDefault[iView].positionOfCM, -elementsDefault[iView].pitch, -elementsDefault[iView].yaw)[1:],
            "BeamRadiusDefault": statesDefault[elementView].width,
            "Flags": flags            
           }

def resetWidgets():
    systemsFiles = [f for f in os.listdir(workingDir + '/Systems') if os.path.isfile(os.path.join(workingDir + '/Systems', f))]

    #Removes the ".ini" from the end of the files
    #global systems
    systems = [system[:-4] for system in systemsFiles]
    #Alphabetizes the systems
    systems.sort()

    #global systemIndex
    systemIndex = 0
    #global beam
    beam = beamFromFile('Systems/' + systems[systemIndex] + '.ini')
    #global elements
    elements = elementsFromFile('Systems/' + systems[systemIndex] + '.ini')

    #Makes list of the ID of the optical elements.
    #global elementsIDs
    elementsIDs = np.array([element.ID for element in elements])
    #Default element to be viewed (last in the series).
    #global elementView
    elementView = elementsIDs[-1]
    #Default element to be controlled (last in the series).
    #global elementControl
    elementControl = elementsIDs[0]

    #Finds index of the element being viewed.
    #global viewIndex
    viewIndex = np.argwhere(elementsIDs == elementView)[0][0]
    #Finds index of the element being controlled.
    #global controlIndex
    controlIndex = np.argwhere(elementsIDs == elementControl)[0][0]

    #Initialized the "extra Yaws" to be zero for all elements.
    #global extraYaws
    extraYaws = np.array([0.0 for element in elements])
    #Initialized the "extra Pitches" to be zero for all elements.
    #global extraPitches
    extraPitches = np.array([0.0 for element in elements])

    #Defines the System Selection Dropdown menu.
    #global systemSelectionDropdown
    systemSelectionDropdown = widgets.Dropdown(value = systems[systemIndex], options = systems, description = "System Selection", width = '100%', style = {'description_width': '35%'})

    #Dfines Yaw and Pitch text boxes.
    #global yawTextBox
    yawTextBox = widgets.Text(value = str(extraYaws[controlIndex]), description='Extra Yaw', layout = widgets.Layout(width = '15%'))
    #global pitchTextBox
    pitchTextBox = widgets.Text(value = str(extraPitches[controlIndex]), description='Extra Pitch', layout = widgets.Layout(width = '15%'))
    
    #Defines element view and controll selection dropdown menus.
    #global elementViewDropdown
    elementViewDropdown = widgets.Dropdown(options = elementsIDs, value = elementView, layout = widgets.Layout(width = '20%'), description = 'View Element: ', style = {'description_width': '50%'})
    #global elementControlDropdown
    elementControlDropdown = widgets.Dropdown(options = elementsIDs, value = elementControl, layout = widgets.Layout(width = '20%'), description = 'Control Element: ', style = {'description_width': '50%'})

    #Defines "Show Default Beam" checkbox.
    #global showDefaultBeamCheckbox
    showDefaultBeamCheckbox.value = False

    #Defines "Show displacement" checkbox.
    #global showDisplacementCheckbox
    showDisplacementCheckbox.value = False

    #Defines Flag Font Size text box.
    #global flagFontSizeTextBox
    flagFontSizeTextBox.value = '12'

    #Defines Plot Font Size text box.
    #global plotFontSizeTextBox
    plotFontSizeTextBox.value = '12'

    #Defines Axes Ticks text box.
    #global axesTicksTextBox
    axesTicksTextBox.value = '7'

    #Defines Nummber of Circles text box
    #global numberOfCirclesTextBox
    numberOfCirclesTextBox.value = '10'

    #Defines Nummber of Lines text box
    #global numberOfLinesTextBox
    numberOfLinesTextBox.value = '16'

    #Defines Beam Color dropdown menu
    #global beamColorDropdown
    beamColorDropdown.value = 'Red'

#global beam
beam = beamFromFile('Systems/' + systems[systemIndex] + '.ini')
#global elements
elements = elementsFromFile('Systems/' + systems[systemIndex] + '.ini')

#Makes list of the ID of the optical elements.
#global elementsIDs
elementsIDs = np.array([element.ID for element in elements])
#Default element to be viewed (last in the series).
#global elementView
elementView = elementsIDs[-1]
#Default element to be controlled (last in the series).
#global elementControl
elementControl = elementsIDs[0]

#Finds index of the element being viewed.
#global viewIndex
viewIndex = np.argwhere(elementsIDs == elementView)[0][0]
#Finds index of the element being controlled.
#global controlIndex
controlIndex = np.argwhere(elementsIDs == elementControl)[0][0]

#Initialized the "extra Yaws" to be zero for all elements.
#global extraYaws
extraYaws = np.array([0.0 for element in elements])
#Initialized the "extra Pitches" to be zero for all elements.
#global extraPitches
extraPitches = np.array([0.0 for element in elements])
                   
#Defines the System Selection Dropdown menu.
#global systemSelectionDropdown
systemSelectionDropdown = widgets.Dropdown(value = systems[systemIndex], options = systems, description = "System Selection", width = '100%', style = {'description_width': '35%'})
   
#Defines Set System Button
#global selectSystemButton
selectSystemButton = widgets.Button(description = 'Select System!', layout = widgets.Layout(width = '99%'), button_style = 'primary')

#Defines the "Show Plot!" button.
#global showPlotButton
showPlotButton = widgets.Button(description = 'Show Plot!', layout = widgets.Layout(width = '99%'), button_style = 'success')

#Dfines Yaw and Pitch text boxes.
#global yawTextBox
yawTextBox = widgets.Text(value = str(extraYaws[controlIndex]), description='Extra Yaw', layout = widgets.Layout(width = '15%'))
#global pitchTextBox
pitchTextBox = widgets.Text(value = str(extraPitches[controlIndex]), description='Extra Pitch', layout = widgets.Layout(width = '15%'))
#Defines Yaw and Pitch unit selection dropdown menus.
#global yawUnitsDropdown
yawUnitsDropdown = widgets.Dropdown(options = ['Radians', 'Milli Radians', 'Micro Radians', 'Nano Radians'], layout = widgets.Layout(width = '10%'))
#global pitchUnitsDropdown
pitchUnitsDropdown = widgets.Dropdown(options = ['Radians', 'Milli Radians', 'Micro Radians', 'Nano Radians'], layout = widgets.Layout(width = '10%'))

#Defines element view and controll selection dropdown menus.
#global elementViewDropdown
elementViewDropdown = widgets.Dropdown(options = elementsIDs, value = elementView, layout = widgets.Layout(width = '20%'), description = 'View Element: ', style = {'description_width': '50%'})
#global elementControlDropdown
elementControlDropdown = widgets.Dropdown(options = elementsIDs, value = elementControl, layout = widgets.Layout(width = '20%'), description = 'Control Element: ', style = {'description_width': '50%'})

#Definds "Set Element" button.
#global setElementButton
setElementButton = widgets.Button(description = 'Set Elements!', layout = widgets.Layout(width = '99%'), button_style = 'primary')

#Defines "Show Default Beam" checkbox.
#global showDefaultBeamCheckbox
showDefaultBeamCheckbox = widgets.Checkbox(description = 'Show Default Beam', layout = widgets.Layout(width = '13%'), style = {'description_width': '0%'})

#Defines "Show displacement" checkbox.
showDisplacementCheckbox = widgets.Checkbox(description = 'Show Displacement', layout = widgets.Layout(width = '13%'), style = {'description_width': '0%'})

viewXYPlaneCheckbox = widgets.Checkbox(description = 'View x-y plane', layout = widgets.Layout(width = '13%'), style = {'description_width': '0%'})

#Defines Flag Font Size text box.
flagFontSizeTextBox = widgets.Text(value = '12', description = 'Flag Font Size', layout = widgets.Layout(width = '20%'), style = {'description_width': '50%'})

#Defines Plot Font Size text box.
#global plotFontSizeTextBox
plotFontSizeTextBox = widgets.Text(value = '12', description = 'Plot Font Size', layout = widgets.Layout(width = '20%'), style = {'description_width': '50%'})

#Defines Axes Ticks text box.
#global axesTicksTextBox
axesTicksTextBox = widgets.Text(value = '7', description = 'Axes Ticks', layout = widgets.Layout(width = '15%'))

#Defines Nummber of Circles text box
#global numberOfCirclesTextBox
numberOfCirclesTextBox = widgets.Text(value = '10', description = 'Number Of Circles', layout = widgets.Layout(width = '15%'), style = {'description_width': '70%'})

#Defines Nummber of Lines text box
#global numberOfLinesTextBox
numberOfLinesTextBox = widgets.Text(value = '16', description = 'Number Of Lines', layout = widgets.Layout(width = '15%'), style = {'description_width': '70%'})

#Defines Beam Color dropdown menu
#global beamColorDropdown
beamColorDropdown = widgets.Dropdown(options = ['Red', 'Green', 'Blue', 'Cyan', 'Magenta', 'Yellow'], value = 'Red', description = 'Beam Color', layout = widgets.Layout(width = '15%'))

#global mainLabel
mainLabel = widgets.Label(value = 'Initial Allignment Assistant', width = '100%')

#global optionsLabel
optionsLabel = widgets.Label(value = 'Formating Options', width = '100%')


