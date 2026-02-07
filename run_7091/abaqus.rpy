# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2025 replay file
# Internal Version: 2024_09_20-14.00.46 RELr427 198590
# Run by lpacheco on Mon Oct 20 17:47:44 2025
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=312.472900390625, 
    height=207.680557250977)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o1 = session.openOdb(
    name='/home/lpacheco/biopolymer-networks/3_affnet_cl_ai_mcs/test_in_abaqus/runs_rheometer_03_nemix/run_7091/gel_rheometer.odb', 
    readOnly=False)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/lpacheco/biopolymer-networks/3_affnet_cl_ai_mcs/test_in_abaqus/runs_rheometer_03_nemix/run_7091/gel_rheometer.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     3
#: Number of Meshes:             3
#: Number of Element Sets:       3
#: Number of Node Sets:          6
#: Number of Steps:              1
leaf = dgo.LeafFromPartInstance(partInstanceName=('GEL-1', ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.replace(leaf=leaf)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=NONE)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=FREE)
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap=session.viewports['Viewport: 1'].colorMappings['Part instance']
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].enableMultipleColors()
session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
cmap = session.viewports['Viewport: 1'].colorMappings['Part instance']
cmap.updateOverrides(overrides={'GEL-1':(True, '#008080', 'Default', 
    '#008080')})
session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
session.viewports['Viewport: 1'].disableMultipleColors()
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(session.views['Bottom'])
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=33147, 
    farPlane=48094.5, width=15938.3, height=7794.54, cameraPosition=(30803.6, 
    34340.5, -1792.6), cameraUpVector=(-0.104356, -0.103626, -0.989127), 
    cameraTarget=(9255.89, -4334.65, 4532.56))
session.viewports['Viewport: 1'].view.setValues(nearPlane=33603.2, 
    farPlane=47638.3, width=16157.6, height=7901.81, cameraPosition=(31062.6, 
    34045.4, -2714.38), cameraUpVector=(-0.181826, -0.0592242, -0.981546), 
    cameraTarget=(9514.88, -4629.7, 3610.78))
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(nearPlane=42764.3, 
    farPlane=46680.2, width=20562.7, height=10056.1, cameraPosition=(8203.15, 
    44792.2, -7596.6), cameraUpVector=(-0.935854, 0, -0.352388), cameraTarget=(
    8203.15, 70, -7596.6))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(54722.2, 
    70.0039, 5000), cameraUpVector=(0, 1, 0))
session.viewports['Viewport: 1'].view.setValues(session.views['User-1'])
session.viewports['Viewport: 1'].view.setValues(session.views['User-2'])
session.viewports['Viewport: 1'].view.setValues(session.views['User-3'])
session.viewports['Viewport: 1'].view.setValues(session.views['User-4'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    referencePoints=OFF)
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    patternRotationAxis=YAXIS, patternNumCircular=2)
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    referencePoints=ON, mirrorCsysName=GLOBAL, patternCsysName=GLOBAL, 
    patternRotationAxis=ZAXIS, patternNumCircular=1)
session.viewports['Viewport: 1'].view.setValues(nearPlane=42764.3, 
    farPlane=46680.1, width=20562.7, height=10056.1, cameraPosition=(7718.76, 
    44792.2, -8088.31), cameraUpVector=(-0.955815, 0, -0.293969), 
    cameraTarget=(7718.76, 70, -8088.31))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(10000, 70, 
    49722.2), cameraUpVector=(0, 1, 0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=36336.5, 
    farPlane=53107.9, width=38536.2, height=18845.9, cameraPosition=(65.6739, 
    10000, 49722.2), cameraUpVector=(-0.999908, 0.0135667, 0), cameraTarget=(
    65.6739, 10000, 5000))
session.viewports['Viewport: 1'].view.setValues(nearPlane=45894.5, 
    farPlane=53423.3, width=48672.8, height=23803.1, cameraPosition=(-427.922, 
    49940.4, -8836.01), cameraUpVector=(-0.999599, -0.011726, -0.0257862), 
    cameraTarget=(123.424, 5233.76, -9879.11))
session.viewports['Viewport: 1'].view.setValues(cameraPosition=(44845.6, 
    5233.76, -9879.11), cameraUpVector=(0, 1, 0))
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    translucency=ON)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    translucency=OFF)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(session.views['Bottom'])
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=OFF, 
    legend=OFF, title=OFF, state=OFF, annotations=OFF, compass=OFF)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=NONE)
session.printOptions.setValues(vpDecorations=OFF)
session.printToFile(fileName='final1_0', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='final1_0', format=SVG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='final1_0', format=EPS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=ON)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=OFF)
session.printToFile(fileName='final1_1', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='final1_1', format=SVG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='final1_1', format=EPS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(
    patternRotationAxis=YAXIS, patternNumCircular=12)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=EXTERIOR, edgeColorFillShade='#000000')
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    edgeColorFillShade='Black')
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=NONE)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=EXTERIOR)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=NONE)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.fitView()
session.printToFile(fileName='final2_0', format=PNG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='final2_0', format=SVG, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.printToFile(fileName='final2_0', format=EPS, canvasObjects=(
    session.viewports['Viewport: 1'], ))
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=ON)
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(session.views['Back'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].view.setValues(session.views['Top'])
