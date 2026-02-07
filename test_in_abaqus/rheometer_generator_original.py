# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2018 replay file
# Internal Version: 2017_11_07-17.21.41 127140
# Run by jferreira on Tue Apr 21 13:32:11 2020
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from driverUtils import executeOnCaeStartup
from caeModules import *
import os
from abaqus import *
from abaqusConstants import *
#material
kk = 1.0
k_bulk = 1000.0*kk
d_bulk = 1/(2*k_bulk)
# sweep angle
sweepangle = 360.0
#mesh seed
seed_size = 5.0
#disk dimensions
rd=100.0
t=10.0
h3=30.0
hb=10.0
#
#gel
h1 = 1.0
h2 = 2.0
#rotation
rotation = 0.1 #radians
#parallel computing
num = 16
# locations and directories
cwd = os.getcwd()
pth = cwd+'/'
output = 'gel_rheometer'
output_odb = output+'.odb'
output_cae = output+'.cae'
pthname = pth+output
routine = 'umat_.f'
pthrout = pth+routine
pthroutodb = pth+output_odb
pthoutput_cae = pth+output_cae
# user material
dsdv = 10
d = 0.0001
c = 1.32442
cc = 0.0
######################## GEOMETRY ########################
## BOTTOM ##
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=257.439575195312,
                 height=127.375930786133)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
                                            sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -rd), point2=(0.0, rd))
s.FixedConstraint(entity=g[2])
s.rectangle(point1=(0.0, 0.0), point2=(rd, -hb))
s.undo()
s.Line(point1=(0.0, 0.0), point2=(rd, 0.0))
s.HorizontalConstraint(entity=g[3], addUndoState=False)
s.Line(point1=(rd, 0.0), point2=(rd, -hb))
s.VerticalConstraint(entity=g[4], addUndoState=False)
s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
s.Line(point1=(rd, -hb), point2=(0.0, -hb))
s.HorizontalConstraint(entity=g[5], addUndoState=False)
s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D,
                               type=ANALYTIC_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['Part-1']
p.AnalyticRigidSurfRevolve(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
mdb.models['Model-1'].parts.changeKey(fromName='Part-1', toName='bottom')
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
                                             sheetSize=200.0)
## TOP ##
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.ConstructionLine(point1=(0.0, -rd), point2=(0.0, rd))
#s1.FixedConstraint(entity=g[2])
s1.Line(point1=(0.0, h1), point2=(t, h1))
s1.Line(point1=(t, h1), point2=(rd, h2))
s1.Line(point1=(rd, h2), point2=(rd, h2+hb))
#s1.VerticalConstraint(entity=g[4], addUndoState=False)
s1.Line(point1=(rd, h2+hb), point2=(h3, h2+hb))
#s1.HorizontalConstraint(entity=g[5], addUndoState=False)
#s1.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
s1.Line(point1=(h3, h2+hb), point2=(h3, h2+hb+h3))
#s1.VerticalConstraint(entity=g[6], addUndoState=False)
#s1.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
s1.Line(point1=(h3, h2+hb+h3), point2=(0.0, h2+hb+h3))
#s1.HorizontalConstraint(entity=g[7], addUndoState=False)
#s1.PerpendicularConstraint(entity1=g[6], entity2=g[7], addUndoState=False)
p = mdb.models['Model-1'].Part(name='top', dimensionality=THREE_D,
                               type=ANALYTIC_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['top']
p.AnalyticRigidSurfRevolve(sketch=s1)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['top']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=447.421,
                                                farPlane=791.349, width=351.084, height=148.918, cameraPosition=(311.862,
                                                                                                                 -48.5864, 531.471), cameraUpVector=(-0.383801, 0.920543, -0.0727827),
                                                cameraTarget=(5.84661, 3.30666, 5.84665))
p = mdb.models['Model-1'].parts['top']
## GEL ##
v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v1[0])
p = mdb.models['Model-1'].parts['bottom']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['bottom']
v2, e1, d2, n1 = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v2[1])
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
                                            sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -rd), point2=(0.0, rd))
s.FixedConstraint(entity=g[2])
s.Line(point1=(0.0, 0.0), point2=(rd, 0.0))
s.HorizontalConstraint(entity=g[3], addUndoState=False)
s.Line(point1=(rd, 0.0), point2=(rd, h2))
s.VerticalConstraint(entity=g[4], addUndoState=False)
s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
s.Line(point1=(rd, h2), point2=(t, h1))
s.Line(point1=(t, h1), point2=(0.0, h1))
s.Line(point1=(0.0, h1), point2=(0.0, 0.0))
s.VerticalConstraint(entity=g[7], addUndoState=False)
p = mdb.models['Model-1'].Part(name='gel', dimensionality=THREE_D,
                               type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['gel']
p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['gel']
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['gel']
c = p.cells
pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
v1, e, d1 = p.vertices, p.edges, p.datums
p.PartitionCellByPlaneThreePoints(point3=v1[0], cells=pickedCells,
                                  point1=p.InterestingPoint(
                                      edge=e[1], rule=MIDDLE),
                                  point2=p.InterestingPoint(edge=e[3], rule=MIDDLE))
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
mdb.saveAs(
    pathName=pthoutput_cae)
######################## MATERIAL ########################
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON,
                                                       engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
mdb.models['Model-1'].Material(name='Material-1')
##mdb.models['Model-1'].materials['Material-1'].Elastic(table=((rd, 0.45), ))
mdb.models['Model-1'].materials['Material-1'].Hyperelastic(materialType=ISOTROPIC,
                                                    testData=OFF, type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA, table=((
                                                        kk, d_bulk), ))
mdb.models['Model-1'].materials.changeKey(fromName='Material-1', toName='gel')
mdb.models['Model-1'].HomogeneousSolidSection(name='gel', material='gel',
                                              thickness=None)
p = mdb.models['Model-1'].parts['gel']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
region = p.Set(cells=cells, name='gel')
p = mdb.models['Model-1'].parts['gel']
p.SectionAssignment(region=region, sectionName='gel', offset=0.0,
                    offsetType=MIDDLE_SURFACE, offsetField='',
                    thicknessAssignment=FROM_SECTION)
p1 = mdb.models['Model-1'].parts['gel']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].parts['gel']
c = p.cells
cells = c.getSequenceFromMask(mask=('[#3 ]', ), )
p.Set(cells=cells, name='gel')
mdb.save()
######################## ASSEMBLY ########################
session.viewports['Viewport: 1'].view.setValues(width=345.381, height=146.499,
                                                viewOffsetX=-1.2352, viewOffsetY=-0.360678)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['gel']
a.Instance(name='gel-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['bottom']
a.Instance(name='bottom-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['top']
a.Instance(name='top-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('top-1', ), vector=(0.0, 0, 0.0))
######################## BOUNDARY CONDITIONS ########################
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].StaticStep(name='rotation', previous='Initial',
                                 nlgeom=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='rotation')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON,
                                                           constraints=ON, connectors=ON, engineeringFeatures=ON,
                                                           adaptiveMeshConstraints=OFF)
i1 = mdb.models['Model-1'].rootAssembly.allInstances['bottom-1']
leaf = dgm.LeafFromInstance(instances=(i1, ))
session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
    leaf=leaf)
i1 = mdb.models['Model-1'].rootAssembly.allInstances['gel-1']
leaf = dgm.LeafFromInstance(instances=(i1, ))
session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
    leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=451.107,
                                                farPlane=800.527, width=366.36, height=155.397, cameraPosition=(383.754,
                                                                                                                -169.651, 460.637), cameraUpVector=(-0.138174, 0.989958, 0.0298515),
                                                cameraTarget=(5.85443, -1.70882, 5.85442))
session.viewports['Viewport: 1'].view.setValues(nearPlane=458.547,
                                                farPlane=795.363, width=372.402, height=157.96, cameraPosition=(369.903,
                                                                                                                -244.048, 437.932), cameraUpVector=(-0.050425, 0.99135, 0.121173),
                                                cameraTarget=(5.60809, -3.03192, 5.45063))
 ######################## INTERACTION ########################
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['bottom-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#7 ]', ), )
region1 = a.Surface(side1Faces=side1Faces1, name='master_bottom')
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['gel-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#4 ]', ), )
region2 = a.Surface(side1Faces=side1Faces1, name='s_Surf-1')
mdb.models['Model-1'].Tie(name='tie_bottom', master=region1, slave=region2,
                          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=446.886,
                                                farPlane=796.673, width=362.932, height=153.943, cameraPosition=(438.768,
                                                                                                                 143.725, 419.837), cameraUpVector=(-0.462454, 0.836723, -0.293312),
                                                cameraTarget=(6.95566, 4.5561, 5.09654))
i1 = mdb.models['Model-1'].rootAssembly.allInstances['top-1']
leaf = dgm.LeafFromInstance(instances=(i1, ))
session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
    leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=454.864,
                                                farPlane=802.376, width=369.411, height=156.691, cameraPosition=(474.335,
                                                                                                                 -152.189, 373.025), cameraUpVector=(-0.172866, 0.981995, 0.0761738),
                                                cameraTarget=(7.16826, 2.78726, 4.81672))
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
i1 = mdb.models['Model-1'].rootAssembly.allInstances['gel-1']
leaf = dgm.LeafFromInstance(instances=(i1, ))
session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
    leaf=leaf)
session.viewports['Viewport: 1'].view.setValues(nearPlane=479.759,
                                                farPlane=767.954, width=389.629, height=165.267, cameraPosition=(619.713,
                                                                                                                 47.2131, 61.562), cameraUpVector=(-0.340371, 0.8027, -0.489715),
                                                cameraTarget=(9.53021, 6.02692, -0.243579))
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['top-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#1f ]', ), )
region1 = a.Surface(side1Faces=side1Faces1, name='master_top')
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['gel-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region2 = a.Surface(side1Faces=side1Faces1, name='slave_top')
mdb.models['Model-1'].Tie(name='tie_top', master=region1, slave=region2,
                          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON,
                                                           predefinedFields=ON, interactions=OFF, constraints=OFF,
                                                           engineeringFeatures=OFF)
leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
    leaf=leaf)
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['bottom-1'].referencePoints
refPoints1 = (r1[2], )
region = a.Set(referencePoints=refPoints1, name='bottom_rp')
mdb.models['Model-1'].DisplacementBC(name='fix_bottom',
                                     createStepName='rotation', region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0,
                                     ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
                                     fieldName='', localCsys=None)
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['top-1'].referencePoints
refPoints1 = (r1[2], )
region = a.Set(referencePoints=refPoints1, name='top_rp')
mdb.models['Model-1'].DisplacementBC(name='top_rotation',
                                     createStepName='rotation', region=region, u1=UNSET, u2=UNSET, u3=UNSET,
                                     ur1=UNSET, ur2=0.1, ur3=UNSET, amplitude=UNSET, fixed=OFF,
                                     distributionType=UNIFORM, fieldName='', localCsys=None)
mdb.models['Model-1'].boundaryConditions['fix_bottom'].move('rotation',
                                                            'Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF,
                                                           bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['gel']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF,
                                                       engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['Model-1'].parts['top']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['gel']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['bottom']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['gel']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['gel']
p.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['gel']
p.generateMesh()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, loads=ON,
                                                           bcs=ON, predefinedFields=ON, connectors=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.save()
#: The model database has been saved
mdb.Job(name=output, model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=2,
        numDomains=2, numGPUs=0)
