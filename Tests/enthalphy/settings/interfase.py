#### import the simple module from the paraview
from paraview.simple import *

Source = 'results/graph/graph.pvd'
End = 'results/data.txt'

# create a new 'PVD Reader'
Data = PVDReader(registrationName='Field', FileName=Source)
Data.CellArrays = ['attribute']
Data.PointArrays = ['Temperature']
Data.ColumnArrays = []

# create a new 'Contour'
Phase = Contour(registrationName='Phase', Input=Data)
Phase.ContourBy = ['POINTS', 'Temperature']
Phase.ComputeNormals = 1
Phase.ComputeGradients = 0
Phase.ComputeScalars = 1
Phase.OutputPointsPrecision = 'Same as input'
Phase.GenerateTriangles = 1
Phase.Isosurfaces = [0.0]
Phase.PointMergeMethod = 'Uniform Binning'
Phase.PointMergeMethod.Divisions = [50, 50, 50]
Phase.PointMergeMethod.Numberofpointsperbucket = 8

# create a new 'Plot Data Over Time'
Plot = PlotDataOverTime(registrationName='Plot', Input=Phase)
Plot.FieldAssociation = 'Points'
Plot.OnlyReportSelectionStatistics = 1

# save data
SaveData(End, proxy=Plot, WriteTimeSteps=0,
    Filenamesuffix='_%d',
    ChooseArraysToWrite=1,
    PointDataArrays=[],
    CellDataArrays=[],
    FieldDataArrays=[],
    VertexDataArrays=[],
    EdgeDataArrays=[],
    RowDataArrays=['Time', 'avg(Y)', 'std(Y)'],
    Precision=5,
    UseScientificNotation=0,
    FieldAssociation='Row Data',
    AddMetaData=1,
    AddTime=0)
