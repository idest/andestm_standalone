import numpy as np
import matplotlib.cm as cm
from src.colormaps import jet_white_r_2

colors = jet_white_r_2(np.arange(256))
print(colors)
name = 'jet_white_r_2'
path = 'vtk_files/jet_white_r_2.json'

def writePyPlotCMapToParaView(colors,name,path):
    N = len(colors)
    fid = open(path, 'w')
    fid.write('[\n')
    fid.write('  {\n')
    fid.write('    "ColorSpace" : "RGB",\n')
    fid.write('    "Name" : "{}",\n'.format(name))
    fid.write('    "RGBPoints" : \n')
    fid.write('    [')
    for i in range(N):
        fid.write('      {},\n'.format(-1 + 2*(i/N)))
        fid.write('      {},\n'.format(colors[i,0]))
        fid.write('      {},\n'.format(colors[i,1]))
        fid.write('      {},\n'.format(colors[i,2]))
    fid.write('    ]\n')
    fid.write('  }\n')
    fid.write(']\n')
    fid.close()

writePyPlotCMapToParaView(colors, name, path)
