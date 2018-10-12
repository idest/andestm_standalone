import pandas as pd
import numpy as np
import shapely
from shapely.geometry import Point, LineString, Polygon
from descartes.patch import PolygonPatch
from pyproj import Proj
from utils import module_from_file
import matplotlib.pyplot as plt
#from src import setup
#from src import compute

compute = module_from_file('src', 'src/compute.py')
setup = module_from_file('src', 'src/setup.py')

data = compute.Data(*setup.data_setup(), *setup.input_setup())
cs = compute.CoordinateSystem(data.get_cs_data(), 0.2, 1)
gm = compute.GeometricModel(data.get_gm_data(), cs)

#eq1 = pd.read_excel('data/earthquakes/CSN.xlsx', sheet_name='Sheet1')
#eq2 = pd.read_excel('data/earthquakes/CSN.xlsx', sheet_name='Sheet2')

utm = Proj(proj='utm', zone=19, ellps='WGS84')
for idx, lat in enumerate(cs.get_y_axis()[:-1]):
    df = pd.DataFrame(
        {'lon': cs.get_x_axis().copy(),
         'topo': gm.get_topo().cross_section(latitude=lat).copy(),
         'slab_lab': gm.get_slab_lab().cross_section(latitude=lat).copy()}
    )
    df['lon'] = df['lon'].apply(lambda lon: utm(lon,lat)[0])
    first_lon = df['lon'].iloc[0]
    df['lon'] = df['lon'].apply(lambda lon: abs(first_lon-lon))
    # Remove nan values at margins
    first_valid = max(df['topo'].first_valid_index(), df['slab_lab'].first_valid_index())
    last_valid = min(df['topo'].last_valid_index(), df['slab_lab'].last_valid_index())
    df = df[first_valid:last_valid+1].reset_index(drop=True)
    # Remove nan values at middle of slab_lab
    df['slab_lab'] = df['slab_lab'].interpolate(method='nearest')
    # Remove values to the left of last intersecion between topo and slab_lab
    intersections = df['topo'] == df['slab_lab']
    if intersections.any():
        last_intersection = intersections.index[intersections][-1]
        df = df[last_intersection:].reset_index(drop=True)
    # Create shapely polygon
    polygon_x_points = df['lon'].append(df['lon'][::-1])
    polygon_y_points = df['topo'].append(df['slab_lab'][::-1])
    polygon = Polygon(zip(polygon_x_points,polygon_y_points))
    print(polygon.is_valid)
    patch = PolygonPatch(polygon)
    # Plot
    if True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(df['lon'], df['topo'])
        ax.plot(df['lon'], df['slab_lab'])
        ax.add_patch(patch)
        print(lat)
        plt.show()
        #break
    #a = df['topo'] == df['slab_lab']
    #if a.any() is True:
    #    print(a)
    #print('###')
    #print(df['slab_lab'].shape)
    #print(df['topo'].shape)
    #print(df['lon'].shape)
    #print('###')
    #print(topo)
