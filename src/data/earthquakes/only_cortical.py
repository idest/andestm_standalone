import pandas as pd
import numpy as np
import shapely
from shapely.geometry import Point, LineString, Polygon, MultiPoint
from descartes.patch import PolygonPatch
from pyproj import Proj
from utils import module_from_file
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
#from src import setup
#from src import compute

compute = module_from_file('src', 'src/compute.py')
setup = module_from_file('src', 'src/setup.py')

data = compute.Data(*setup.data_setup(), *setup.input_setup())
cs = compute.CoordinateSystem(data.get_cs_data(), 0.2, 1)
gm = compute.GeometricModel(data.get_gm_data(), cs)

#Earthquakes CSN
#eqs_csn = pd.read_excel('data/earthquakes/CSN_2000_2018.xlsx', sheet_name='Sheet1')
#Earthquakes USGS
eqs_usgs = pd.read_csv('data/earthquakes/USGS_1900_2018.csv')

eqs_new = pd.DataFrame()

utm = Proj(proj='utm', zone=19, ellps='WGS84')
for idx, lat in enumerate(cs.get_y_axis()[:-1]):
    df = pd.DataFrame(
        {'lon': cs.get_x_axis().copy(),
         'topo': gm.get_topo().cross_section(latitude=lat).copy(),
         'slab_lab': gm.get_slab_lab().cross_section(latitude=lat).copy()}
    )
    df['sli'] = False
    df['sli'].iloc[gm.get_slab_lab_int_index()[idx]] = True
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
    # SA Plate Polygon 
    polygon_x_points = df['lon'].append(df['lon'][::-1])
    polygon_y_points = df['topo'].append(df['slab_lab'][::-1])
    polygon = Polygon(zip(polygon_x_points,polygon_y_points))
    #print(polygon.is_valid)
    patch = PolygonPatch(polygon)
    # SLAB Line Buffer
    line_x_points = np.asarray(
        df['lon'].iloc[:np.argwhere(df['sli'] == True)[0][0]+1])
    line_y_points = np.asarray(
        df['slab_lab'].iloc[:np.argwhere(df['sli'] == True)[0][0]+1])
    line_last_slope = ((line_y_points[-1] - line_y_points[-5])
        /(line_x_points[-1] - line_x_points[-5]))
    line_last_y = line_y_points[-1] - 20.
    line_last_x = ((line_last_y - line_y_points[-1]) 
            + line_last_slope*line_x_points[-1])/line_last_slope
    line_y_points = np.append(line_y_points, line_last_y)
    line_x_points = np.append(line_x_points, line_last_x)
    line = LineString(zip(line_x_points, line_y_points))
    lineBuffer = line.buffer(20, cap_style=1, mitre_limit=100)
    patch2 = PolygonPatch(lineBuffer, fc='green')

    # Earthquakes CSN
    #eq = eqs_csn[
    #    (eqs_csn['Lat.'] >= lat - 0.1 ) &
    #    (eqs_csn['Lat.'] < lat + 0.1)]
    #eq['x'] = eq['Lon.'].apply(lambda lon: abs(first_lon - utm(lon,lat)[0]))
    #eq['ISA'] = False
    #eq['OSB'] = True
    #print(lat)
    #for i in np.arange(len(eq.index)):
    #    p = Point(eq['x'].iloc[i], -eq['Prof.'].iloc[i])
    #    #print(polygon.contains(p))
    #    # Inside South America
    #    eq['ISA'].iloc[i] = polygon.contains(p)
    #    # Outside Slab Buffer
    #    eq['OSB'].iloc[i] = not lineBuffer.contains(p)
    #eqs_new = eqs_new.append(eq)

    #Earthquakes USGS
    eq = eqs_usgs[
        (eqs_usgs['latitude'] >= lat - 0.1) &
        (eqs_usgs['latitude'] < lat + 0.1)]
    eq['x'] = eq['longitude'].apply(lambda lon: abs(first_lon - utm(lon,lat)[0]))
    eq['ISA'] = False
    eq['OSB'] = True
    print(lat)
    for i in np.arange(len(eq.index)):
        p = Point(eq['x'].iloc[i], -eq['depth'].iloc[i])
        #print(polygon.contains(p))
        # Inside South America
        eq['ISA'].iloc[i] = polygon.contains(p)
        # Outside Slab Buffer
        eq['OSB'].iloc[i] = not lineBuffer.contains(p)
    eqs_new = eqs_new.append(eq)

"""
    # Plot
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Lat {:.1f}'.format(lat))
        ax.plot(df['lon'], df['topo'])
        ax.plot(df['lon'], df['slab_lab'])
        ax.add_patch(patch)
        ax.add_patch(patch2)
        eq['color'] = 'red'
        #eq['color'] = eq['color'].where(eq['ISA'] == True, other='black')
        # Earthquakes CSN
        depth = -eq['Prof.']
        mag = eq['m1mag']
        # Earthquakes USGS
        #depth = -eq['depth']
        #mag = eq['mag']
        #ax.scatter(
        #    eq['x'], depth, color=eq['color'], s=mag**2.+20., zorder=-1000)
        ax.scatter(eq['x'], depth, color=eq['color'], s=1., zorder=1000)
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
"""

# Earthquakes CSN
#writer = pd.ExcelWriter('data/earthquakes/CSN_2000_2018_C.xlsx')
# Earthquakes USGS
writer = pd.ExcelWriter('data/earthquakes/USGS_1900_2018_C.xlsx')
eqs_new.drop(columns=['x'])
eqs_new = eqs_new.sort_index()
eqs_new.to_excel(writer, 'Sheet1')
writer.save()

