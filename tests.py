#SHF Map plot test
import os
from plot2 import shf_map

def shf_map_test(shf, data, diff, data_coords, data_types, rmse):
    save_dir = 'Tests/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    # base map
    shf_map(save_dir=save_dir, name='map_base')
    # shf
    shf_map(shf=shf,
            rmse=rmse,
            save_dir=save_dir, name='map_shf')
    # data
    shf_map(data=data,
            data_coords=data_coords, data_types=data_types, rmse=rmse,
            save_dir=save_dir, name='map_data')
    # diff
    shf_map(diff=diff,
            data_coords=data_coords, data_types=data_types, rmse=rmse,
            save_dir=save_dir, name='map_diff')
    # shf + data
    shf_map(shf=shf, data=data,
            data_coords=data_coords, data_types=data_types, rmse=rmse,
            save_dir=save_dir, name='map_shf_data')
    # shf + diff
    shf_map(shf=shf, diff=diff,
            data_coords=data_coords, data_types=data_types, rmse=rmse,
            save_dir=save_dir, name='map_shf_diff')
    # data + diff
    shf_map(data=data, diff=diff,
            data_coords=data_coords, data_types=data_types, rmse=rmse,
            save_dir=save_dir, name='map_data_diff')
    # shf + data + diff
    shf_map(shf=shf, data=data, diff=diff,
            data_coords=data_coords, data_types=data_types, rmse=rmse,
            save_dir=save_dir, name='map_shf_data_diff')
