
import re
from iris.unit import Unit

import matplotlib.pyplot as plt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import iris
import iris.plot as iplt

#
# {GRIB_ID: (standard_name, long_name, units)
#
standard_name_map = {
    1: ('atmosphere_horizontal_streamfunction', None, 'm2 s-1'), # Stream function, strf, m2 s-1
    34: ('sea_surface_temperature', None, 'K'), # Sea surface temperature, sst, K
    59: ('atmosphere_specific_convective_available_potential_energy', None, 'J kg-1'), # Convective available potential energy, cape, J kg-1
    60: ('ertel_potential_vorticity', None, 'K m2 kg-1 s-1'),
    130: ('air_temperature', None, 'K'), # Temperature, t, K
    131: ('eastward_wind', None, 'm s-1'), # U component of wind, u, m s-1
    132: ('northward_wind', None, 'm s-1'), # V component of wind, v, m s-1
    133: ('specific_humidity', None, 'kg kg-1'), # Specific humidity, q, kj kj-1
    135: ('lagrangian_tendency_of_air_pressure', None, 'Pa s-1'), # Vertical velocity, w, Pa s-1
    138: ('atmosphere_relative_vorticity', None, 's-1'), # Vorticity (relative), vo, s-1
    142: ('lwe_thickness_of_large_scale_precipitation_amount', None, 'm'), # Large-scale precipitation, lsp, m
    143: ('lwe_thickness_of_convective_precipitation_amount', None, 'm'), # Convective precipitation, cp, m
    155: ('divergence_of_wind', None, 's-1'), # Divergence, d, s-1
    156: ('geopotential_height', None, None), # Geopotential Height, gh, gpm
    157: ('relative_humidity', None, '%'), # Relative humidity, r, %
    159: ('atmospheric_boundary_layer_thickness', None, 'm'), # Boundary layer height, blh, m
    186: ('cloud_area_fraction', None, '1'), # Low cloud cover, lcc, (0-1)
    151: ('air_pressure_at_sea_level', None, 'Pa'), # Top net radiation, tnr, W m-2
    165: (None, '10 METER U WIND COMPONENT', 'm s-1'), # 10 meter U wind component, 10u, m s-1
    166: (None, '10 METER V WIND COMPONENT', 'm s-1'), # 10 meter V wind component, 10v, m s-1
    167: (None, '2 METER TEMPERATURE', 'K'), # 2 meter temperature, 2t, K
    168: (None, '2 METER DEWPOINT TEMPERATURE', 'K'), # 2 meter dwepoint temperature, 2d, K
    187: (None, 'MEDIUM CLOUD COVER', '1'), # Medium cloud cover, lcc, (0-1)
    188: (None, 'HIGH CLOUD COVER', '1'), # High cloud cover, hcc, (0-1)
    206: (None, 'TOTAL COLUMN OZONE', 'kg m-2'), # Total column ozone, tco3, kg m-2
    228: (None, 'HUMIDITY TENDENCY BY LARGE SCALE CONDENSATION', 'kg kg-1'), # Humidity tendency by large-scale condensation, htlc, kg kg-1
}

param_rexp = re.compile(r'UNKNOWN LOCAL PARAM (\d+)\.(\d+)')

def _metadata_callback(cube, field, filename):
    match = param_rexp.match(cube.long_name)
    if match:
        parameter, table_version = (int(x) for x in match.groups())
        if parameter in standard_name_map:
            standard_name, long_name, unit = standard_name_map[parameter]
            if standard_name:
                cube.long_name = standard_name
            if long_name:
                cube.long_name = long_name
            if unit:
                cube.units = Unit(unit)



def load_cubes(filename):
    cubes = iris.load(filename, callback=_metadata_callback)

    return cubes


def plot_cube(cube, pressure_level=None):

    # Detect whether the cube is 2D or 3D
    coord_names = [d.name() for d in cube.dim_coords]
    if len(coord_names) == 3:
        assert 'pressure' in coord_names
        is_3d = True
    else:
        assert len(coord_names) == 2
        is_3d = False

    # Select slice
    if is_3d:
        if pressure_level:
            p = cube.coord('pressure')
            cslice = cube.subset(p[p==pressure_level])[0]
        else:
            cslice = cube[0]
            level = cube.coord('pressure').points[0]
    else:
        cslice = cube
        level = 'surface'
    

    #plot = qplt.contourf(cslice)
    #ax = plot.ax

    ax = plt.axes([0.1, 0.1, 0.6, 0.6], projection=ccrs.PlateCarree())

    ax.set_aspect('auto')
    ax.set_xlim(80, 240)
    ax.set_ylim(-40, 40)
    ax.gridlines(draw_labels=True)
    ax.coastlines()
    iplt.contourf(cslice, axes=ax)

    #!TODO: adjust viewport
    ax.autoscale(False)

    return ax
