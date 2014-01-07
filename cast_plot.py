
import re

import matplotlib
matplotlib.use('Agg')

from iris.unit import Unit

import matplotlib.pyplot as plt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import iris
import iris.plot as iplt

import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

import logging
log = logging.getLogger(__name__)

import os, sys, getopt

class GribMapping(object):
    def __init__(self):
        self._grib_id_map = {}

    def add(self, grib_id, standard_name, long_name, units):
        self._grib_id_map[grib_id] = (standard_name, long_name, units)
    
    def _metadata_callback(self, cube, field, filename):
        # Use IRIS values when present
        if cube.standard_name:
            return

        match = re.match(r'UNKNOWN LOCAL PARAM (\d+)\.(\d+)', cube.long_name)
        if match:
            parameter, table_version = (int(x) for x in match.groups())
            log.debug('Cube "{0}" --> {1}.{2}'.format(cube.long_name, parameter, table_version))
            if parameter in self._grib_id_map:
                standard_name, long_name, unit = self._grib_id_map[parameter]
                if standard_name:
                    cube.long_name = standard_name
                if long_name:
                    cube.long_name = long_name
                if unit:
                    cube.units = Unit(unit)
            else:
                log.error('Grib Parameter {0} not in grib mapping'.format(parameter))

    def load_cubes(self, filename):
        cubes = iris.load(filename, callback=self._metadata_callback)
        
        return cubes


grib_mapping = GribMapping()

grib_mapping.add(1, 'atmosphere_horizontal_streamfunction', None, 'm2 s-1') # Stream function, strf, m2 s-1
grib_mapping.add(34, 'sea_surface_temperature', None, 'K') # Sea surface temperature, sst, K
grib_mapping.add(59, 'atmosphere_specific_convective_available_potential_energy', None, 'J kg-1') # Convective available potential energy, cape, J kg-1
grib_mapping.add(60, 'ertel_potential_vorticity', None, 'K m2 kg-1 s-1')
grib_mapping.add(130, 'air_temperature', None, 'K') # Temperature, t, K
grib_mapping.add(131, 'eastward_wind', None, 'm s-1') # U component of wind, u, m s-1
grib_mapping.add(132, 'northward_wind', None, 'm s-1') # V component of wind, v, m s-1
grib_mapping.add(133, 'specific_humidity', None, 'kg kg-1') # Specific humidity, q, kj kj-1
grib_mapping.add(135, 'lagrangian_tendency_of_air_pressure', None, 'Pa s-1') # Vertical velocity, w, Pa s-1
grib_mapping.add(138, 'atmosphere_relative_vorticity', None, 's-1') # Vorticity (relative), vo, s-1
grib_mapping.add(142, 'lwe_thickness_of_large_scale_precipitation_amount', None, 'm') # Large-scale precipitation, lsp, m
grib_mapping.add(143, 'lwe_thickness_of_convective_precipitation_amount', None, 'm') # Convective precipitation, cp, m
grib_mapping.add(155, 'divergence_of_wind', None, 's-1') # Divergence, d, s-1
grib_mapping.add(156, 'geopotential_height', None, None) # Geopotential Height, gh, gpm
grib_mapping.add(157, 'relative_humidity', None, '%') # Relative humidity, r, %
grib_mapping.add(159, 'atmospheric_boundary_layer_thickness', None, 'm') # Boundary layer height, blh, m
grib_mapping.add(186, 'cloud_area_fraction', None, '1') # Low cloud cover, lcc, (0-1)
grib_mapping.add(151, 'air_pressure_at_sea_level', None, 'Pa') # Top net radiation, tnr, W m-2
grib_mapping.add(165, None, '10 METER U WIND COMPONENT', 'm s-1') # 10 meter U wind component, 10u, m s-1
grib_mapping.add(166, None, '10 METER V WIND COMPONENT', 'm s-1') # 10 meter V wind component, 10v, m s-1
grib_mapping.add(167, None, '2 METER TEMPERATURE', 'K') # 2 meter temperature, 2t, K
grib_mapping.add(168, None, '2 METER DEWPOINT TEMPERATURE', 'K') # 2 meter dwepoint temperature, 2d, K
grib_mapping.add(187, None, 'MEDIUM CLOUD COVER', '1') # Medium cloud cover, lcc, (0-1)
grib_mapping.add(188, None, 'HIGH CLOUD COVER', '1') # High cloud cover, hcc, (0-1)
grib_mapping.add(206, None, 'TOTAL COLUMN OZONE', 'kg m-2') # Total column ozone, tco3, kg m-2
grib_mapping.add(228, None, 'HUMIDITY TENDENCY BY LARGE SCALE CONDENSATION', 'kg kg-1') # Humidity tendency by large-scale condensation, htlc, kg kg-1





def plot_and_save(cube, out_filename, pressure_level=None, format=None):
    fig = plot_cube(cube, pressure_level)

    # Reasonable values for balancing plot and colourbar
    fig.set_size_inches(10, 7)
    fig.savefig(out_filename, format=format)
    log.info('Saving figure %s' % out_filename)


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
        p = cube.coord('pressure')
        if pressure_level:
            cslice = cube.subset(p[p==pressure_level])[0]
            level = pressure_level
        else:
            cslice = cube[0]
            level = p.points[0]
    else:
        cslice = cube
        level = 'surface'
        

    fig = plt.figure()
    
    log.info('Preparing plot for %s at level %s' % (cslice.name(), level))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines()
    fig.add_axes(ax)
             
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator([100, 140, 180, -140])
    
    contours = iplt.contourf(cslice, axes=ax)
    cb = plt.colorbar(contours, orientation='horizontal')
    cb.ax.set_xlabel(cslice.units)
    cb.ax.set_aspect(0.03)
             
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    ax.set_title('%s (level: %s)' % (cslice.name(), level))
             
    return fig


CAST_STANDARD_NAMES = [
    'relative_humidity',
    'air_temperature',
    'atmosphere_relative_vorticity',
    'geopotential_height',
    'eastward_wind',
    'northward_wind',
    'x_wind',
    'y_wind',
]
CAST_LONG_NAMES = [
  '10 METER U WIND COMPONENT',
  '10 METER V WIND COMPONENT',
  '2 METER TEMPERATURE',
]

def plot_cast_file(filename, outdir='.'):
    cubes = grib_mapping.load_cubes(filename)

    for i, cube in enumerate(cubes):
        if ((cube.standard_name not in CAST_STANDARD_NAMES) and
            (cube.long_name not in CAST_LONG_NAMES)):
            log.info('Skipping cube {0} ({1})'.format(cube.name(), cube.long_name))
            continue

        plot_name = '{0}_{1}.png'.format(cube.name().lower().replace(' ', '_'), 
                                         i)

        #!TODO: Configure pressure level
        fig = plot_and_save(cube, os.path.join(outdir,plot_name))

def usage():
    print '''usage: python2.7 cast_plot.py [-h | --help] [-d | --debug] [-o | --outdir OUTDIR ] file
    positional arguments:
        file        A GRIB file to be processed

    Optional arguments:
        -h, --help              show this help message and exit
        -d, --debug             show debug logs
        -o, --outdir OUTDIR     Output directory. Must be writable by the current user.
    '''


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],'hdo:',['help','debug','outdir='])
    except getopt.GetoptError as err:
        print str(err);
        usage()
        sys.exit(2)
    #defaults
    outdir = '.'
    level=logging.WARNING
   
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-d", "--debug"):
            level=logging.DEBUG 
        elif o in ("-o", "--outdir"):
            outdir = a
        else:
            assert False, "unhandled option" 
    
    logging.basicConfig(level=level)
    #os.access is technically vulnerable to TOCTOU errors, but...
    if(not os.access(a, os.W_OK)):
        log.error("directory '" + a +"' is not writable")
        sys.exit(2)
    cubes_file, = args[0:]

    plot_cast_file(cubes_file, outdir)
