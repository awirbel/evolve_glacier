from __future__ import division, print_function
from osgeo import gdal
'''
    Create geometry input files
    
    Evolve
    
    Copyright (C) 2019  Alexander Jarosch 
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

''' open a geotif and sample '''
# define the artificial ice Layer
ail = 1.0

# set up the geotiff stuff
driver = gdal.GetDriverByName('GTiff')

surf_filename = 'surf.tif'
surf_dataset = gdal.Open(surf_filename)
surf_band = surf_dataset.GetRasterBand(1)
surf_cols = surf_dataset.RasterXSize
surf_rows = surf_dataset.RasterYSize
surf_transform = surf_dataset.GetGeoTransform()
surf_xOrigin = surf_transform[0]
surf_yOrigin = surf_transform[3]
surf_pixelWidth = surf_transform[1]
surf_pixelHeight = -surf_transform[5]
surf_data = surf_band.ReadAsArray(0, 0, surf_cols, surf_rows)

bed_filename = 'bed.tif'
bed_dataset = gdal.Open(bed_filename)
bed_band = bed_dataset.GetRasterBand(1)
bed_cols = bed_dataset.RasterXSize
bed_rows = bed_dataset.RasterYSize
bed_transform = bed_dataset.GetGeoTransform()
bed_xOrigin = bed_transform[0]
bed_yOrigin = bed_transform[3]
bed_pixelWidth = bed_transform[1]
bed_pixelHeight = -bed_transform[5]
bed_data = bed_band.ReadAsArray(0, 0, bed_cols, bed_rows)

# open the input and output files
filepath = 'surface_shell_domain.stl'
filepath_out = 'glacier_start.stl'
fp = open(filepath, 'r')
fp_out = open(filepath_out, 'w')

for line in fp:
    linesplit = line.split()
    if linesplit[0] == 'vertex':
        if linesplit[3] == '0.000':
            xcoord = float(linesplit[1])
            ycoord = float(linesplit[2])
            point = (xcoord, ycoord)
            # get the z coord
            scol = int((point[0] - surf_xOrigin) / surf_pixelWidth)
            srow = int((surf_yOrigin - point[1]) / surf_pixelHeight)
            zcoord = surf_data[srow][scol]
            fp_out.write("    vertex %.03f %.03f %.03f\n" % (xcoord, ycoord,
                                                       zcoord))

        elif linesplit[3] == '-1000.000':
            xcoord = float(linesplit[1])
            ycoord = float(linesplit[2])
            point = (xcoord, ycoord)
            # get the z coord
            bcol = int((point[0] - bed_xOrigin) / bed_pixelWidth)
            brow = int((bed_yOrigin - point[1]) / bed_pixelHeight)
            zcoord = bed_data[brow][bcol] - ail
            fp_out.write("    vertex %.03f %.03f %.03f\n" % (xcoord, ycoord,
                                                       zcoord))

    else:
        fp_out.write(line)

fp_out.close()
fp.close()

print('all done')
