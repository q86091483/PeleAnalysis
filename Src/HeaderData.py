#%%
import os
import shutil
import linecache
import numpy as np
from tqdm import tqdm

class HeaderData(object):

    def __init__(self, plotfile, limit_level=None):
        """
        Parse the header data and save as attributes
        """
        self.pfile = plotfile
        filepath = os.path.join(plotfile, 'Header')
        with open(filepath) as hfile:
            # Hyperclaw version
            self.version = hfile.readline().strip()
            # Field dict, e.g., fields["x_velocity"] = 0 
            self.nvars = int(hfile.readline())
            self.fields = {}
            for i in range(self.nvars):
                self.fields[hfile.readline().replace('\n', '')] = i
            # Space dimension
            self.ndims = int(hfile.readline())
            # Current time
            self.time = float(hfile.readline())
            # Low and high end of domain
            self.max_level = int(hfile.readline())
            self.geo_low = [float(n) for n in hfile.readline().split()]
            self.geo_high = [float(n) for n in hfile.readline().split()]
            # Refinement between level i/i+1 = factors[i]
            self.refinement_factors = [int(n) for n in hfile.readline().split()]
            # Domain box
            self.grid_sizes = []
            str = hfile.readline() 
            for block in str.split()[1::3]:
                grid_size = np.array(block.replace('(', '').replace(")", '').split(','), dtype=int)
                self.grid_sizes.append(grid_size + 1)
            self.box_types = []
            for block in str.split()[2::3]:
                box_type = np.array(block.replace('(', '').replace(")", '').split(','), dtype=int)
                self.box_types.append(box_type)
            # Step #
            self.nsteps = [int(n) for n in hfile.readline().split()]
            # Grid resolutions
            resolutions = []
            for i in range(self.max_level + 1):
                resolutions.append([float(n) for n in hfile.readline().split()])
            self.dx = resolutions
            # Coordinate system
            self.sys_coord = hfile.readline()
            # Sanity check
            assert 0 == int(hfile.readline())

print("test")
fn = "/scratch/b/bsavard/zisen347/scopingRuns/Burke9_Re6000_2J6_nonreactive/plt_07600"
dh = HeaderData(fn)
#%%