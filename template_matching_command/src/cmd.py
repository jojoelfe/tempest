# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc      # Command description
from chimerax.atomic import AtomsArg            # Collection of atoms argument
from chimerax.core.commands import BoolArg      # Boolean argument
from chimerax.core.commands import FileNameArg  
from chimerax.core.commands import FloatArg
from chimerax.core.commands import ModelArg
from chimerax.core.commands import ColorArg     # Color argument
from chimerax.core.commands import IntArg       # Integer argument
from chimerax.core.commands import EmptyArg     # (see below)
from chimerax.core.commands import Or, Bounded  # Argument modifiers
from chimerax.geometry import Place, Places, translation, rotation
from chimerax.map_data import mrc
from chimerax.map.volume import Volume
import numpy as np
from scipy.spatial import KDTree
from chimerax.open_command import cmd as open_cmd


# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================



def loadtm(session, mip, phi, theta, psi, defocus, template, threshold, pixelsize):
    """Load template matching results."""

    if type(template) == str:
        template = open_cmd.provider_open(session,[template])[0]

    

    mip_data = mrc.open(mip)[0].matrix()
    template.tm_mip_data = mip_data
    phi_data = mrc.open(phi)[0].matrix()
    template.tm_phi_data = phi_data
    theta_data = mrc.open(theta)[0].matrix()
    template.tm_theta_data = theta_data
    psi_data = mrc.open(psi)[0].matrix()
    template.tm_psi_data = psi_data
    defocus_data = mrc.open(defocus)[0].matrix()
    template.tm_defocus_data = defocus_data

    template.tm_pixelsize = pixelsize

    changethreshold(session,template,threshold)



    
    



loadtm_desc = CmdDesc(required=[("mip", FileNameArg),
                              ("phi", FileNameArg),
                              ("theta", FileNameArg),
                              ("psi", FileNameArg),
                              ("defocus", FileNameArg),
                              ("template", Or(ModelArg, FileNameArg)),
                              ("threshold", FloatArg),
                              ("pixelsize", FloatArg)],
                    optional=[])



def changethreshold(session,template,threshold):

    if not hasattr(template,"tm_mip_data"):
        session.logger.error("Model is not a Template Matching model")
        return

    # Transform needed to center template around (0,0,0)
    if type(template) == Volume:
        origin_transform = translation(-0.5 * np.array(template.data.size) * template.tm_pixelsize)
    else:
        session.logger.error("Only Volumes support for now")
        return

    pixel_coordinates = (template.tm_mip_data>threshold).nonzero()

    tree = KDTree(np.transpose(pixel_coordinates))

    phi_of_peaks = template.tm_phi_data[pixel_coordinates]
    theta_of_peaks = template.tm_theta_data[pixel_coordinates]
    psi_of_peaks = template.tm_psi_data[pixel_coordinates]
    defocus_of_peaks = template.tm_defocus_data[pixel_coordinates]

    placements = np.transpose(np.concatenate((np.array(pixel_coordinates),
                                             np.array(phi_of_peaks).reshape(1,-1),
                                             np.array(theta_of_peaks).reshape(1,-1),
                                             np.array(psi_of_peaks).reshape(1,-1)),
                                             0))

    


    t = Places([translation(np.array((x[2],x[1],x[0]))*template.tm_pixelsize) # Translation to correct coordinates
                * rotation(np.array((0,0,1)),-x[5]) # Psi
                * rotation(np.array((0,1,0)),-x[4]) # Theta
                * rotation(np.array((0,0,1)),-x[3]) # Phi
                * origin_transform for x in placements])
    template.set_positions(t)


changethreshold_desc = CmdDesc(required=[("template", ModelArg),
                    ("threshold", FloatArg)],
                    optional=[])
