# vim: set expandtab shiftwidth=4 softtabstop=4:

from numpy.core.fromnumeric import sort
from chimerax.core.commands import CmdDesc      # Command description
from chimerax.atomic import AtomsArg            # Collection of atoms argument
from chimerax.core.commands import BoolArg      # Boolean argument
from chimerax.core.commands import FileNameArg
from chimerax.core.commands import FloatArg
from chimerax.core.commands import ModelArg
from chimerax.core.commands import StringArg     # Color argument
from chimerax.core.commands import IntArg       # Integer argument
from chimerax.core.commands import EmptyArg     # (see below)
from chimerax.core.commands import Or, Bounded  # Argument modifiers
from chimerax.geometry import Place, Places, translation, rotation
from chimerax.map_data import mrc
from chimerax.map.volume import Volume
import numpy as np
from scipy.spatial import KDTree, distance_matrix
from chimerax.open_command import cmd as open_cmd
from chimerax.core.colors import BuiltinColormaps
from chimerax.core.commands import run
# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================


def loadtm(session, mip, unscaled_mip, phi, theta, psi, defocus, template, threshold, pixelsize):
    """Load template matching results."""

    if type(template) == str:
        template = open_cmd.provider_open(session, [template])[0]

    mip_data = mrc.open(mip)[0].matrix()

    template.tm_mip_data = mip_data
    unscaled_mip_data = mrc.open(unscaled_mip)[0].matrix()

    template.tm_unscaled_mip_data = unscaled_mip_data
    phi_data = mrc.open(phi)[0].matrix()
    template.tm_phi_data = phi_data
    theta_data = mrc.open(theta)[0].matrix()
    template.tm_theta_data = theta_data
    psi_data = mrc.open(psi)[0].matrix()
    template.tm_psi_data = psi_data
    defocus_data = mrc.open(defocus)[0].matrix()
    template.tm_defocus_data = defocus_data

    template.tm_pixelsize = pixelsize

    changethreshold(session, template, threshold)


loadtm_desc = CmdDesc(required=[("mip", FileNameArg),
                                ("phi", FileNameArg),
                                ("theta", FileNameArg),
                                ("psi", FileNameArg),
                                ("defocus", FileNameArg),
                                ("template", Or(ModelArg, FileNameArg)),
                                ("threshold", FloatArg),
                                ("pixelsize", FloatArg)],
                      optional=[])


def changethreshold(session, template, threshold, unscaled_mip=False):

    if not hasattr(template, "tm_mip_data"):
        session.logger.error("Model is not a Template Matching model")
        return

    # Transform needed to center template around (0,0,0)
    if type(template) == Volume:
        origin_transform = translation(-0.5 *
                                       np.array(template.data.size) * template.tm_pixelsize)
    else:
        session.logger.error("Only Volumes support for now")
        return

    

    if unscaled_mip:
        pixel_coordinates = (template.tm_unscaled_mip_data > threshold).nonzero()
        mips = template.tm_unscaled_mip_data[pixel_coordinates]
    else:
        pixel_coordinates = (template.tm_mip_data > threshold).nonzero()
        mips = template.tm_mip_data[pixel_coordinates]

    sort_index = np.argsort(mips)

    pixel_coordinates = np.array(pixel_coordinates)[:, sort_index]
    pixel_coordinates = (
        pixel_coordinates[0], pixel_coordinates[1], pixel_coordinates[2])

    if unscaled_mip:
        mips = template.tm_unscaled_mip_data[pixel_coordinates]
    else:
        mips = template.tm_mip_data[pixel_coordinates]

    phi_of_peaks = template.tm_phi_data[pixel_coordinates]
    theta_of_peaks = template.tm_theta_data[pixel_coordinates]
    psi_of_peaks = template.tm_psi_data[pixel_coordinates]
    defocus_of_peaks = template.tm_defocus_data[pixel_coordinates]
    pc = (template.tm_defocus_data[pixel_coordinates],
          pixel_coordinates[1], pixel_coordinates[2])
    coors = np.transpose((pixel_coordinates[1], pixel_coordinates[2]))
    dm = distance_matrix(coors, coors)

    display = np.zeros(mips.shape[0])

    for i, peak in enumerate(mips):
        if i == mips.shape[0]-1 or np.min(dm[i][i+1:]) > 10.0:
            display[i] = 1
        if i == mips.shape[0]-1:
            continue

    placements = np.transpose(np.concatenate((np.array(pc),
                                             np.array(phi_of_peaks).reshape(
                                                 1, -1),
                                             np.array(theta_of_peaks).reshape(
                                                 1, -1),
                                             np.array(psi_of_peaks).reshape(1, -1),
                                             np.array(mips).reshape(1, -1)),
                                             0))
    placements = placements[display > 0]
    template.tm_placements = placements
    t = Places([translation(np.array((x[2]*template.tm_pixelsize, x[1]*template.tm_pixelsize, x[0])))  # Translation to correct coordinates
                * rotation(np.array((0, 0, 1)), -x[5])  # Psi
                * rotation(np.array((0, 1, 0)), -x[4])  # Theta
                * rotation(np.array((0, 0, 1)), -x[3])  # Phi
                * origin_transform for x in placements])
    template.tm_positions = t
    for surface in template.surfaces:
        surface.set_positions(t)


changethreshold_desc = CmdDesc(required=[("template", ModelArg),
                                         ("threshold", FloatArg)],
                               optional=[("unscaled_mip", BoolArg)])


def loadtm_project(session, cistem_database, tm_index=None, image_asset=None, volume_asset=None, job_number=None):
    import sqlite3

    con = sqlite3.connect(cistem_database)
    cur = con.cursor()

    if tm_index is None:
        cur.execute(
            f"SELECT IMAGE_ASSET_ID,FILENAME FROM IMAGE_ASSETS WHERE IMAGE_ASSET_ID={image_asset}")
        results = cur.fetchall()

        if len(results) > 0:
            micrograph = open_cmd.provider_open(session, [results[0][1]])[0]

            cur.execute(
                f"SELECT SCALED_MIP_OUTPUT_FILE, PHI_OUTPUT_FILE, THETA_OUTPUT_FILE, PSI_OUTPUT_FILE, DEFOCUS_OUTPUT_FILE, USED_PIXEL_SIZE, USED_THRESHOLD, MIP_OUTPUT_FILE, REFERENCE_VOLUME_ASSET_ID FROM TEMPLATE_MATCH_LIST WHERE IMAGE_ASSET_ID={image_asset} AND REFERENCE_VOLUME_ASSET_ID={volume_asset}")
            resultstm = cur.fetchall()
            if len(resultstm) > 1:
                session.logger.info("Multiple results found, using first")
            elif len(resultstm) == 1:
                session.logger.info("Results found")
            else:
                session.logger.error(
                    "No Template Match results for this volume found!")
                con.close()
                return
        else:
            session.logger.error(
            "No Template Match results for this volume found!")

            con.close()
            return
    else:
        cur.execute(
                f"SELECT SCALED_MIP_OUTPUT_FILE, PHI_OUTPUT_FILE, THETA_OUTPUT_FILE, PSI_OUTPUT_FILE, DEFOCUS_OUTPUT_FILE, USED_PIXEL_SIZE, USED_THRESHOLD, MIP_OUTPUT_FILE, REFERENCE_VOLUME_ASSET_ID FROM TEMPLATE_MATCH_LIST WHERE TEMPLATE_MATCH_ID={tm_index}")
        resultstm = cur.fetchall()
        volume_asset = resultstm[0][8]
    

    cur.execute(
        f"SELECT VOLUME_ASSET_ID,FILENAME FROM VOLUME_ASSETS WHERE VOLUME_ASSET_ID={volume_asset}")
    vol_results = cur.fetchall()
    # session.logger.info(f"{resultstm[0][0]},{resultstm[0][1]},{resultstm[0][2]},{resultstm[0][3]},{resultstm[0][4]},{vol_results[0][1]},{resultstm[0][6]},{resultstm[0][5]}")
    loadtm(session, resultstm[0][0], resultstm[0][7], resultstm[0][1], resultstm[0][2],
            resultstm[0][3], resultstm[0][4], vol_results[0][1], resultstm[0][6], resultstm[0][5])

       
    con.close()


color_by_score_desc = CmdDesc(required=[("template", ModelArg),
                                        ("colormap", StringArg)
                                        ],
                              optional=[])

def color_by_score(session, template, colormap):
    cm1 = BuiltinColormaps[colormap]
    mips = np.array([x[6] for x in template.tm_placements],dtype=np.float32)
    min = np.min(np.floor(mips))
    max = np.max(np.ceil(mips))
    cmn = cm1.rescale_range(min,max)
    colors = cmn.interpolated_rgba8(mips)

    for surface in template.surfaces:   
        surface.set_colors(colors)
    run(session,f"key {colormap} {' '.join([':'+str(x) for x in np.arange(min,max+1)])}")


loadtm_project_desc = CmdDesc(required=[("cistem_database", FileNameArg)],
                                        
                                       
                              optional=[("job_number", IntArg),
                                        ("tm_index", IntArg),
                                        ("image_asset", IntArg),
                                        ("volume_asset", IntArg)])


def add_molecule(session, template, molecule):
    positions = template.get_positions()
    molecule.set_positions(positions)

def show_tm(session, template, show):

    if show:
        template.set_positions(template.tm_positions)
    else:
        template.set_positions()
