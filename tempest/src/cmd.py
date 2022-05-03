# vim: set expandtab shiftwidth=4 softtabstop=4:

from numpy.core.fromnumeric import sort
from chimerax.core.commands import CmdDesc      # Command description
from chimerax.atomic import AtomsArg            # Collection of atoms argument
from chimerax.core.commands import BoolArg      # Boolean argument
from chimerax.core.commands import FileNameArg
from chimerax.core.commands import FloatArg
from chimerax.core.commands import ModelArg, ModelsArg
from chimerax.core.commands import StringArg     # Color argument
from chimerax.core.commands import ColorArg
from chimerax.core.commands import IntArg       # Integer argument
from chimerax.core.commands import EmptyArg     # (see below)
from chimerax.core.commands import Or, Bounded  # Argument modifiers
from chimerax.geometry import Place, Places, translation, rotation
from chimerax.map_data import mrc
from chimerax.map.volume import Volume
import numpy as np
from scipy.spatial import KDTree, distance_matrix
from chimerax.open_command import cmd as open_cmd
from chimerax.core.colors import BuiltinColormaps, BuiltinColors
from chimerax.core.commands import run
from chimerax.std_commands.view import view
from chimerax.core.objects import Objects
# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================


def loadtm(session, mip, unscaled_mip, phi, theta, psi, defocus, template, threshold, pixelsize, image=None):
    """Load template matching results."""
    from chimerax.map.volumecommand import volume
    from chimerax.std_commands.camera import camera
    if type(image) == str:
        image = open_cmd.provider_open(session, [image])[0]
    if type(template) == str:
        template = open_cmd.provider_open(session, [template])[0]
   
    #from chimerax.core.commands import run
    # run(session, f"volume #{image.id[0]} color black color white color white")

    volume(session, volumes=[image], origin=(0.0, 0.0, -2000))
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
    run(session, "view orient")
    im_mean = image.data.full_matrix().mean()
    im_std = np.std(image.data.full_matrix())
    im_min = np.min(image.data.full_matrix())
    im_max = np.max(image.data.full_matrix())
    volume(session, [image], level=[
        (im_min, 1.0),
        (im_mean-im_std, 1.0),
        (im_mean+im_std, 1.0),
        (im_max, 0.999)],
        color=[
        BuiltinColors["black"],
        BuiltinColors["black"],

        BuiltinColors["white"],
        BuiltinColors["white"],
    ])
    camera(session, type='ortho')
    session.image_obj = image
    session.template_obj = template

    


loadtm_desc = CmdDesc(required=[("mip", FileNameArg),
                                ("phi", FileNameArg),
                                ("theta", FileNameArg),
                                ("psi", FileNameArg),
                                ("defocus", FileNameArg),
                                ("template", Or(ModelArg, FileNameArg)),
                                ("threshold", FloatArg),
                                ("pixelsize", FloatArg)],
                      optional=[
    ("image", FileNameArg)
]
)


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
        pixel_coordinates = (
            template.tm_unscaled_mip_data > threshold).nonzero()
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
                                             np.array(psi_of_peaks).reshape(
                                                 1, -1),
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
    template_obj.orig_position = template_obj.surfaces[0].get_positions()

    for surface in template.surfaces:
        surface.set_positions(t)


changethreshold_desc = CmdDesc(required=[("template", ModelArg),
                                         ("threshold", FloatArg)],
                               optional=[("unscaled_mip", BoolArg)])


def loadtm_star(session, filename):
    # Load tempalte matches from a starfile
    import starfile
    from chimerax.map.volumecommand import volume
    from chimerax.std_commands.camera import camera
    matches = starfile.read(filename)
    images = matches["image_filename"].unique()
    templates = matches["template_filename"].unique()
    for image in images:
        for template in templates:
            image_obj = open_cmd.provider_open(session, [image])[0]

            template_obj = open_cmd.provider_open(session, [template])[0]
            if type(template_obj) == Volume:
                origin_transform = translation(-0.5 *
                                               np.array(template_obj.data.size) * 1.5)
            my_matches = matches[(matches["image_filename"] == image) & (
                matches["template_filename"] == template)]
            pixel_size = my_matches["pixel_size"].loc[0]
            if pixel_size != 0.0:
                volume(session, volumes=[image_obj], voxel_size=(
                    pixel_size, pixel_size, pixel_size))
            volume(session, volumes=[image_obj], origin=(0.0, 0.0, -2000))
            im_mean = image_obj.data.full_matrix().mean()
            im_std = np.std(image_obj.data.full_matrix())
            im_min = np.min(image_obj.data.full_matrix())
            im_max = np.max(image_obj.data.full_matrix())
            volume(session, [image_obj], change="image",level=[
                (im_min, 1.0),
                (im_mean-im_std, 1.0),
                (im_mean+im_std, 1.0),
                (im_max, 0.999)],
                color=[
                BuiltinColors["black"],
                BuiltinColors["black"],

                BuiltinColors["white"],
                BuiltinColors["white"],
            ])
            camera(session, type='ortho')
            session.image_obj = image_obj
            session.template_obj = template_obj
            placements = np.transpose(np.array([np.array(my_matches["defocus"]),
                                                np.array(my_matches["y"]),
                                                np.array(my_matches["x"]),
                                                np.array(my_matches["phi"]),
                                                np.array(my_matches["theta"]),
                                                np.array(my_matches["psi"]),
                                                np.array(
                                                    my_matches["peak_value"]),
                                                ]))
            if "display" in my_matches:
                placements = placements[my_matches["display"]]
            # image_obj.data.set_step((pixel_size,pixel_size,pixel_size))
            template_obj.tm_placements = placements
            t = Places([translation(np.array((x[2], x[1], x[0])))  # Translation to correct coordinates
                        * rotation(np.array((0, 0, 1)), -x[5])  # Psi
                        * rotation(np.array((0, 1, 0)), -x[4])  # Theta
                        * rotation(np.array((0, 0, 1)), -x[3])  # Phi
                        * origin_transform for x in placements])
            template_obj.tm_positions = t
            template_obj.orig_position = template_obj.surfaces[0].get_positions()
            for surface in template_obj.surfaces:
                surface.set_positions(t)
    run(session, "view orient")


loadtm_star_desc = CmdDesc(required=[("filename", FileNameArg)
                                     ],
                           optional=[])


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
                f"SELECT SCALED_MIP_OUTPUT_FILE, PHI_OUTPUT_FILE, THETA_OUTPUT_FILE, PSI_OUTPUT_FILE, DEFOCUS_OUTPUT_FILE, USED_PIXEL_SIZE, USED_THRESHOLD, MIP_OUTPUT_FILE, REFERENCE_VOLUME_ASSET_ID, IMAGE_ASSET_ID FROM TEMPLATE_MATCH_LIST WHERE IMAGE_ASSET_ID={image_asset} AND REFERENCE_VOLUME_ASSET_ID={volume_asset}")
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
            f"SELECT SCALED_MIP_OUTPUT_FILE, PHI_OUTPUT_FILE, THETA_OUTPUT_FILE, PSI_OUTPUT_FILE, DEFOCUS_OUTPUT_FILE, USED_PIXEL_SIZE, USED_THRESHOLD, MIP_OUTPUT_FILE, REFERENCE_VOLUME_ASSET_ID, IMAGE_ASSET_ID  FROM TEMPLATE_MATCH_LIST WHERE TEMPLATE_MATCH_ID={tm_index}")
        resultstm = cur.fetchall()
        volume_asset = resultstm[0][8]

    cur.execute(
        f"SELECT VOLUME_ASSET_ID,FILENAME FROM VOLUME_ASSETS WHERE VOLUME_ASSET_ID={volume_asset}")
    vol_results = cur.fetchall()
    cur.execute(
        f"SELECT IMAGE_ASSET_ID,FILENAME FROM IMAGE_ASSETS WHERE IMAGE_ASSET_ID={resultstm[0][9]}")
    img_results = cur.fetchall()
    # session.logger.info(f"{resultstm[0][0]},{resultstm[0][1]},{resultstm[0][2]},{resultstm[0][3]},{resultstm[0][4]},{vol_results[0][1]},{resultstm[0][6]},{resultstm[0][5]}")
    loadtm(session, resultstm[0][0], resultstm[0][7], resultstm[0][1], resultstm[0][2],
           resultstm[0][3], resultstm[0][4], vol_results[0][1], resultstm[0][6], resultstm[0][5], img_results[0][1])

    con.close()



color_by_distance_desc = CmdDesc(required=[("model_to_color", ModelArg),
                                           ("distance_threshold", FloatArg),
                                           ("model_distance_from", ModelArg), 
                                           
                                           ],
                                 optional=[("color_far", ColorArg),("color_close", ColorArg)])

def color_by_distance(session, model_to_color, model_distance_from, distance_threshold, color_far=None, color_close=None):
    
    if not hasattr(model_to_color, "tm_placements"):
        session.logger.error("Model does not have tm_placements")
        return
    
    session.markers = model_distance_from
    points = model_distance_from.atoms.scene_coords
    colors = model_to_color.surfaces[0].get_colors()
    kd = KDTree([[x[2],x[1],x[0]] for x in model_to_color.tm_placements])
    near = kd.query_ball_point(points, distance_threshold)
    if color_close is not None:
        for res in near:
            colors[res] = color_close.uint8x4()
    if color_far is not None:
        for res in near:
            colors[np.setdiff1d(np.arange(colors.shape[0]),res)] = color_far.uint8x4()
    for surface in model_to_color.surfaces:
        surface.set_colors(colors)


color_by_score_desc = CmdDesc(required=[("template", ModelArg),
                                        ("colormap", StringArg)
                                        ],
                              optional=[])


def color_by_score(session, template, colormap):
    cm1 = BuiltinColormaps[colormap]
    mips = np.array([x[6] for x in template.tm_placements], dtype=np.float32)
    min = 7.0
    max = 14.0
    cmn = cm1.rescale_range(min, max)
    colors = cmn.interpolated_rgba8(mips)

    for surface in template.surfaces:
        surface.set_colors(colors)
    run(session,
        f"key {colormap} {' '.join([':'+str(x) for x in np.linspace(min,max+1,num=len(cm1.colors))])}")


loadtm_project_desc = CmdDesc(required=[("cistem_database", FileNameArg)],


                              optional=[("job_number", IntArg),
                                        ("tm_index", IntArg),
                                        ("image_asset", IntArg),
                                        ("volume_asset", IntArg)])


def transfer_instancing(session, template, target):
    
    if hasattr(template, "surfaces") and isinstance(template.surfaces, list):
        positions = template.surfaces[0].get_positions()
    else:
        positions = template.get_positions()

   
    # If the target has already multiple positions, revert to orig_position, otherwise use current positions
    model = target
    if len(model.get_positions()) > 1:
        model.set_positions(model.orig_position)
    orig_position = model.get_positions()
    print(orig_position.array())
    if not hasattr(model, "orig_position"):
        model.orig_position = orig_position
    model.tm_positions = positions
    model.set_positions(positions * orig_position)

transfer_instancing_desc = CmdDesc(required=[("template", ModelArg)],keyword=[("target", ModelArg)],required_arguments=['target'])

def toggle_instancing(session, template):
    if not hasattr(template, "tm_positions"):
        session.logger.error("Model appears not to be loaded by tempest")
        return
    
    if hasattr(template, "surfaces") and isinstance(template.surfaces, list):
        current_positions = template.surfaces[0].get_positions()
        if len(current_positions) > 1:
            for surface in template.surfaces:
                surface.set_positions(template.orig_position)
            run(session, f"view #{template.id[0]}")
        else:
            for surface in template.surfaces:
                surface.set_positions(template.tm_positions)
            run(session, "view orient")
    else:
        current_positions = template.get_positions()
        if len(current_positions) > 1:
            template.set_positions(template.orig_position)
            run(session, f"view #{template.id[0]}")
        else:
            template.set_positions(template.tm_positions)
            run(session, "view orient")

toggle_instancing_desc = CmdDesc(required=[("template", ModelArg)])