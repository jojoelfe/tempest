# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.toolshed import BundleAPI


# Subclass from chimerax.core.toolshed.BundleAPI and
# override the method for registering commands,
# inheriting all other methods from the base class.
class _MyAPI(BundleAPI):

    api_version = 1     # register_command called with BundleInfo and
                        # CommandInfo instance instead of command name
                        # (when api_version==0)

    # Override method
    @staticmethod
    def start_tool(session, bi, ti):
        # session is an instance of chimerax.core.session.Session
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ti is an instance of chimerax.core.toolshed.ToolInfo

        # This method is called once for each time the tool is invoked.

        # We check the name of the tool, which should match one of the
        # ones listed in bundle_info.xml (without the leading and
        # trailing whitespace), and create and return an instance of the
        # appropriate class from the ``gui`` module.
        from . import tm_gui
        if ti.name == "Tempest":
            return tm_gui.TemplateMatchingTool(session, ti.name)
        raise ValueError("trying to start unknown tool: %s" % ti.name)

    @staticmethod
    def register_command(bi, ci, logger):
        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ci is an instance of chimerax.core.toolshed.CommandInfo
        # logger is an instance of chimerax.core.logger.Logger

        # This method is called once for each command listed
        # in bundle_info.xml.  Since we list two commands,
        # we expect two calls to this method.

        # We check the name of the command, which should match
        # one of the ones listed in bundle_info.xml
        # (without the leading and trailing whitespace),
        # and import the function to call and its argument
        # description from the ``cmd`` module.
        # If the description does not contain a synopsis, we
        # add the one in ``ci``, which comes from bundle_info.xml.
        from . import cmd
        if ci.name == "tempest load_manual":
            func = cmd.loadtm
            desc = cmd.loadtm_desc
        elif ci.name == "tempest load_project":
            func = cmd.loadtm_project
            desc = cmd.loadtm_project_desc
        elif ci.name == "tempest load_star":
            func = cmd.loadtm_star
            desc = cmd.loadtm_star_desc
        elif ci.name == "tempest change_threshold":
            func = cmd.changethreshold
            desc = cmd.changethreshold_desc
        elif ci.name == "tempest color_by_score":
            func = cmd.color_by_score
            desc = cmd.color_by_score_desc
        elif ci.name == "tempest color_by_distance":
            func = cmd.color_by_distance
            desc = cmd.color_by_distance_desc
        elif ci.name == "tempest filter_by_distance":
            func = cmd.filter_by_distance
            desc = cmd.filter_by_distance_desc
        elif ci.name == "tempest toggle_instancing":
            func = cmd.toggle_instancing
            desc = cmd.toggle_instancing_desc
        elif ci.name == "tempest transfer_instancing":
            func = cmd.transfer_instancing
            desc = cmd.transfer_instancing_desc
        else:
            raise ValueError("trying to register unknown command: %s" % ci.name)
        if desc.synopsis is None:
            desc.synopsis = ci.synopsis

        # We then register the function as the command callback
        # with the chimerax.core.commands module.
        from chimerax.core.commands import register
        register(ci.name, desc, func)


# Create the ``bundle_api`` object that ChimeraX expects.
bundle_api = _MyAPI()
