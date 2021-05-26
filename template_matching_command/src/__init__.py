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
        if ci.name == "tm load_manual":
            func = cmd.loadtm
            desc = cmd.loadtm_desc
        elif ci.name == "tm load_project":
            func = cmd.loadtm_project
            desc = cmd.loadtm_project_desc
        elif ci.name == "tm change_threshold":
            func = cmd.changethreshold
            desc = cmd.changethreshold_desc
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
