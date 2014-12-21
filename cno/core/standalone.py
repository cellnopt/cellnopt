# -*- python -*-
#
#  This file is part of CNO software
#
#  Copyright (c) 2014 - EBI-EMBL
#
#  File author(s): Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: http://github.com/cellnopt/cellnopt
#
##############################################################################


class Standalone(object):
    """Common class for all standalone applications"""
    def __init__(self, args, user_options):

        # stores the arguments
        self.args = args
        self.user_options = user_options

        if len(args) == 1:
            # shows the help message if no arguments provided
            self.help()
        else:
            # The user values should be used to update the
            # user_options
            options = self.user_options.parse_args(args[1:])
            # Should update the CNOConfig file with the provided options
            for key in self.user_options.config.keys():
                for option in self.user_options.config[key]._get_names():
                    value = getattr(options, option)
                    setattr(getattr( getattr(self.user_options.config, key), option ), 'value', value)
            self.options = options

    def help(self):
        self.user_options.parse_args(["prog", "--help"])

    def report(self):
        """Create report and shows report (or not)"""
        if self.options.onweb is True:
            self.trainer.report(show=True)
        elif self.options.report is True:
            self.trainer.report(show=False)
        else:
            from easydev.console import red
            print(red("No report requested; nothing will be saved or shown"))
            print("use --on-web or --report options")
