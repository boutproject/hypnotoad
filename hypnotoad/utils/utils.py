# Copyright 2019 J.T. Omotani
#
# Contact John Omotani john.omotani@ukaea.uk
#
# This file is part of Hypnotoad 2.
#
# Hypnotoad 2 is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Hypnotoad 2 is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# Hypnotoad 2.  If not, see <http://www.gnu.org/licenses/>.

import sys


def with_default(value, default):

    if value is not None:
        return value

    return default


def list_loaded_modules():
    module_list = []
    for module in sys.modules:
        module_split = module.split(".")
        if not (len(module_split) > 1 and module_split[0] in sys.modules):
            # Only get versions of 'top-level' modules, e.g. from scipy but not
            # scipy.interpolate, etc.
            module_list.append(module)
    module_list.sort()

    return module_list


def module_versions_formatted():
    module_list = list_loaded_modules()

    module_string = "{\n"
    for module in module_list:
        try:
            module_string += f"{module}: {sys.modules[module].__version__},\n"
        except AttributeError:
            module_string += f"{module}: unknown,\n"
    module_string += "}"

    return module_string
