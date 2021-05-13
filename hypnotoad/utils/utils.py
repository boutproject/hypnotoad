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
    try:
        # Available in Python>=3.10
        stdlib_modules = sys.stdlib_module_names
    except AttributeError:
        stdlib_modules = [
            "string",
            "re",
            "difflib",
            "textwrap",
            "unicodedata",
            "stringprep",
            "readline",
            "rlcompleter",
            "struct",
            "codecs",
            "datetime",
            "zoneinfo",
            "calendar",
            "collections",
            "heapq",
            "bisect",
            "array",
            "weakref",
            "types",
            "copy",
            "pprint",
            "reprlib",
            "enum",
            "graphlib",
            "numbers",
            "math",
            "cmath",
            "decimal",
            "fractions",
            "random",
            "statistics",
            "itertools",
            "functools",
            "operator",
            "pathlib",
            "os",
            "fileinput",
            "stat",
            "filecmp",
            "tempfile",
            "glob",
            "fnmatch",
            "linecache",
            "shutil",
            "pickle",
            "copyreg",
            "shelve",
            "marshal",
            "dbm",
            "sqlite3",
            "zlib",
            "gzip",
            "bz2",
            "lzma",
            "zipfile",
            "tarfile",
            "csv",
            "configparser",
            "netrc",
            "xdrlib",
            "plistlib",
            "hashlib",
            "hmac",
            "secrets",
            "io",
            "time",
            "argparse",
            "getopt",
            "logging",
            "getpass",
            "curses",
            "platform",
            "errno",
            "ctypes",
            "threading",
            "multiprocessing",
            "concurrent",
            "subprocess",
            "sched",
            "queue",
            "contextvars",
            "_thread",
            "asyncio",
            "socket",
            "ssl",
            "select",
            "selectors",
            "asyncore",
            "asynchat",
            "signal",
            "mmap",
            "email",
            "json",
            "mailcap",
            "mailbox",
            "mimetypes",
            "base64",
            "binhex",
            "binascii",
            "quopri",
            "uu",
            "html",
            "xml",
            "webbrowser",
            "cgi",
            "cgitb",
            "wsgiref",
            "urllib",
            "http",
            "ftplib",
            "poplib",
            "imaplib",
            "nntplib",
            "smtplib",
            "smtpd",
            "telnetlib",
            "uuid",
            "socketserver",
            "xmlrpc",
            "ipaddress",
            "audioop",
            "aifc",
            "sunau",
            "wave",
            "chunk",
            "colorsys",
            "imghdr",
            "sndhdr",
            "ossaudiodev",
            "gettext",
            "locale",
            "turtle",
            "cmd",
            "shlex",
            "tkinter",
            "typing",
            "pydoc",
            "doctest",
            "unittest",
            "test",
            "bdb",
            "faulthandler",
            "pdb",
            "timeit",
            "trace",
            "tracemalloc",
            "distutils",
            "ensurepip",
            "venv",
            "zipapp",
            "sys",
            "sysconfig",
            "builtins",
            "__main__",
            "warnings",
            "dataclasses",
            "contextlib",
            "abc",
            "atexit",
            "traceback",
            "__future__",
            "gc",
            "inspect",
            "site",
            "code",
            "codeop",
            "zipimport",
            "pkgutil",
            "modulefinder",
            "runpy",
            "importlib",
            "parser",
            "ast",
            "symtable",
            "symbol",
            "token",
            "keyword",
            "tokenize",
            "tabnanny",
            "pyclbr",
            "py_compile",
            "compileall",
            "dis",
            "pickletools",
            "formatter",
            "msilib",
            "msvcrt",
            "winreg",
            "winsound",
            "posix",
            "pwd",
            "spwd",
            "grp",
            "crypt",
            "termios",
            "tty",
            "pty",
            "fcntl",
            "pipes",
            "resource",
            "nis",
            "syslog",
            "optparse",
            "imp",
        ]
    system_modules = list(stdlib_modules) + list(sys.builtin_module_names)

    module_list = []
    for module in sys.modules:
        module_split = module.split(".")
        if not (
            (len(module_split) > 1 and module_split[0] in sys.modules)
            or module_split[0] in system_modules
            or module[0] == "_"
        ):
            # Only get versions of 'top-level' modules, e.g. from scipy but not
            # scipy.interpolate, etc.
            # Ignore stdlib and built-in modules
            # Ignore internal modules (beginning with "_")
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
