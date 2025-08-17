import os
import re
import warnings
from pathlib import Path

class NufebScriptUtils:
    """
    Collection of utility functions for dealing with Nufeb Scripts

    Most of the time scrips will be generated via NufebProject calls, but this is useful for post-hoc management.
    Generally piecemeal functionality as needs arise. For example, this was started off as a simple way to strip
    comments and undo line-continuations to create less brittle tests for the example scripts.

    For now, the design choice is to not have this maintain state (e.g. no init with a specific script) but rather
    follow a no-side effects functional approach. E.g. ``strip_comments(file)`` will return a a copy of the file with
    comments stripped.

    Parsing rules are based on the guidance at: https://docs.lammps.org/Commands_parse.html
    """
    def reduce_to_commands(self, file_path: Path  =None, contents: str = None) -> str:
        """
        Removes multilines, comments, and lines which are just whitespace resulting in a file where, unless triple quotes
         are in play, each line is a single command
        command
        :param file_path:
        :param contents:
        :return:
        """

    def _strip_comments(self, file_path: Path  =None, contents: str = None) -> str:
        """
        Remove comments from a NUFEB inputscript

        Exactly one of `file_path` or `contents` must be provided.

        IMPORTANT For multiline logic to work, first called ``_concat_multiline``

        Comments are defined as per LAMMPS parsing rules. Briefly:
        * any line starting with a # is a comment
        * any multiline starting with a # is completely commented
        * # within a line begins a comment UNLESS it is enclosed in a single/double/triple quote section
           * We have tried to consider pathological cases of weirdly nested quote blocks,
             but cannot guarantee correctness
           * You probably should avoid complex nested quotes anyway.

        Args:
            file_path (Path, optional): Path to the file to read. Defaults to None.
            contents (str, optional): File contents as a string. Defaults to None.

        Returns:
            str: The transformed contents with comments removed.

        Raises:
            ValueError: If neither or both of `file_path` and `contents` are provided.
            TypeError: If `file_path` is not a Path object.
            FileNotFoundError: If `file_path` is given and the file does not exist.
        """

        if (file_path is None) and (contents is None):
            raise ValueError('Must specify either a file_path or contents')

        if (file_path is not None) and (contents is not None):
            raise ValueError('Must specify only a file_path or contents, not both')

        if not isinstance(file_path, Path):
            raise TypeError(f"''file_path' must be a Path (e.g. pathlib) object, not: {type(file_path).__name__}")

        lines = []
        with file_path.open() as f:
            for line in f.readlines():
                if not re.match(r"^\s*#", line):
                    lines.append(line.split("#", 1)[0].rstrip('\r\n'))
                if line.strip() != ''  and line.strip()[-1] == "&":
                    warnings.warn("Warning: _strip_comments encountered a probable line continuation. Run _concat_multiline first", UserWarning)

        stripped = os.linesep.join(lines)
        return stripped+os.linesep

    def _concat_multiline(self, file_path: Path  =None, contents: str = None) -> str:
        """
        Remove comments from a NUFEB inputscript

        Exactly one of `file_path` or `contents` must be provided.

        IMPORTANT For multiline logic to work, first called ``_concat_multiline``

        Comments are defined as per LAMMPS parsing rules. Briefly:
        * any line starting with a # is a comment
        * any multiline starting with a # is completely commented
        * # within a line begins a comment UNLESS it is enclosed in a single/double/triple quote section
           * We have tried to consider pathological cases of weirdly nested quote blocks,
             but cannot guarantee correctness
           * You probably should avoid complex nested quotes anyway.

        Args:
            file_path (Path, optional): Path to the file to read. Defaults to None.
            contents (str, optional): File contents as a string. Defaults to None.

        Returns:
            str: The transformed contents with comments removed.

        Raises:
            ValueError: If neither or both of `file_path` and `contents` are provided.
            TypeError: If `file_path` is not a Path object.
            FileNotFoundError: If `file_path` is given and the file does not exist.
        """
        lines = []
        with file_path.open() as f:
            linebuffer = ''
            for line in f.readlines():
                if line.strip() != '' and line.strip()[-1] == "&":
                    if linebuffer == '':
                        linebuffer = line.rstrip()[0:-1].rstrip()
                    else:
                        linebuffer = linebuffer + ' ' + line.strip()[0:-1].strip()

                else:
                    if linebuffer != '':
                        linebuffer = linebuffer + ' ' + line.strip()
                        lines.append(linebuffer)
                        linebuffer=''
                    else:
                        lines.append(line.rstrip('\r\n'))

        stripped = os.linesep.join(lines)
        return stripped + os.linesep