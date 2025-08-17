import pytest
import warnings
from pathlib import Path
from nufebmgr.NufebScriptUtils import NufebScriptUtils


def test_raise_neither_path_nor_conents():
    nsu = NufebScriptUtils()

    with pytest.raises(ValueError) as excinfo:
        nsu._strip_comments()
    assert f'Must specify either a file_path or contents' in str(excinfo.value)

def test_raise_both_path_nor_conents():
    nsu = NufebScriptUtils()

    file_path = Path("./tests/data/test1.inputscript")
    file_contents = "foo bar baz"
    with pytest.raises(ValueError) as excinfo:
        nsu._strip_comments(file_path, file_contents)
    assert f'Must specify only a file_path or contents, not both' in str(excinfo.value)

def test_raise_path_type_error():
    nsu = NufebScriptUtils()

    "./tests/data/test1.inputscript"
    with pytest.raises(TypeError) as excinfo:
        nsu._strip_comments("./tests/data/test1.inputscript")
    assert f"'file_path' must be a Path (e.g. pathlib) object, not: str" in str(excinfo.value)

def test_raise_file_not_found_error():
    nsu = NufebScriptUtils()

    file_path = Path("./tests/data/does_not_exist.inputscript")
    with pytest.raises(FileNotFoundError) as excinfo:
        nsu._strip_comments(file_path)

def test_warns_on_suspected_multiline():
    nsu = NufebScriptUtils()

    file_path = Path("./tests/data/script_parsing/simple_oneline_comments.inputscript")
    with pytest.warns(UserWarning) as record:
        nsu._strip_comments(file_path)
    assert any("Warning" in str(w.message) for w in record)

def test_strip_simple_oneline_comments():
    nsu = NufebScriptUtils()

    file_path = Path("./tests/data/script_parsing/simple_oneline_comments.inputscript")
    expected = Path("./tests/data/script_parsing/simple_oneline_comments.inputscript.stripped").read_text()
    stripped = nsu._strip_comments(file_path)
    assert stripped == expected

def test_strip_simple_with_inline_comments():
    nsu = NufebScriptUtils()

    file_path = Path("./tests/data/script_parsing/simple_with_inline_comments.inputscript")
    expected = Path("./tests/data/script_parsing/simple_with_inline_comments.inputscript.stripped").read_text()
    stripped = nsu._strip_comments(file_path)
    assert stripped == expected


# cancelling multi-line
# a comment between single, double, or triple quotes is not a comment:
#     " this # is not a comment"
#     ' this # is not a comment'
#
#     """ this
#     is # not a comment
#     """
#
#     "" this "" # is  a comment
#     ' "" thi '  # is  a comment
#     "" thi "" foo  ""# isnot  a comment""


###multiline
def test_simple_multiline_concat():
    nsu = NufebScriptUtils()

    file_path = Path("./tests/data/script_parsing/simple_with_inline_comments.inputscript")
    expected = Path("./tests/data/script_parsing/simple_with_inline_comments.inputscript.no_multiline").read_text()
    stripped = nsu._concat_multiline(file_path)
    assert stripped == expected

# raise for above
# ends in &.  remove & and line break (detects even with whitespace at end of &. if not in quotes, removes leading whitespace on next line
# acts the same in single and double quotes except only the & is removed, no whitespace contraction
# does not occur in triple quotes