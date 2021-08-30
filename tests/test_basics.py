import pytest
from subprocess import check_call
import sys

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf

@pytest.mark.parametrize("command", ["-version", "-this-is-not-a-valid-parameter"])
def test_basis(command):
    if command == "-this-is-not-a-valid-parameter":
        with pytest.raises(Exception):
            check_call([tcf.WSCLEAN, command])
    else:
        check_call([tcf.WSCLEAN, command])
