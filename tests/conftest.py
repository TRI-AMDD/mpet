import pytest
import os.path as osp
from os import walk


def pytest_addoption(parser):
    parser.addoption(
        "--modDir",
        action="store",
        default=None,
        help="The directory with the modified test results"
    )
    parser.addoption(
        "--baseDir",
        action="store",
        default=None,
        help="The directory with the base test results"
    )
    parser.addoption(
        "--tolerance",
        action="store",
        default=1e-4,
        help="tolerance for test cases"
    )
    parser.addoption(
        "--slowdown_tolerance",
        action="store",
        default=2,
        help="tolerance for slowdown of test cases"
    )
    parser.addoption(
        "--tests",
        action="store",
        default=[],
        nargs='*',
        help="test cases to check"
    )
    parser.addoption(
        "--skip-analytic", action="store_true", default=False, help="skip analytic tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "analytic: mark test as analytic")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--skip-analytic"):
        skip_analytic = pytest.mark.skip(reason="--skip-analytic defined")
        for item in items:
            if "analytic" in item.keywords:
                item.add_marker(skip_analytic)


def pytest_generate_tests(metafunc):
    if "Dirs" in metafunc.fixturenames:
        dir_t = metafunc.config.getoption("modDir")
        dir_b = metafunc.config.getoption("baseDir")
        if metafunc.config.getoption("tests") == []:
            _, directories, _ = next(walk(dir_b))
            directories.sort()
            try:
                directories.remove('plots')  # Directories to ignore
            except Exception:
                pass

            metafunc.parametrize("Dirs", [(osp.join(dir_b, dir),
                                           osp.join(dir_t, dir)) for dir in directories])
        else:
            # Sometimes the tests option is a list, other times it is
            # an array with one element (a string of tests).
            # This handles both.
            tests_str = " ".join(metafunc.config.getoption("tests"))
            tests_lst = tests_str.split()
            metafunc.parametrize("Dirs", [(dir_b + "/" + test, dir_t + "/" + test)
                                          for test in tests_lst])
    if "tol" in metafunc.fixturenames:
        metafunc.parametrize("tol", [float(metafunc.config.getoption("tolerance"))])
    if "testDir" in metafunc.fixturenames:
        dir_t = metafunc.config.getoption("modDir")
        metafunc.parametrize("testDir", [dir_t])
    pytest.slowdown_tolerance = float(metafunc.config.getoption("slowdown_tolerance"))
