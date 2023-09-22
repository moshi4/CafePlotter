from pathlib import Path

import pytest


@pytest.fixture
def cafe_base_result_dir() -> Path:
    """CAFE5 Base result directory"""
    testdata_dir = Path(__file__).parent / "testdata"
    return testdata_dir / "singlelambda"


@pytest.fixture
def cafe_gamma_result_dir() -> Path:
    """CAFE5 Gamma result directory"""
    testdata_dir = Path(__file__).parent / "testdata"
    return testdata_dir / "gamma_dist"
