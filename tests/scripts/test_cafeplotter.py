import subprocess as sp
from pathlib import Path

from cafeplotter import __version__


def test_cli_base_run(cafe_base_result_dir: Path, tmp_path: Path):
    """Test CLI Base(singlelambda) run"""
    cmd = f"cafeplotter -i {cafe_base_result_dir} -o {tmp_path} "
    result = sp.run(cmd, shell=True)
    assert result.returncode == 0


def test_cli_gamma_run(cafe_gamma_result_dir: Path, tmp_path: Path):
    """Test CLI Base(gamma_dist) run"""
    cmd = f"cafeplotter -i {cafe_gamma_result_dir} -o {tmp_path} --format pdf "
    cmd += "--p_label_size 6 --innode_label_size 6 --ignore_branch_length"
    result = sp.run(cmd, shell=True)
    assert result.returncode == 0


def test_cli_show_help():
    """Test CLI show help"""
    result = sp.run("cafeplotter -h", shell=True)
    assert result.returncode == 0


def test_cli_show_version():
    """Test CLI show version"""
    result = sp.run(
        "cafeplotter -v",
        shell=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.rstrip() == f"v{__version__}"
    assert result.returncode == 0
