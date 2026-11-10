"""Case files and the command-line interface."""

import json

import pytest

from dewi_lanjar import builtin_names, case_from_dict, load_case, run_case
from dewi_lanjar.cli import main
from dewi_lanjar.explain import TOPICS, explain


def test_builtin_cases_load_and_validate():
    names = builtin_names()
    assert len(names) >= 10
    for n in names:
        c = load_case(n)
        assert c.name == n
        c.grid(), c.params(), c.forcing()
        assert c.description


def test_overrides_and_errors():
    c = load_case("stommel_1948", ["grid.nx=40", "physics.bc=\"no-slip\""])
    assert c.grid().nx == 40 and c.params().bc == "no-slip"
    with pytest.raises(ValueError):
        load_case("stommel_1948", ["grid.nonsense=1"])
    with pytest.raises(ValueError):
        load_case("no_such_case")
    with pytest.raises(ValueError):
        case_from_dict({"case": {"kind": "magic"}})


def test_toml_roundtrip(tmp_path):
    c = load_case("munk_1950")
    p = tmp_path / "m.toml"
    p.write_text(c.to_toml())
    c2 = load_case(p)
    assert c2.spec == c.spec


@pytest.mark.parametrize("kind,extra", [
    ("steady", {}),
    ("transient", {"run": {"t_end": 4.0, "n_frames": 3}}),
    ("sweep", {"sweep": {"parameter": "beta_factor", "start": 0.0, "stop": 1.0, "n": 3}}),
])
def test_run_case_small(tmp_path, kind, extra):
    spec = {"case": {"name": f"tiny_{kind}", "kind": kind},
            "grid": {"nx": 24, "ny": 24}, "physics": {"delta_s": 0.08},
            "floats": {"n": 20}, "animation": {"dpi": 40, "width": 6, "height": 4}}
    spec["run"] = {"n_frames": 3}
    for k, v in extra.items():
        spec.setdefault(k, {}).update(v)
    res = run_case(case_from_dict(spec), tmp_path, verbose=False)
    for key in ("summary", "json", "npz", "toml", "png", "gif"):
        assert res.files[key].exists(), key
    json.loads(res.files["json"].read_text())


def test_cli_commands(tmp_path, capsys):
    assert main(["list"]) == 0
    assert "stommel_1948" in capsys.readouterr().out
    assert main(["show", "munk_1950"]) == 0
    assert main(["steady", "--delta-s", "0.05", "--nx", "32"]) == 0
    out = capsys.readouterr().out
    assert "exact separable" in out
    assert main(["widths", "--delta-m", "0.03"]) == 0
    assert main(["scales", "--L", "5000", "--lat", "30", "--r", "1e-6"]) == 0
    assert main(["explain"]) == 0
    assert main(["explain", "stommel"]) == 0
    assert main(["explain", "nonsense"]) == 2
    f = tmp_path / "mine.toml"
    assert main(["init", str(f), "--template", "f_plane_control"]) == 0
    assert f.exists()
    assert main(["run", str(f), "--out", str(tmp_path), "--no-gif", "--quiet",
                 "--set", "grid.nx=24", "--set", "grid.ny=24", "--set", "floats.n=0"]) == 0
    assert main(["steady", "--delta-s", "0", "--nx", "16"]) == 2      # error path


def test_validate_quick():
    assert main(["validate", "--quick"]) == 0


def test_every_explain_topic_renders():
    for t in TOPICS:
        assert len(explain(t)) > 200
