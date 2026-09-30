"""
Execute the "Updating Code to the New Gait Module" example from the `skdh.gait`
module docstring on the repository test data, so the documented gait_old -> gait
migration path cannot silently break.
"""

import textwrap
import warnings

import skdh.gait

START = "Updating Code to the New Gait Module"
END = "Adding Custom Gait Endpoints"


def _example_source():
    section = skdh.gait.__doc__.split(START, 1)[1].split(END, 1)[0]
    code = section.split(".. code-block :: python", 1)[1]
    return textwrap.dedent(code)


def test_compatibility_example_runs(gait_input_50):
    t, acc = gait_input_50

    ns = {"t": t, "acc": acc}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        exec(compile(_example_source(), "skdh.gait.__doc__", "exec"), ns)

    res = ns["res"]
    assert ns["bouts"]["Gait Bout Start Index"].size > 0
    assert res["IC Timestamp"].size > 0
    assert "gait speed m1" in res
