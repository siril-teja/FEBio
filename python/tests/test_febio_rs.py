"""Tests for the febio_rs driver. Run tests need febio4 (FEBIO4=... or a repo build)."""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import febio_rs as fr  # noqa: E402

GEN = os.path.join(HERE, "..", "..", "rust", "bench", "gen_model.py")


def make_model(d, steps=2):
    feb = os.path.join(d, "block.feb")
    subprocess.run([sys.executable, GEN, "--n", "3", "--steps", str(steps), "--solver", "rust", "-o", feb],
                   check=True, capture_output=True)
    return feb


try:
    FEBIO = fr.find_febio()
except FileNotFoundError:
    FEBIO = None


class TestFeb(unittest.TestCase):
    def test_set_and_get(self):
        with tempfile.TemporaryDirectory() as d:
            feb = make_model(d)
            out = os.path.join(d, "edited.feb")
            fr.set_values(feb, out, {"Control/time_steps": 7, "Control/solver@type": "solid"})
            self.assertEqual(fr.get_value(out, "Control/time_steps"), "7")
            self.assertEqual(fr.get_value(out, "Control/solver@type"), "solid")
            with self.assertRaises(KeyError):
                fr.set_values(feb, out, {"Control/no_such_thing": 1})


@unittest.skipIf(FEBIO is None, "febio4 not found")
class TestRun(unittest.TestCase):
    def test_run_summary_plot_determinism(self):
        with tempfile.TemporaryDirectory() as d:
            feb = make_model(d)
            a = fr.run(feb, solver="rust", threads=2, workdir=os.path.join(d, "a"), check=True)
            b = fr.run(feb, solver="rust", threads=1, workdir=os.path.join(d, "b"), check=True)
            s = a.summary()
            self.assertTrue(s["ok"])
            self.assertGreater(s["iterations"], 0)
            self.assertEqual(s["iterations"], b.summary()["iterations"])
            pa, pb = a.plot(last_only=False), b.plot(last_only=False)
            self.assertEqual(len(pa.states), 3)          # initial state + 2 steps
            self.assertEqual(pa.states, pb.states)       # bitwise equal for 1 and 2 threads
            self.assertIn("displacement", [v[0] for v in pa.vars])

    def test_failure_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            bad = os.path.join(d, "bad.feb")
            open(bad, "w").write("<febio_spec version='4.0'><Nonsense/></febio_spec>")
            r = fr.run(bad, workdir=d)
            self.assertFalse(r.ok)
            with self.assertRaises(RuntimeError):
                fr.run(bad, workdir=d, check=True)


if __name__ == "__main__":
    unittest.main()
