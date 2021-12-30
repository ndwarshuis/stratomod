import unittest as ut
import subprocess as sp


class TestPipeline(ut.TestCase):
    # super lame test, just test of the pipeline runs and doesn't blow up
    def test_pipeline_success(self):
        args = [
            "snakemake",
            "--printshellcmds",
            "--verbose",
            "--reason",
            "-j 1",
            "--use-conda",
        ]
        p = sp.run(args)

        self.assertEqual(0, p.returncode)


if __name__ == "__main__":
    ut.main()
