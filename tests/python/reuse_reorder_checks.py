import os
import sys
from subprocess import STDOUT, check_output

from astropy.io.fits import FITSDiff
from utils import check_and_remove_files, compare_rms_fits

# Append current directory to system path in order to import testconfig variables
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


class TestReuseReorder:
    def test_reuse_reorder_multi_spw(self):
        # Run WSClean with save-reordered and save the temporary reordered files
        save_cmd = f"""{tcf.WSCLEAN} -no-update-model-required -verbose -save-reordered \
            -size 768 768 -scale 0.05arcsec -pol QU -mgain 0.85 -niter 1000 \
            -auto-threshold 3 -join-polarizations -squared-channel-joining \
            -no-mf-weighting -name reference_output {tcf.WORKING_DIR}/{tcf.JVLA_MS}"""
        print("Running: " + save_cmd)
        check_output(save_cmd.split(), stderr=STDOUT)

        # Run WSClean with reuse-reordered and reuse previously generated files
        reuse_reorder_cmd = f"""{tcf.WSCLEAN} -no-update-model-required -verbose -reuse-reordered \
            -size 768 768 -scale 0.05arcsec -pol QU -mgain 0.85 -niter 1000 \
            -auto-threshold 3 -join-polarizations -squared-channel-joining \
            -no-mf-weighting -name expected_output {tcf.WORKING_DIR}/{tcf.JVLA_MS}"""
        print("Running: " + reuse_reorder_cmd)
        check_output(reuse_reorder_cmd.split(), stderr=STDOUT)

        # Compare the generated reference with outputs generated from reuse-reordered files
        reference_filenames = [
            "reference_output-Q-dirty.fits",
            "reference_output-Q-image.fits",
            "reference_output-Q-model.fits",
            "reference_output-Q-residual.fits",
            "reference_output-U-dirty.fits",
            "reference_output-U-image.fits",
            "reference_output-U-model.fits",
            "reference_output-U-residual.fits",
            "reference_output-psf.fits",
        ]
        expected_filenames = [
            "expected_output-Q-dirty.fits",
            "expected_output-Q-image.fits",
            "expected_output-Q-model.fits",
            "expected_output-Q-residual.fits",
            "expected_output-U-dirty.fits",
            "expected_output-U-image.fits",
            "expected_output-U-model.fits",
            "expected_output-U-residual.fits",
            "expected_output-psf.fits",
        ]

        # Compare the RMS of the generated FITS files.
        for i in range(len(reference_filenames)):
            print(
                "Comparing "
                + reference_filenames[i]
                + " & "
                + expected_filenames[i]
            )
            compare_rms_fits(
                reference_filenames[i], expected_filenames[i], threshold=1.0e-6
            )

        # Test teardown
        print("Test teardown: Deleting temporary files and FITS files.")
        file_path = f"{tcf.WORKING_DIR}/"
        files = [
            f for f in os.listdir(file_path) if os.path.isfile(file_path + f)
        ]

        for file in files:
            file_extension = file.rsplit(".", 1)[1]
            if file_extension == "tmp" or file_extension == "fits":
                os.remove(file_path + file)
