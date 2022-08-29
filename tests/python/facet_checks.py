import pytest
import shutil
import sys
from utils import (
    assert_taql,
    basic_image_check,
    check_and_remove_files,
    compare_rms_fits,
    validate_call,
)

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def gridders():
    return {"wstacking": "", "wgridder": "-use-wgridder", "idg": "-use-idg"}


def predict_full_image(ms, gridder):
    """Predict full image"""
    s = f"{tcf.WSCLEAN} -predict {gridder} -name point-source {ms}"
    validate_call(s.split())


def predict_facet_image(ms, gridder, apply_facet_beam):
    name = "point-source"
    facet_beam = "-apply-facet-beam -mwa-path ." if apply_facet_beam else ""
    if apply_facet_beam:
        shutil.copyfile(f"{name}-model.fits", f"{name}-model-pb.fits")

    # Predict facet based image
    s = f"{tcf.WSCLEAN} -predict {gridder} {facet_beam} -facet-regions {tcf.FACETFILE_4FACETS} -name {name} {ms}"
    validate_call(s.split())


def deconvolve_facets(ms, gridder, reorder, mpi, apply_beam=False):
    nthreads = 4
    mpi_cmd = f"mpirun -tag-output -np {nthreads} {tcf.WSCLEAN_MP}"
    thread_cmd = f"{tcf.WSCLEAN} -parallel-gridding {nthreads}"
    reorder_ms = "-reorder" if reorder else "-no-reorder"
    s = [
        mpi_cmd if mpi else thread_cmd,
        gridder,
        tcf.DIMS_SMALL,
        reorder_ms,
        "-niter 1000000 -auto-threshold 5 -mgain 0.8",
        f"-facet-regions {tcf.FACETFILE_4FACETS}",
        f"-name facet-imaging{reorder_ms}",
        "-mwa-path . -apply-facet-beam" if apply_beam else "",
        "-v",
        ms,
    ]
    validate_call(" ".join(s).split())


@pytest.mark.usefixtures("prepare_mock_ms", "prepare_model_image")
class TestFacets:
    def test_makepsfonly(self):
        """
        Test that wsclean with the -make-psf-only flag exits gracefully and
        that the psf passes basic checks.
        """
        s = [
            tcf.WSCLEAN,
            "-make-psf-only -name facet-psf-only",
            tcf.DIMS_SMALL,
            f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.MWA_MOCK_MS}",
        ]
        validate_call(" ".join(s).split())

        basic_image_check("facet-psf-only-psf.fits")

    # Test assumes that IDG and EveryBeam are installed
    @pytest.mark.parametrize("gridder", gridders().items())
    def test_stitching(self, gridder):
        """Test stitching of the facets"""
        prefix = f"facet-stitch-{gridder[0]}"
        s = [
            tcf.WSCLEAN,
            "-quiet",
            gridder[1],
            tcf.DIMS_SMALL,
            "" if (gridder[0] == "idg") else "-pol XX,YY",
            f"-facet-regions {tcf.FACETFILE_2FACETS}",
            f"-name {prefix}",
            tcf.MWA_MOCK_MS,
        ]
        validate_call(" ".join(s).split())
        fpaths = (
            [prefix + "-dirty.fits", prefix + "-image.fits"]
            if (gridder[0] == "idg")
            else [
                prefix + "-XX-dirty.fits",
                prefix + "-YY-dirty.fits",
                prefix + "-XX-image.fits",
                prefix + "-YY-image.fits",
            ]
        )
        check_and_remove_files(fpaths, remove=True)

    # FIXME: we should test wstacking and wgridder here too
    # but they fail on the taql assertion
    @pytest.mark.parametrize("gridder", ["-use-wgridder"])
    @pytest.mark.parametrize("apply_facet_beam", [False, True])
    def test_predict(self, gridder, apply_facet_beam):
        """
        Test predict only run

        Parameters
        ----------
        gridder : str
            wsclean compatible description of gridder to be used.
        """

        predict_facet_image(tcf.MWA_MOCK_FACET, gridder, apply_facet_beam)

        # A numerical check can only be performed in case no DD effects were applied.
        if not apply_facet_beam:
            predict_full_image(tcf.MWA_MOCK_FULL, gridder)
            taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tcf.MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
            assert_taql(taql_command)

    @pytest.mark.parametrize("gridder", ["-use-wgridder"])
    @pytest.mark.parametrize("reorder", [False, True])
    @pytest.mark.parametrize("mpi", [False, True])
    def test_facetdeconvolution(self, gridder, reorder, mpi):
        """
        Test facet-based deconvolution

        Parameters
        ----------
        gridder : str
            wsclean compatible description of gridder to be used.
        reorder : bool
            Reorder MS?
        mpi : bool
            True: Use MPI for parallel gridding.
            False: Use multi-threading for parallel gridding.
        """
        # Parametrization causes some overhead in that predict of full image is run for
        # every parametrization
        predict_full_image(tcf.MWA_MOCK_FULL, gridder)

        # Make sure old versions of the facet mock ms are removed
        shutil.rmtree(tcf.MWA_MOCK_FACET)

        # Copy output to new MS, swap DATA column, and remove MODEL_DATA
        validate_call(
            f"cp -r {tcf.MWA_MOCK_FULL} {tcf.MWA_MOCK_FACET}".split()
        )
        assert shutil.which("taql") is not None, "taql executable not found!"

        validate_call(
            [
                "taql",
                "-noph",
                f"UPDATE {tcf.MWA_MOCK_FACET} SET DATA=MODEL_DATA",
            ]
        )
        validate_call(
            [
                "taql",
                "-noph",
                f"ALTER TABLE {tcf.MWA_MOCK_FACET} DROP COLUMN MODEL_DATA",
            ]
        )
        taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tcf.MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.DATA, 4e-3))"
        assert_taql(taql_command)

        deconvolve_facets(tcf.MWA_MOCK_FACET, gridder, reorder, mpi)

        taql_command = f"select from {tcf.MWA_MOCK_FACET} where not all(near(DATA,MODEL_DATA, 4e-3))"
        assert_taql(taql_command)

    @pytest.mark.parametrize("mpi", [False, True])
    def test_facetbeamimages(self, mpi):
        """
        Basic checks of the generated images when using facet beams. For each image,
        test that the pixel values are valid (not NaN/Inf) and check the percentage
        of zero pixels.
        """

        deconvolve_facets(tcf.MWA_MOCK_FACET, "-use-wgridder", True, mpi, True)

        basic_image_check("facet-imaging-reorder-psf.fits")
        basic_image_check("facet-imaging-reorder-dirty.fits")

    def test_parallel_gridding(self):
        # Compare serial, threaded and mpi run for facet based imaging
        # with h5 corrections. Number of used threads/processes is
        # deliberately chosen smaller than the number of facets.
        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        names = ["facets-h5-serial", "facets-h5-threaded", "facets-h5-mpi"]
        wsclean_commands = [
            tcf.WSCLEAN,
            f"{tcf.WSCLEAN} -parallel-gridding 3",
            f"mpirun -np 3 {tcf.WSCLEAN_MP}",
        ]
        for name, command in zip(names, wsclean_commands):
            # -j 1 to ensure deterministic iteration over visibilities
            s = f"{command} -j 1 -use-wgridder -name {name} -apply-facet-solutions mock_soltab_2pol.h5 ampl000,phase000 -pol xx,yy -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -join-polarizations -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {tcf.MWA_MOCK_MS}"
            validate_call(s.split())

        # Compare images, the threshold is chosen relatively large since the difference
        # seems to fluctuate somewhat between runs.
        threshold = 1e-6
        compare_rms_fits(
            f"{names[0]}-YY-image.fits", f"{names[1]}-YY-image.fits", threshold
        )
        compare_rms_fits(
            f"{names[0]}-YY-image.fits", f"{names[2]}-YY-image.fits", threshold
        )

    @pytest.mark.parametrize("beam", [False, True])
    @pytest.mark.parametrize(
        "h5file",
        [
            None,
            ["mock_soltab_2pol.h5"],
            ["mock_soltab_2pol.h5", "mock_soltab_2pol.h5"],
        ],
    )
    def test_multi_ms(self, beam, h5file):
        """
        Check that identical images are obtained in case multiple (identical) MSets and H5Parm
        files are provided compared to imaging one MSet
        """

        h5download = f"wget -N -q {tcf.WSCLEAN_DATA_URL}/mock_soltab_2pol.h5"
        validate_call(h5download.split())

        # Make a new copy of tcf.MWA_MOCK_MS into two MSets
        validate_call(f"cp -r {tcf.MWA_MOCK_MS} {tcf.MWA_MOCK_COPY_1}".split())
        validate_call(f"cp -r {tcf.MWA_MOCK_MS} {tcf.MWA_MOCK_COPY_2}".split())

        names = ["facets-single-ms", "facets-multiple-ms"]
        commands = [
            f"{tcf.MWA_MOCK_MS}",
            f"{tcf.MWA_MOCK_COPY_1} {tcf.MWA_MOCK_COPY_2}",
        ]

        if beam:
            commands = [
                "-mwa-path . -apply-facet-beam " + command
                for command in commands
            ]

        if h5file is not None:
            commands[0] = (
                f"-apply-facet-solutions {h5file[0]} ampl000,phase000 "
                + commands[0]
            )
            commands[1] = (
                f"-apply-facet-solutions {','.join(h5file)} ampl000,phase000 "
                + commands[1]
            )

        # Note: -j 1 enabled to ensure deterministic iteration over visibilities
        for name, command in zip(names, commands):
            s = f"{tcf.WSCLEAN} -j 1 -nmiter 2 -use-wgridder -name {name} -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {command}"
            validate_call(s.split())

        # Compare images.
        threshold = 1.0e-6
        compare_rms_fits(
            f"{names[0]}-image.fits", f"{names[1]}-image.fits", threshold
        )

        # Model data columns should be equal
        taql_commands = [
            f"select from {tcf.MWA_MOCK_MS} t1, {tcf.MWA_MOCK_COPY_1} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,1e-6))"
        ]
        taql_commands.append(
            f"select from {tcf.MWA_MOCK_COPY_1} t1, {tcf.MWA_MOCK_COPY_2} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,1e-6))"
        )
        # assert_taql(taql_command for taql_command in taql_commands)
        for taql_command in taql_commands:
            assert_taql(taql_command)
