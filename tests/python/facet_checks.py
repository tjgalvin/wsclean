import shutil
import sys

import casacore.tables
import h5py
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from utils import (
    assert_taql,
    basic_image_check,
    check_and_remove_files,
    compare_rms_fits,
    compute_rms,
    validate_call,
)

# Append current directory to system path in order to import testconfig
sys.path.append(".")

# Import configuration variables as test configuration (tcf)
import config_vars as tcf


def gridders():
    return {
        "wstacking": "-gridder wstacking",
        "wgridder": "-gridder wgridder",
        "idg": "-gridder idg",
    }


def predict_full_image(ms, gridder):
    """Predict full image"""
    s = f"{tcf.WSCLEAN} -predict -gridder {gridder} -name point-source {ms}"
    validate_call(s.split())


def predict_facet_image(
    ms, gridder="wgridder", apply_beam=False, wsclean_command=tcf.WSCLEAN
):
    name = "point-source"
    facet_beam = "-apply-facet-beam -mwa-path ." if apply_beam else ""
    if apply_beam:
        shutil.copyfile(f"{name}-model.fits", f"{name}-model-fpb.fits")

    # Predict facet based image
    s = (
        f"{wsclean_command} -predict -gridder {gridder} {facet_beam} "
        f"-facet-regions {tcf.FACETFILE_4FACETS} -name {name} {ms}"
    )
    validate_call(s.split())


def deconvolve_facets(ms, gridder, reorder, mpi, apply_beam=False):
    nthreads = 4
    mpi_cmd = f"mpirun -tag-output -np {nthreads} {tcf.WSCLEAN_MP}"
    thread_cmd = f"{tcf.WSCLEAN} -parallel-gridding {nthreads}"
    reorder_ms = "-reorder" if reorder else "-no-reorder"
    facet_beam = "-mwa-path . -apply-facet-beam" if apply_beam else ""
    s = (
        f"{mpi_cmd if mpi else thread_cmd} -gridder {gridder} {reorder_ms} "
        f"{tcf.DIMS_SMALL} -niter 1000000 -auto-threshold 5 -mgain 0.8 "
        f"-facet-regions {tcf.FACETFILE_4FACETS} {facet_beam} "
        f"-name facet-imaging{reorder_ms} -v {ms}"
    )
    validate_call(s.split())


def create_pointsource_grid_skymodel(
    skymodel_filename, grid_size, nr_pixels, wcs
):
    """
    Writes a skymodel file for a square grid of point sources in a square image.

    Parameters
    ----------
    skymodel_filename: str
    grid_size: int
        Number of point sources (in one direction)
    nr_pixels: int
        Number of pixels in the image (in one direction)
    wcs: astropy.wcs.WCS
        World coordinate system of the image

    Returns
    -------
    list of tuples
        A list of source/pixel positions of length grid_size*grid_size
    """
    source_positions = []
    source_pixel_index_range = (np.arange(grid_size)) * (
        nr_pixels // grid_size
    ) + (nr_pixels // grid_size // 2)
    with open(skymodel_filename, "w") as sky_model_file:
        print(
            "Format = Name, Patch, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, ReferenceFrequency='150000000', MajorAxis, MinorAxis, Orientation",
            file=sky_model_file,
        )
        for i, idx0 in enumerate(source_pixel_index_range):
            for j, idx1 in enumerate(source_pixel_index_range):
                sky = wcs.pixel_to_world(idx0, idx1, 0, 0)
                print(
                    f",direction_{i}{j},,{sky[0].ra.rad},{sky[0].dec.rad},,,,,,,",
                    file=sky_model_file,
                )
                print(
                    f"source-{i}-{j},direction_{i}{j},POINT,{sky[0].ra.rad},{sky[0].dec.rad},1.0,[],false,150000000,,,",
                    file=sky_model_file,
                )
                source_positions.append((idx0, idx1))

    return source_positions


@pytest.mark.usefixtures(
    "prepare_mock_ms", "prepare_model_image", "prepare_mock_soltab"
)
class TestFacets:
    def test_makepsfonly(self):
        """
        Test that wsclean with the -make-psf-only flag exits gracefully and
        that the psf passes basic checks.
        """
        s = (
            f"{tcf.WSCLEAN} -name facet-psf-only -make-psf-only "
            f"-facet-regions {tcf.FACETFILE_4FACETS} "
            f"{tcf.DIMS_SMALL} {tcf.MWA_MOCK_MS}"
        )
        validate_call(s.split())

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

    # FIXME: we should test wstacking here too
    # but it fails on the taql assertion
    @pytest.mark.parametrize("gridder", ["wgridder"])
    @pytest.mark.parametrize("apply_facet_beam", [False, True])
    def test_predict(self, gridder, apply_facet_beam, tmp_mwa_mock_facet):
        """
        Test predict only run

        Parameters
        ----------
        gridder : str
            wsclean compatible description of gridder to be used.
        """

        predict_facet_image(tmp_mwa_mock_facet, gridder, apply_facet_beam)

        # A numerical check can only be performed in case no DD effects were applied.
        if not apply_facet_beam:
            predict_full_image(tcf.MWA_MOCK_FULL, gridder)
            taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tmp_mwa_mock_facet} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
            assert_taql(taql_command)

    @pytest.mark.parametrize("gridder", ["wgridder"])
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

    def test_read_only_ms(self):
        chmod = f"chmod a-w -R {tcf.MWA_MOCK_FULL}"
        validate_call(chmod.split())
        try:
            # When "-no-update-model-required" is specified, processing a read-only measurement set should be possible.
            s = (
                f"{tcf.WSCLEAN} -name facet-readonly-ms -interval 10 20 "
                "-no-update-model-required -auto-threshold 0.5 -auto-mask 3 "
                "-mgain 0.95 -nmiter 2 -multiscale -niter 100000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} "
                f"{tcf.DIMS_SMALL} {tcf.MWA_MOCK_FULL}"
            )
            validate_call(s.split())
        finally:
            chmod = f"chmod u+w -R {tcf.MWA_MOCK_FULL}"
            validate_call(chmod.split())

    @pytest.mark.parametrize("mpi", [False, True])
    def test_facetbeamimages(self, mpi, tmp_mwa_mock_facet):
        """
        Basic checks of the generated images when using facet beams. For each image,
        test that the pixel values are valid (not NaN/Inf) and check the percentage
        of zero pixels.
        """

        deconvolve_facets(tmp_mwa_mock_facet, "wgridder", True, mpi, True)

        basic_image_check("facet-imaging-reorder-psf.fits")
        basic_image_check("facet-imaging-reorder-dirty.fits")

    def test_multi_channel(self):
        # Test for issue 122. Only test if no crash occurs.
        validate_call(
            (
                f"{tcf.WSCLEAN} -name multi-channel-faceting "
                "-parallel-gridding 3 -channels-out 2 "
                "-pol xx,yy -join-polarizations "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} "
                "-interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 "
                f"{tcf.MWA_MOCK_MS}"
            ).split()
        )

    def test_diagonal_solutions(self):
        validate_call(
            (
                f"{tcf.WSCLEAN} -name faceted-diagonal-solutions "
                "-parallel-gridding 3 -channels-out 2 "
                "-diagonal-solutions "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} "
                "-interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 "
                f"{tcf.MWA_MOCK_MS}"
            ).split()
        )

    def test_diagonal_solutions_with_beam(self):
        validate_call(
            (
                f"{tcf.WSCLEAN} -name faceted-diagonal-solutions "
                "-parallel-gridding 3 -channels-out 2 "
                "-diagonal-solutions -mwa-path . -apply-facet-beam "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} "
                "-interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 "
                f"{tcf.MWA_MOCK_MS}"
            ).split()
        )

    def test_parallel_gridding(self):
        """
        Run a single gridding cycle (no deconvolution / degridding).
        Compare serial, threaded and mpi run for facet based imaging
        with h5 corrections. Number of used threads/processes is
        deliberately chosen smaller than the number of facets.
        """
        names = [
            "facets-h5-serial",
            "facets-h5-threaded",
            "facets-h5-mpi",
            "facets-h5-hybrid",
        ]
        # Using only 2 threads/gridder yields relatively stable results.
        wsclean_commands = [
            f"{tcf.WSCLEAN} -j 2",
            f"{tcf.WSCLEAN} -j 6 -parallel-gridding 3",
            f"mpirun -np 3 {tcf.WSCLEAN_MP} -j 2 -max-mpi-message-size 42k",
            f"mpirun -np 3 {tcf.WSCLEAN_MP} -j 6 -parallel-gridding 3",
        ]
        for name, command in zip(names, wsclean_commands):
            s = (
                f"{command} -name {name} "
                "-pol xx,yy -join-polarizations "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} "
                f"-interval 10 14 {tcf.MWA_MOCK_MS}"
            )
            validate_call(s.split())

            # All images will be compared against the first image.
            # For the first image itself, only test whether the image is finite.
            if name == names[0]:
                rms = compute_rms(f"{names[0]}-YY-image.fits")
                assert np.isfinite(rms)
            else:
                # Typical rms difference is about 1.0e-7
                threshold = 3.0e-7
                compare_rms_fits(
                    f"{names[0]}-YY-image.fits",
                    f"{name}-YY-image.fits",
                    threshold,
                )

    @pytest.mark.parametrize("compound_tasks", [False, True])
    def test_parallel_predict(
        self, compound_tasks, tmp_path, tmp_mwa_mock_facet
    ):
        """
        Run a single predict/degridding cycle (no deconvolution / gridding).
        Compare serial, threaded, mpi and hybrid runs.
        Do all parallel runs with and without enabling compound tasks.
        """
        names = ["threaded", "mpi", "hybrid"]
        wsclean_commands = [
            f"{tcf.WSCLEAN} -j 3 -parallel-gridding 3",
            f"mpirun -np 3 {tcf.WSCLEAN_MP} -max-mpi-message-size 42k",
            f"mpirun -np 3 {tcf.WSCLEAN_MP} -j 3 -parallel-gridding 3",
        ]

        # Create reference output using a basic sequential run.
        predict_facet_image(tmp_mwa_mock_facet)

        # Run various alternatives and compare output against the reference.
        for name, command in zip(names, wsclean_commands):
            name = "test_" + name + "_degridding"

            if compound_tasks:
                name += "_compound"
                command += " -compound-tasks"

            ms = tmp_path / name
            shutil.copytree(tcf.MWA_MOCK_FACET, ms)
            predict_facet_image(ms, wsclean_command=command)
            assert_taql(
                f"select from {tmp_mwa_mock_facet} t1, {ms} t2 "
                "where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
            )

    def test_compound_tasks(self):
        """
        Run a single gridding cycle (no deconvolution / degridding).
        Compares a basic serial run without compound tasks to
        runs with compound tasks.
        """
        names = [
            "facets-h5-nocompound-sequential",
            "facets-h5-compound-sequential",
            "facets-h5-compound-threaded",
            "facets-h5-compound-sequential-mpi-local",
            "facets-h5-compound-threaded-mpi-remote",
        ]
        # Because of the static channel-to-node map, using more than
        # 2 processes makes no sense: This test only has a single channel.
        # The MPI tests either run everything 'local'ly or 'remote'ly.
        mpi_cmd = f"mpirun -np 2 {tcf.WSCLEAN_MP}"
        # Using 5 tasks/node makes the main node send the compound tasks for
        # the yy polarization while the task for xx is not yet finished
        # Using only 1 thread/gridder yields very stable results: It allows
        # using zero tolerance when comparing sequential runs (see below).
        pg = "-j 5 -parallel-gridding 5"
        wsclean_commands = [
            f"{tcf.WSCLEAN} -j 1",
            f"{tcf.WSCLEAN} -j 1 -compound-tasks",
            f"{tcf.WSCLEAN} {pg} -compound-tasks",
            f"{mpi_cmd} -j 1 -compound-tasks",
            f"{mpi_cmd} {pg} -compound-tasks -no-work-on-master",
        ]
        for name, command in zip(names, wsclean_commands):
            s = (
                f"{command} -name {name} "
                "-pol xx,yy -join-polarizations "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} "
                f"-interval 10 14 {tcf.MWA_MOCK_MS}"
            )
            validate_call(s.split())

            # All images will be compared against the first image.
            # For the first image itself, only test whether the image is finite.
            if name == names[0]:
                rms = compute_rms(f"{names[0]}-YY-image.fits")
                assert np.isfinite(rms)
            else:
                # Pure sequential tests should produce equal results.
                # In parallel tests, typical RMS difference is about 1.0e-7.
                threshold = 3.0e-7 if pg in command else 0.0
                compare_rms_fits(
                    f"{names[0]}-YY-image.fits",
                    f"{name}-YY-image.fits",
                    threshold,
                )

    @pytest.mark.parametrize("beam", [False, True])
    @pytest.mark.parametrize(
        "h5file",
        [
            None,
            [tcf.MOCK_SOLTAB_2POL],
            [tcf.MOCK_SOLTAB_2POL, tcf.MOCK_SOLTAB_2POL],
        ],
    )
    def test_multi_ms(self, beam, h5file):
        """
        Check that identical images are obtained in case multiple (identical) MSets and H5Parm
        files are provided compared to imaging one MSet
        """
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
            s = f"{tcf.WSCLEAN} -j 1 -nmiter 2 -gridder wgridder -name {name} -facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} -interval 10 14 -niter 1000000 -auto-threshold 5 -mgain 0.8 {command}"
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

    def test_diagonal_solutions(self):
        # Initialize random rumber generator
        rng = np.random.default_rng(1)

        # Strip unused stations from mock measurement set
        s = f"DP3 msin={tcf.MWA_MOCK_MS}  msout=diagonal_solutions.ms msout.overwrite=True steps=[filter] filter.remove=True"
        validate_call(s.split())

        # Fill WEIGHT_SPECTRUM with random values
        with casacore.tables.table(
            "diagonal_solutions.ms", readonly=False
        ) as t:
            weight_spectrum_shape = np.concatenate(
                (
                    np.array([t.nrows()]),
                    t.getcoldesc("WEIGHT_SPECTRUM")["shape"],
                )
            )
            weights = rng.uniform(0, 1, weight_spectrum_shape) + np.array(
                [1, 2, 3, 4], ndmin=3
            )
            t.putcol("WEIGHT_SPECTRUM", weights)

        # Create a template image
        s = (
            f"{tcf.WSCLEAN} -gridder wgridder -name template-diagonal-solutions "
            f"{tcf.DIMS_SMALL} -interval 0 1 diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Use template image to create a sky model consisting of a grid of point sources
        with fits.open("template-diagonal-solutions-image.fits") as f:
            wcs = WCS(f[0].header)
            nr_pixels = f[0].shape[-1]
        pointsource_grid_size = 2
        source_positions = create_pointsource_grid_skymodel(
            "diagonal-solutions-skymodel.txt",
            pointsource_grid_size,
            nr_pixels,
            wcs,
        )

        # Predict (without solutions)
        s = f"DP3 msin=diagonal_solutions.ms msout= steps=[predict] predict.sourcedb=diagonal-solutions-skymodel.txt"
        validate_call(s.split())

        # Image (without solutions)
        s = (
            f"{tcf.WSCLEAN} -name diagonal-solutions-reference -no-reorder "
            f"{tcf.DIMS_SMALL} diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Create template solutions .h5 file
        s = "DP3 msin=diagonal_solutions.ms msout= steps=[ddecal] ddecal.sourcedb=diagonal-solutions-skymodel.txt ddecal.h5parm=diagonal-solutions.h5 ddecal.mode=complexgain"
        validate_call(s.split())

        # Fill the template solutions file with random data
        with h5py.File("diagonal-solutions.h5", mode="r+") as f:
            f["sol000"]["phase000"]["val"][:] = rng.uniform(
                -np.pi, np.pi, f["sol000"]["phase000"]["val"].shape
            )
            f["sol000"]["phase000"]["weight"][:] = 1.0
            f["sol000"]["amplitude000"]["val"][:] = rng.uniform(
                0.5, 3, f["sol000"]["amplitude000"]["val"].shape
            )
            f["sol000"]["amplitude000"]["weight"][:] = 1.0

        # Predict with (random) solutions
        s = (
            "DP3 msin=diagonal_solutions.ms msout= steps=[h5parmpredict] "
            "h5parmpredict.sourcedb=diagonal-solutions-skymodel.txt "
            "h5parmpredict.applycal.parmdb=diagonal-solutions.h5 "
            "h5parmpredict.applycal.steps=[ampl,phase] "
            "h5parmpredict.applycal.ampl.correction=amplitude000 "
            "h5parmpredict.applycal.phase.correction=phase000 "
            "h5parmpredict.applycal.correction=amplitude000"
        )
        validate_call(s.split())

        # Image data predicted with solutions applied,
        # without applying corrections for the solutions while imaging
        s = (
            f"{tcf.WSCLEAN} -name diagonal-solutions-no-correction -no-reorder "
            f"{tcf.DIMS_SMALL} diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Image data predicted with solutions applied,
        # while applying corrections
        s = (
            f"{tcf.WSCLEAN} -name diagonal-solutions -no-reorder "
            "-parallel-gridding 3 "
            f"{tcf.DIMS_SMALL} -mgain 0.8 -threshold 10mJy -niter 10000 "
            f"-facet-regions {tcf.FACETFILE_4FACETS} "
            "-apply-facet-solutions diagonal-solutions.h5 "
            "amplitude000,phase000 -diagonal-solutions "
            "diagonal_solutions.ms"
        )
        validate_call(s.split())

        # Compare reference, uncorrection and corrected fluxes
        reference_image_data = fits.getdata(
            "diagonal-solutions-reference-image.fits"
        )[0, 0]
        no_correction_image_data = fits.getdata(
            "diagonal-solutions-no-correction-image.fits"
        )[0, 0]
        image_data = fits.getdata("diagonal-solutions-image-pb.fits")[0, 0]
        # loop over input sources
        for idx0, idx1 in source_positions:
            # Assert that without corrections less than 5 percent flux is recovered
            assert np.abs(no_correction_image_data[idx0, idx1]) < 5e-2
            # Assert that with corrections the recovered flux is within 2 percent of the reference
            assert np.isclose(
                reference_image_data[idx0, idx1],
                image_data[idx0, idx1],
                rtol=2e-2,
            )

    def test_dd_psfs_with_faceting(self):
        validate_call(
            (
                f"{tcf.WSCLEAN} -name dd-psfs-with-faceting "
                f"-dd-psf-grid 3 3 -parallel-gridding 5 {tcf.DIMS_SMALL} "
                "-parallel-deconvolution 100 -channels-out 2 -join-channels "
                "-niter 100 -mgain 0.8 -apply-facet-beam -mwa-path . "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.MWA_MOCK_MS}"
            ).split()
        )
        import os.path

        basic_image_check("dd-psfs-with-faceting-MFS-image.fits")
        for i in range(9):
            assert os.path.isfile(
                f"dd-psfs-with-faceting-d000{i}-0000-psf.fits"
            )
            assert os.path.isfile(
                f"dd-psfs-with-faceting-d000{i}-0001-psf.fits"
            )
            assert os.path.isfile(
                f"dd-psfs-with-faceting-d000{i}-MFS-psf.fits"
            )
        assert not os.path.isfile(f"dd-psfs-with-faceting-0000-psf.fits")
        assert not os.path.isfile(f"dd-psfs-with-faceting-0001-psf.fits")
        assert not os.path.isfile(f"dd-psfs-with-faceting-MFS-psf.fits")

    def test_predict_with_solutions(self):
        # This is a more advanced prediction run which at some point failed
        shutil.copyfile(
            "point-source-model.fits", "point-source-0000-model-fpb.fits"
        )
        shutil.copyfile(
            "point-source-model.fits", "point-source-0001-model-fpb.fits"
        )
        validate_call(
            (
                f"{tcf.WSCLEAN} -name point-source -v -predict -reorder "
                "-parallel-gridding 4 -channels-out 2 -diagonal-solutions "
                "-apply-facet-beam -facet-beam-update 60 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-mwa-path . {tcf.MWA_MOCK_FACET}"
            ).split()
        )

    def test_facet_continuing(self):
        nthreads = 4
        s = (
            f"{tcf.WSCLEAN} -parallel-gridding {nthreads} "
            f"{tcf.DIMS_SMALL} -niter 100 -auto-threshold 5 -mgain 0.8 -channels-out 2 "
            f"-facet-regions {tcf.FACETFILE_4FACETS} "
            f"-name facet-continuing-a {tcf.MWA_MOCK_FULL}"
        )
        validate_call(s.split())
        s = (
            f"{tcf.WSCLEAN} -reuse-psf facet-continuing-a -reuse-dirty facet-continuing-a "
            f"-parallel-gridding {nthreads} {tcf.DIMS_SMALL} -niter 100 "
            f"-auto-threshold 5 -mgain 0.8 -channels-out 2 -facet-regions {tcf.FACETFILE_4FACETS} "
            f"-name facet-continuing-b -v {tcf.MWA_MOCK_FULL}"
        )
        validate_call(s.split())
        basic_image_check("facet-continuing-b-0000-dirty.fits")
        basic_image_check("facet-continuing-b-0000-image.fits")
        basic_image_check("facet-continuing-b-0000-psf.fits")
        basic_image_check("facet-continuing-b-0000-residual.fits")
        basic_image_check("facet-continuing-b-0001-dirty.fits")
        basic_image_check("facet-continuing-b-0001-image.fits")
        basic_image_check("facet-continuing-b-0001-psf.fits")
        basic_image_check("facet-continuing-b-0001-residual.fits")
