import os
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

# Control whether some expensive tests are skipped (every CI run) or run (weekly/nightly runs)
RUN_EXPENSIVE_TESTS = (
    os.environ.get("RUN_EXPENSIVE_TESTS", "False").lower() == "true"
)


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


def deconvolve_facets(
    ms, gridder, reorder, mpi, apply_beam=False, n_threads=4
):
    """
    Perform deconvolution with facets.

    Parameters
    ----------
    ms: str
        MeasurementSet name
    gridder: str
        Name of the gridder
    reorder: bool
        Enable reordering
    mpi: bool
        Runs WSClean in MPI mode
    apply_beam: bool
        Enables facet beam corrections
    n_threads: int
        Number of threads to use for -parallel-gridding

    Returns
    -------
    image_prefix: str
        Prefix of all the images created by WSClean
    """
    mpi_cmd = f"{tcf.MPIRUN} -tag-output -np {n_threads} {tcf.WSCLEAN_MP}"
    thread_cmd = f"{tcf.WSCLEAN} -parallel-gridding {n_threads}"
    reorder_ms = "-reorder" if reorder else "-no-reorder"
    facet_beam = "-mwa-path . -apply-facet-beam" if apply_beam else ""
    image_prefix = "facet-imaging-{gridder}"

    # Since IDG assumes square facets, the facets described in in FACETFILE_4FACETS
    # will be extended outside of the bound of the default dimensions (DIMS_SMALL)
    # and, hence, will result in segfaults when stitching the facets back together.
    image_dims = (
        "-size 320 320 -scale 4amin"
        if gridder == "facet-idg"
        else tcf.DIMS_SMALL
    )
    s = (
        f"{mpi_cmd if mpi else thread_cmd} -gridder {gridder} {reorder_ms} "
        f"{image_dims} -niter 1000000 -auto-threshold 5 -mgain 0.8 "
        f"-facet-regions {tcf.FACETFILE_4FACETS} {facet_beam} "
        f"-name {image_prefix} {reorder_ms} -v {ms}"
    )
    validate_call(s.split())

    return image_prefix


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
    "prepare_mock_ms",
    "prepare_model_image",
    "prepare_mock_soltab",
    "prepare_large_ms",
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

    def test_intervals(self):
        """
        Test that faceting + intervals works (this produced an error before 2025-11-12 due to a bug when
        -no-dirty was also specified).
        We have to use the full mwa MS in this case, because multiple intervals are needed and the
        mwa mock MSs have only one timestep.
        """
        s = (
            f"{tcf.WSCLEAN} -name facet-with-intervals -intervals-out 4 -no-dirty -pol iquv -join-polarizations "
            f"-facet-regions {tcf.FACETFILE_4FACETS} "
            f"{tcf.DIMS_SMALL} {tcf.MWA_MS}"
        )
        validate_call(s.split())

        basic_image_check("facet-with-intervals-t0000-I-image.fits")
        basic_image_check("facet-with-intervals-t0003-I-image.fits")
        basic_image_check("facet-with-intervals-t0000-V-image.fits")
        basic_image_check("facet-with-intervals-t0003-V-image.fits")

    # Test assumes that IDG and EveryBeam are installed
    @pytest.mark.parametrize(
        "gridder", ["wstacking", "wgridder", "idg", "facet-idg"]
    )
    def test_stitching(self, gridder):
        """Test stitching of the facets"""
        prefix = f"facet-stitch-{gridder}"
        s = [
            tcf.WSCLEAN,
            "-quiet",
            f"-gridder {gridder}",
            tcf.DIMS_SMALL,
            "" if (gridder == "idg") else "-pol XX,YY",
            f"-facet-regions {tcf.FACETFILE_2FACETS}",
            f"-name {prefix}",
            tcf.MWA_MOCK_MS,
        ]
        validate_call(" ".join(s).split())
        fpaths = (
            [prefix + "-dirty.fits", prefix + "-image.fits"]
            if (gridder == "idg")
            else [
                prefix + "-XX-dirty.fits",
                prefix + "-YY-dirty.fits",
                prefix + "-XX-image.fits",
                prefix + "-YY-image.fits",
            ]
        )
        check_and_remove_files(fpaths, remove=True)

    # FIXME: we should test wstacking and facet-idg here too
    # but it fails on the taql assertion
    @pytest.mark.parametrize("gridder", ["wgridder"])
    @pytest.mark.parametrize(
        "apply_facet_beam", ["without_facet_beam", "with_facet_beam"]
    )
    def test_predict(self, gridder, apply_facet_beam, tmp_mwa_mock_facet):
        """
        Test predict only run

        Parameters
        ----------
        gridder : str
            wsclean compatible description of gridder to be used.
        """

        do_apply_facet_beam = apply_facet_beam == "with_facet_beam"

        predict_facet_image(tmp_mwa_mock_facet, gridder, do_apply_facet_beam)

        # A numerical check can only be performed in case no DD effects were applied.
        if not do_apply_facet_beam:
            predict_full_image(tcf.MWA_MOCK_FULL, gridder)
            taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tmp_mwa_mock_facet} t2 where not all(near(t1.MODEL_DATA,t2.MODEL_DATA,5e-3))"
            assert_taql(taql_command)

    @pytest.mark.parametrize("gridder", ["wgridder", "facet-idg"])
    @pytest.mark.parametrize("reorder", ["without_reorder", "with_reorder"])
    @pytest.mark.parametrize("mpi", ["without_mpi", "with_mpi"])
    def test_facetdeconvolution(self, gridder, reorder, mpi):
        """
        Test facet-based deconvolution on model visibilties containing three point sources.

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

        # Copy the predicted visibilities to new MS, put them in the DATA column instead of the MODEL_DATA column, and remove MODEL_DATA.
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

        # if gridder == "wgridder":
        # Check that the predicted visibilities have been properly copied to the new MS.
        taql_command = f"select from {tcf.MWA_MOCK_FULL} t1, {tcf.MWA_MOCK_FACET} t2 where not all(near(t1.MODEL_DATA,t2.DATA, 5e-3))"
        assert_taql(taql_command)

        do_reorder = reorder == "with_reorder"
        use_mpi = mpi == "with_mpi"

        # IDG is not compatible with the parallel-gridding option,
        # hence, we disable it for test variants involving IDG.
        n_threads = 1 if "idg" in gridder else 4
        image_prefix = deconvolve_facets(
            tcf.MWA_MOCK_FACET,
            gridder,
            do_reorder,
            use_mpi,
            n_threads=n_threads,
        )

        # The following checks whether the initially predicted model data and final CLEAN model agree. It fails
        # for facet-idg, likely due to some residual artiacts in the visibilities, even though it was able to
        # deconvolve the three point sources succesfully. This test remains to check for any regressions.
        if gridder == "wgridder":
            taql_command = f"select from {tcf.MWA_MOCK_FACET} where not all(near(DATA,MODEL_DATA, 5e-2))"
            assert_taql(taql_command)

        # To check that (facet) decolvolution was succesful, we verify that the
        # noise in the image is sufficiently low.
        residual_image = f"{image_prefix}-residual.fits"
        rms_residual = compute_rms(residual_image)
        assert rms_residual < 1e-4

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

    @pytest.mark.parametrize("mpi", ["without_mpi", "with_mpi"])
    def test_facetbeamimages(self, mpi, tmp_mwa_mock_facet):
        """
        Basic checks of the generated images when using facet beams. For each image,
        test that the pixel values are valid (not NaN/Inf) and check the percentage
        of zero pixels.
        """

        use_mpi = mpi == "with_mpi"
        image_prefix = deconvolve_facets(
            tmp_mwa_mock_facet, "wgridder", True, use_mpi, True
        )

        basic_image_check(f"{image_prefix}-psf.fits")
        basic_image_check(f"{image_prefix}-dirty.fits")

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

    @pytest.mark.parametrize(
        "apply_facet_beam", ["without_facet_beam", "with_facet_beam"]
    )
    @pytest.mark.parametrize("polarization", ["pol_i", "pol_iquv"])
    @pytest.mark.parametrize(
        "dd_facets", ["with_dd_facets", "without_dd_facets"]
    )
    def test_shared_facet_reads_and_writes(
        self, apply_facet_beam, polarization, dd_facets
    ):
        if not RUN_EXPENSIVE_TESTS and not (
            polarization == "pol_iquv"
            and apply_facet_beam == "with_facet_beam"
            and dd_facets == "with_dd_facets"
        ):
            pytest.skip(reason="Skipping expensive test")

        names = [
            "facets-no-shared",
            "facets-shared-reads",
            "facets-shared-writes",
            "facets-shared-reads-and-writes",
        ]
        do_apply_facet_beam = apply_facet_beam == "with_facet_beam"
        do_all_polarization = polarization == "pol_iquv"
        polarization_settings = "-pol i"
        if do_all_polarization:
            polarization_settings = "-pol iquv -join-polarizations"
        name_suffix = "-no-beam"
        facet_beam = ""
        if do_apply_facet_beam:
            facet_beam = "-mwa-path . -apply-facet-beam"
            name_suffix = "-beam"
        dd_psf_settings = ""
        use_dd_facets = dd_facets == "with_dd_facets"
        if use_dd_facets:
            dd_psf_settings = "-dd-psf-grid 2 2"
            name_suffix = f"{name_suffix}-dd-psf"

        for name in names:
            shared_args = ""
            if name == names[1]:
                shared_args = "-shared-facet-reads"
            if name == names[2]:
                shared_args = "-shared-facet-writes"
            if name == names[3]:
                shared_args = "-shared-facet-reads -shared-facet-writes"

            s = (
                f"{tcf.WSCLEAN} -name {name}{name_suffix} "
                f"{shared_args} "
                f"{facet_beam} "
                f"{polarization_settings} "
                f"{dd_psf_settings} "
                "-parallel-gridding 3 "
                "-channels-out 3 -join-channels "
                "-no-update-model-required "
                f"-apply-facet-solutions {tcf.MOCK_SOLTAB_2POL} ampl000,phase000 "
                f"-facet-regions {tcf.FACETFILE_4FACETS} {tcf.DIMS_SMALL} "
                "-nmiter 3 -niter 20000 -auto-threshold 5 -mgain 0.8 "
                f"{tcf.MWA_MOCK_MS}"
            )
            validate_call(s.split())

            if name != names[0]:
                if not do_all_polarization:
                    threshold = 5.0e-6
                    if use_dd_facets:
                        threshold = 6.0e-3
                    if do_apply_facet_beam:
                        threshold = 9.0e-3
                    compare_rms_fits(
                        f"{names[0]}{name_suffix}-MFS-image.fits",
                        f"{name}{name_suffix}-MFS-image.fits",
                        threshold,
                    )
                else:
                    threshold = 9.0e-3
                    if do_apply_facet_beam:
                        threshold = 6.0e-2
                    for pol in ["I", "Q", "U", "V"]:
                        compare_rms_fits(
                            f"{names[0]}{name_suffix}-MFS-{pol}-image.fits",
                            f"{name}{name_suffix}-MFS-{pol}-image.fits",
                            threshold,
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
            f"{tcf.MPIRUN} -np 3 {tcf.WSCLEAN_MP} -j 2 -max-mpi-message-size 42k",
            f"{tcf.MPIRUN} -np 3 {tcf.WSCLEAN_MP} -j 6 -parallel-gridding 3",
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

    @pytest.mark.parametrize(
        "compound_tasks", ["without_compound_tasks", "with_compound_tasks"]
    )
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
            f"{tcf.MPIRUN} -np 3 {tcf.WSCLEAN_MP} -max-mpi-message-size 42k",
            f"{tcf.MPIRUN} -np 3 {tcf.WSCLEAN_MP} -j 3 -parallel-gridding 3",
        ]

        # Create reference output using a basic sequential run.
        predict_facet_image(tmp_mwa_mock_facet)
        use_compound_tasks = compound_tasks == "with_compound_tasks"

        # Run various alternatives and compare output against the reference.
        for name, command in zip(names, wsclean_commands):
            name = "test_" + name + "_degridding"

            if use_compound_tasks:
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
        mpi_cmd = f"{tcf.MPIRUN} -np 2 {tcf.WSCLEAN_MP}"
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

    @pytest.mark.parametrize("beam", ["without_beam", "with_beam"])
    @pytest.mark.parametrize(
        "h5file",
        [
            "no_h5file",
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
        use_beam = beam == "with_beam"
        if use_beam:
            commands = [
                "-mwa-path . -apply-facet-beam " + command
                for command in commands
            ]

        if h5file != "no_h5file":
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

    def test_time_frequency_smearing(self):
        das6_ms_path = "/var/scratch/offringa/Raw-RFI-Test-Set/processed/L2014581_SAP000_SB079_uv.ms"
        if not os.path.isdir(das6_ms_path):
            if '"das6"' in os.environ.get("CI_RUNNER_TAGS", []):
                pytest.fail("MS not available while running on das6")
            else:
                pytest.skip("MS not available (not on das6?)")

        # Create a small test data set from a high time frequency resolution data set
        s = [
            "DP3",
            "msin=" + das6_ms_path,
            "msin.baseline=[CR]S*&",
            "msin.ntimes=8",
            "msin.nchan=32",
            "msout=time-frequency-smearing.ms",
            "msout.overwrite=true",
            "steps=[]",
        ]
        validate_call(s)

        # Make template model image
        s = f"{tcf.WSCLEAN} -size 4800 4800 -scale 5asec time-frequency-smearing.ms"
        validate_call(s.split())

        # Fill model images with grid of point sources
        f_image = fits.open("wsclean-image.fits")
        image_size = f_image[0].data.shape[-1]
        GRID_SIZE_1D = 3
        point_source_spacing = image_size // GRID_SIZE_1D
        position_range_1d = (
            point_source_spacing // 2
            + point_source_spacing * np.arange(GRID_SIZE_1D)
        )
        f_image[0].data[:] = 0.0
        for i in position_range_1d:
            for j in position_range_1d:
                f_image[0].data[0, 0, i, j] = 1.0
        f_image.writeto("wsclean-model.fits", overwrite=True)

        # Predict visibilities for the 3x3 point source grid at high time frequency resolution
        s = f"{tcf.WSCLEAN} -predict -model-column DATA -size 4800 4800 -scale 5asec time-frequency-smearing.ms"
        validate_call(s.split())

        # Average to a lower time frequency resolution
        s = [
            "DP3",
            "msin=time-frequency-smearing.ms",
            "msout=time-frequency-smearing-averaged.ms",
            "msout.overwrite=true",
            "steps=[average]",
            "average.timestep=8",
            "average.freqstep=32",
        ]
        validate_call(s)

        # Predict at low time frequency resolution, without taking time frequency smearing into account
        s = f"{tcf.WSCLEAN} -predict -size 4800 4800 -scale 5asec time-frequency-smearing-averaged.ms"
        validate_call(s.split())

        # Check that the error is relatively high
        with casacore.tables.table("time-frequency-smearing-averaged.ms") as t:
            d1 = t.getcol("DATA")[:, :, 0]
            d2 = t.getcol("MODEL_DATA")[:, :, 0]
            r = d1 - d2
            # Threshold is rather arbitrary, determined by running the test and rounding
            # the result downwards to a 'nice' number
            assert np.sum(np.abs(r**2)) > 450

        # Create a 3x3 grid of facets, each facet covers one point source in the model image
        with open("facets.reg", "w") as f:
            print(
                "# Region file format: DS9 version 4.1\n"
                'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1\n'
                "fk5\n"
                "\n"
                "polygon(322.44872,22.02513,322.51426,19.37362,325.32867,19.42002,325.30977,22.07500)\n"
                "polygon(325.30977,22.07500,325.32867,19.42002,325.33176,19.41714,328.16107,19.41714,328.16416,19.42002,328.18306,22.07500)\n"
                "polygon(331.04483,22.02511,328.18306,22.07500,328.16416,19.42002,330.97927,19.37360)\n"
                "polygon(322.57796,16.71169,325.35012,16.74595,325.35312,16.74887,325.33176,19.41714,325.32867,19.42002,322.51426,19.37362)\n"
                "polygon(328.13971,16.74887,328.16107,19.41714,325.33176,19.41714,325.35312,16.74887)\n"
                "polygon(330.97927,19.37360,328.16416,19.42002,328.16107,19.41714,328.13971,16.74887,328.14271,16.74595,330.91557,16.71168)\n"
                "polygon(322.63967,14.05948,325.37381,14.09035,325.35012,16.74595,322.57796,16.71169)\n"
                "polygon(328.11903,14.09035,328.14271,16.74595,328.13971,16.74887,325.35312,16.74887,325.35012,16.74595,325.37381,14.09035)\n"
                "polygon(330.85384,14.05947,330.91557,16.71168,328.14271,16.74595,328.11903,14.09035)\n",
                file=f,
            )

        # Run predict again on the low resolution data, this time including time frequency smearing
        s = f"{tcf.WSCLEAN} -predict -size 4800 4800 -scale 5asec -facet-regions facets.reg -apply-time-frequency-smearing time-frequency-smearing-averaged.ms"
        validate_call(s.split())

        # Check that the error is relatively low this time
        with casacore.tables.table("time-frequency-smearing-averaged.ms") as t:
            d1 = t.getcol("DATA")[:, :, 0]
            d2 = t.getcol("MODEL_DATA")[:, :, 0]
            r = d1 - d2
            # Threshold is rather arbitrary, determined by running the test and rounding
            # the result upwards to a 'nice' number
            assert np.sum(np.abs(r**2)) < 10

    @pytest.mark.parametrize("gridder", ["wgridder", "facet-idg"])
    def test_predict_with_solutions(self, gridder):
        # This is a more advanced prediction run which at some point failed
        shutil.copyfile(
            "point-source-model.fits", "point-source-0000-model-fpb.fits"
        )
        shutil.copyfile(
            "point-source-model.fits", "point-source-0001-model-fpb.fits"
        )
        # The parallel-gridding option is not compatible with IDG.
        parallel_gridding_option = (
            "" if gridder == "facet-idg" else "-parallel-gridding 4"
        )
        validate_call(
            (
                f"{tcf.WSCLEAN} -gridder {gridder} -name point-source -v -predict -reorder "
                f"{parallel_gridding_option} -channels-out 2 -diagonal-solutions "
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
