if [[ "$1" == "" ]] ; then
    echo Syntax: various-tests.sh \<ms\>
    exit 1
fi

ms="$1"
dims="-size 1024 1024 -scale 1amin"
rectdims="-size 1536 1024 -scale 1amin"

# Make dirty image
wsclean -name test-dirty ${dims} ${ms}

# Clean a rectangular unpadded image
wsclean -name clean-rectangular -padding 1 -auto-threshold 5 -mgain 0.8 \
-niter 100000 ${rectdims} ${ms}

# Auto-masked multi-scale clean
wsclean -name multiscale-automasked -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 \
  -multiscale -niter 100000 ${rectdims} ${ms}

# Multiple intervals
wsclean -name intervals -intervals-out 3 ${rectdims} ${ms}

# Multiple intervals + multiple channels with some cleaning
wsclean -name intervals-and-channels -intervals-out 3 -channels-out 2 -niter 1000 -mgain 0.8 ${dims} ${ms}

# Multi-frequency Högbom clean, no parallel gridding
wsclean -name mfhogbom -channels-out 4 -join-channels -auto-threshold 3 -mgain 0.8 \
-niter 1000000 ${rectdims} ${ms}

# Multi-frequency Högbom clean with spectral fitting
wsclean -name mfhogbom-fitted -channels-out 4 -join-channels -parallel-gridding 4 -fit-spectral-pol 2 \
-auto-threshold 3 -mgain 0.8 -niter 1000000 ${rectdims} ${ms}

# Multi-frequency multi-scale clean with spectral fitting, pallel gridding & cleaning
wsclean -name mfms-fitted -channels-out 4 -join-channels -parallel-gridding 4 \
-parallel-deconvolution 1000 -fit-spectral-pol 2 -multiscale \
-auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 ${rectdims} ${ms}

# Save the list of components
wsclean -name mfms-components -save-source-list -channels-out 4 -join-channels -parallel-gridding 4 \
-fit-spectral-pol 2 -auto-threshold 0.5 -auto-mask 3 -mgain 0.8 -niter 1000000 -multiscale \
-parallel-deconvolution 1000 \
${dims} ${ms}

# Linear joined polarizations with 4 joined channels
wsclean -name linearpol -niter 1000000 -auto-threshold 3.0 -pol XX,YY,XY,YX \
-join-polarizations -join-channels -mgain 0.85 -channels-out 4 -parallel-gridding 16 \
${dims} ${ms}

# Image two timesteps
wsclean -name two-timesteps -niter 1000000 -auto-threshold 3.0  \
-intervals-out 2 -interval 20 22 -mgain 0.85 ${rectdims} ${ms}

# Stop on negative components
wsclean -name stop-on-negatives -stop-negative -niter 100000 ${rectdims} ${ms}

# Shift the image with w-stacking gridder
wsclean -mgain 0.8 -auto-threshold 5 -niter 1000000 -make-psf ${rectdims} -shift 08h09m20s -39d06m54s -name shift-ws -no-update-model-required ${ms}

# Shift the image with w-gridder
wsclean -mgain 0.8 -use-wgridder -auto-threshold 5 -niter 1000000 -make-psf ${rectdims} -shift 08h09m20s -39d06m54s -name shift-wg -no-update-model-required ${ms}
