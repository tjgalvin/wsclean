Diverging clean
===============

It may occassionally happen that CLEAN diverges. This might cause it to produce extremely high or NaN model flux values, or it might get into a 'feedback' loop and keep doing major cleaning iterations without progressing (in this case, 'not converging' is probably a better term than 'diverging').

The problem is more common during :doc:`multi-scale clean <multiscale_cleaning>`, but can ocassionally occur during single-scale cleaning.

During multi-scale cleaning, in the command line output of WSClean this typically looks like the following:

.. code-block:: text

    Iteration 5978, scale 761 px : 923.32 mJy at 2148,2892
    Iteration 6100, scale 1521 px : -930.26 mJy at 4238,2897
    Iteration 708620, scale 48 px : -nan Jy at 0,0
    
What this means is that **before** iteration 708620, the deconvolution process has diverged. At the beginning of iteration 708620, the peak in the image is "``-NaN``", according to the log. When, instead of reporting NaN, it reports suddenly a very high value, e.g. 1e7 KJy, this has the same cause. Be aware that the reported position (in this case 0,0) is not the source of the problem; the iteration that caused the problem happened before 708620, but is not in the output. These components somehow became coupled and diverged.

Causes & solutions
------------------

Here are possible causes:

1. When using multi-scale clean, a common cause is that (very) large scales cause problems. As can be seen in the above output, this is likely here the case, since the iteration that caused problems (**Iteration 6100**) used a large scale (**1521 px**). You can use these methods to avoid the issue: 
    A. Don't clean with large scales. Using ``-multiscale-max-scales`` the number of scales can be limitted. Alternatively (but less preferred), using ``-multiscale-scales`` a custom list can be provided. If the image does not contain flux at these scales, you can leave the largest scales out without loss in image fidelity. Fewer scales will also speed op multi-scale.
    B. Stop before the problem occurs, by adding one or more stopping criteria. ``-nmiter`` can especially be helpful for this, but using ``-niter``, ``-auto-threshold`` etc can help as well.
2. Cleaning too deep. It occasionally occurs that cleaning into the noise causes divergence. This is simply fixed by using an appropriate auto-threshold level (e.g. ``-auto-threshold 5``, or ``-auto-threshold 1 -auto-mask 5``).
3. Having bright emission near/against the edge of the image. This can cause a ripple at the edge, and produce weird results. If this is a likely cause, then it is best to increase the image size. Be aware that when using the ``parallel-deconvolution`` option, the image has many more edges. Make sure that the parallel deconvolution "size" parameter is not too small.
4. Using a value of ``mgain`` that's too high. The default of 0.8 is quite stable, but if your PSF has high sidelobes, this value might still be too large. In that case, try lowering it to 0.7 or 0.6. Note that this increases the nr of major iterations, so it can slow down the imaging considerably.
