# pyexr
pure python exr reader and converter

Currently supported are uncompressed .exr files with multiple channels, these can be converted to matlab/octave .mat files. You might be asking why this is any better than the OpenEXR bindings for python? The answer is that OpenEXR is barely supported on windows operating systems (at least I have failed in building these), and it was much easier to write this little python script.

<b>Dependencies</b>: scipy and numpy

<b>Usage</b>: python exr_to_mat.py [-show] <i>exrfile</i> [<i>matfile</i>]<br>
- <i>exrfile</i> might contain a glob pattern<br>
- if <i>matfile</i> is not given, it will be derived from the <i>exrfile</i> name.
