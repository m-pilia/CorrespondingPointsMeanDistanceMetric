## Description

This is a groupwise metric for the [Elastix](http://elastix.isi.uu.nl/)
registration tool that, given a set of landmarks present in a collection of
images, aligns the occurrences of each landmark to their mean position. It
targets the following energy
```
E(mu) = 1/n Σ_i 1/m Σ_j [ p_ij - T( 1/m Σ_k p_ik ) ]
```
where `p_ij` is the position of the i-th landmark in the j-th image, `n` is the
number of landmarks, `m` the number of images, `T` the transform with
parameters `mu`.

This metric assumes that a [`BSplineStackTransform`](http://elastix.isi.uu.nl/doxygen/classelastix_1_1BSplineStackTransform.html) is used in the registration.

## Structure of the moving points file

Only the moving points file is used, the fixed must be provided only
as a dummy to prevent errors, but it is not used in this metric.

It is important to store the moving points in the correct order.
Given `m` moving images (stored in `m` channels in the time axis),
the points must be sorted in the file as
```
x00 y00 z00 0  // image 0, landmark 0
x01 y01 z01 0  // image 0, landmark 1
     ...
x0n y0n z0n 0  // image 0, landmark n

x10 y10 z10 1  // image 1, landmark 0
x11 y11 z11 1  // image 1, landmark 1
     ...
x1n y1n z1n 1  // image 1, landmark n

  .........

xm0 ym0 zm0 m  // image m, landmark 0
xm1 ym1 zm1 m  // image m, landmark 1
    ...
xmn ymn zmn m  // image m, landmark n
```

## Build

In order to use this metric, it is necessary to include it in an Elastix build.
This can be done by passing to CMake the configuration parameter
`-DELASTIX_USER_COMPONENT_DIRS=$dir`, where `$dir` is the root directory of
this repository, when building Elastix.

## License

This software is licensed under the Apache License, Version 2.0. For the full
license text, see the attached
[LICENSE](https://github.com/m-pilia/CorrespondingPointsMeanDistanceMetric/blob/master/LICENSE)
file.

