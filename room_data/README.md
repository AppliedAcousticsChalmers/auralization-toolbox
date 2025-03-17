# The Input Data

The data to be auralized need to be in a format that depends on whether the data are available on a volumetric or on a surface grid. In either case, the origin of the coordinate system needs to be located in the center of the grid, and the coordinate of the sampling points need to be specified in meters in a right hand coordinate system where the `x`-`y`-plane is the horizontal plane and `z` points upwards. Our data use the following file naming conventions:

* A file name that starts with `sound_field_p_` indicates that the pressure is available on a volumetric grid.
* A file name that starts with `sound_field_pv_` indicates that pressure and velocity are available on a surface layer.
* A file name that starts with `sound_field_pp_` indicates that the pressure is available on a surface double layer.

In the present folder, we provide example sound field data for a reverberant room with a reverb decay time of around 1 s ("big hall") as well as for a much dryer living-room sized room ("listening lab") in the files. We computed our example data with the scripts from [this repository](https://github.com/AppliedAcousticsChalmers/acoustic-room-responses).

Note that the definition of the sampling grid (for example the file `resources/grid_spherical_surface_L25.mat`) needs to comprise the surface normals (the variable `normal_vector`) in the case of single-layer cubical grids. The surface normals are computed automatically for all other grid types if they are required. 

## Note on the Coordinate system

As mentioned above, the computation of the auralization matrix requires that the sampling grid is centered around the coordinate origin. The positive `x`-axis corresponds to the direction 'straight ahead'. The coordinate system based on which the simulation data to be auralized were computed may be such that the center point of the auralization sampling grid does not coincide with the coordinate origin. It may be shifted and rotated. This is not a problem so long as the relative positions of the sampling points are the same like in the computation of the auralization matrix and as long as the sampling points are specified in the same order. 

It may be useful store the coordinates of the sampling points together with the acoustic simulation data (even though these coordinates are not used at the current stage).  

## Volumetric Data

Type `data = load('sound_field_p_cubical_volume_listening_lab_L125.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
  
                 fs: 48000
         grid_shape: 'cubical_volume'
           pressure: [7024×125 single]
    sampling_points: [3×125 double]
```

`fs` is the sampling frequency in Hertz of the `pressure`, which are the pressure signals at the 125 points the coordinates of which are specified by `sampling points`. The sampled sound field happens to be a room impulse response with a length of 7024 samples.

## Surface Data

Surface data can be processed for two different geometries (`spherical_surface` and `cubical_surface`) and in two formats. The two formats are explained below.

### Single-Layer Pressure and Normal Particle Velocity

Type `data = load('sound_field_pv_spherical_surface_big_hall_L25.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
                 fs: 48000
         grid_shape: 'spherical_surface'
      normal_vector: [3×25 double]
           pressure: [51024×25 single]
    sampling_points: [3×25 double]
           velocity: [51024×25 single]
```

Here, both sound `pressure` and particle `velocity` are available at the 25 `sampling points`. The outward facing `normal_vector` at the sampling points is redudant for spherical surfaces but not for cubical ones. 

### Double-Layer Pressure

Type `data = load('sound_field_pp_spherical_surface_listening_lab_L25.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
  
                       fs: 48000
               grid_shape: 'spherical_surface'
           pressure_inner: [7024×25 single]
           pressure_outer: [7024×25 single]
    sampling_points_inner: [3×25 double]
    sampling_points_outer: [3×25 double]
```

Here, the processing pipeline takes as input the sound pressure signals on a double-layer surface with an inner and an outer surface. Both the coordinates of the sampling points as well as the pressure signals need to be specified in the same order. I.e., the 3rd column of `pressure_inner` and the 3rd column of `pressure_outer` form a pair of points that are specified by the 3rd columns of `sampling_points_inner` and `sampling_points_outer`. The inner and outer sampling points have to be such that a vector that points from an inner point to the corresponding outer point represents the local surface normal.

Informal studies that we conducted suggest that a double-leyer distance not larger than 5 mm makes the binaural output signals identical to the single-layer data. Larger distances may also be ok if some deviations of the binaural signals at very high frequencies (i.e. > 15 kHz) are tolerable.

