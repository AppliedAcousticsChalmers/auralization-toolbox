# The Input Data

The data to be auralized need to be in a format that depends on whether the data are available on a volumetric or on a surface grid. In either case, the origin of the coordinate system needs to be located in the center of the grid, and the coordinate of the sampling points need to be specified in meters in a right hand coordinate system where the `x`-`y`-plane is the horizontal plane and `z` points upwards. 

* A file name that starts with `sound_field_p_` indicates that the pressure is available on a volumetric grid.
* A file name that starts with `sound_field_pv_` indicates that pressure and velocity are available on a surface layer.
* A file name that starts with `sound_field_pp_` indicates that the pressure is available on a surface double layer.



## Volumetric Data

Type `data = load('sound_field_p_cubical_volume_big_hall.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
  
                 fs: 48000
         grid_shape: 'cubical_volume'
           pressure: [50000×343 single]
               room: 'big_hall'
    sampling_points: [3×343 double]
```

`fs` is the sampling frequency in Hertz of the `pressure`, which are the pressure signals at the 343 points the coordinates of which are specified by `sampling points`. The sampled sound field happens to be a room impulse response with a length of 50,000 samples.

## Surface Data

Surface data can be processed in two formats:

### Single-Layer Pressure and Normal Particle Velocity

Type `data = load('sound_field_pv_spherical_surface_living_room.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
                 fs: 48000
         grid_shape: 'spherical_surface'
      normal_vector: [3×100 double]
           pressure: [4828×100 single]
               room: 'living_room'
    sampling_points: [3×100 double]
           velocity: [4828×100 single]
```

Here, both sound `pressure` and particle `velocity` are available at the `sampling points`. The outward facing `normal_vector` at the sampling points is redudant for spherical surfaces but not for cubical ones. 

### Double-Layer Pressure

Type `data = load('sound_field_pp_spherical_surface_living_room.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
  
                       fs: 48000
               grid_shape: 'spherical_surface'
           pressure_inner: [50000×100 single]
           pressure_outer: [50000×100 single]
                     room: 'big_hall'
    sampling_points_inner: [3×100 double]
    sampling_points_outer: [3×100 double]
```

Here, the processing pipeline takes as input the sound pressure signals on a double-layer surface with an inner and an outer surface. Both the coordinates of the sampling points as well as the pressure signals need to be specified in the same order. I.e., the 3rd column of `pressure_inner` and the 3rd column of `pressure_outer` form a pair of points that are specified by the 3rd columns of `sampling_points_inner` and `sampling_points_outer`. The inner and outer sampling points have to be such that a vector that points from an inner point to the corresponding outer point represents the local surface normal.

We do not have an robust data yet on what the distances between the layers can be meaningful. Anything up to a couple of mm will certainly work well. Large distances may work well, too.

