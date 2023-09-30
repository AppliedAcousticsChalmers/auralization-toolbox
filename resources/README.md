# The Input Data

The data to be auralized need to be in a format that depends on whether the data are available on a volumetric or on a surface grid. 

## Volumetric data:

Type `data = load('sound_field_cubical_volume_big_hall.mat')` in the MATLAB command window to see the contents of the file:

```
data = 
  struct with fields:
  
                 fs: 48000
         grid_shape: 'cubical_volume'
           pressure: [50000×343 single]
               room: 'big_hall'
    sampling_points: [3×343 double]
```

`fs` is the sampling frequency in Hertz of `pressure`, which are the pressure signals at the 343 points the coordinates of which are specified in meters by `sampling points` in right hand coordinate system where the `x`-`y`-plane is the horizontal plane and `z` points upwards. 

## Surface data:

Type `data = load('sound_field_spherical_surface_big_hall.mat')` in the MATLAB command window to see the contents of the file:

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

As of now, it is not possible to use particle velocity as input data. Instead, the processing pipeline takes as input the sound pressure signals on a double-layer surface with an inner and an outer surface. Both the coordinates of the sampling points as well as the pressure signals need to be specified in the same order. I.e., the 3rd column of `pressure_inner` and the 3rd column of `pressure_outer` form a pair of points that are specified by the 3rd columns of `sampling_points_inner` and `sampling_points_outer`.

Specifying the particle velocity as input will be available in future releases.
