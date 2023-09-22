## Graphika

Using SDL2, MSVC, OpenMP, and SIMD, for rasterization of an .OBJ model. The list of possible raster methods used are listed bellow.

<p align="center">
  <img width="30%" src="images/plane_phong.PNG">
</p>

### Raster Methods
- [edging.h](src/graphics/raster_mthods/edging.h) - using the edge equation to determine if a pixel is inside or outside of the triangle
- [stepping.h](src/graphics/raster_mthods/stepping.h) - same idea as edging, but we dont calculate the edge equation for every pixel, instead, we calculate it for a starting pixel and increment it in the x and y direction to cover the entire triangle
- [avx.h](src/graphics/raster_mthods/avx.h) - using SIMD instructions to improve the performance of the stepping method
- [matrix.h](src/graphics/raster_mthods/matrix.h) - [Incremental and Hierarchical Hilbert Order Edge Equation Polygon Rasterization](https://doi.org/10.1145/383507.383528)
- [phong.h](src/graphics/raster_mthods/phong.h) - Phong lighting implemented using the avx method for rasterising
- [normal_map.h](src/graphics/raster_mthods/normal_map.h) - Normal Mapping implemented using the avx method for rasterising
- [parallax.h](src/graphics/raster_mthods/edging.h)- Different Parallax mapping techniques implemented using the avx method; Offset, Steep, Relief and Parallax Occlusion Mapping (POM)  
