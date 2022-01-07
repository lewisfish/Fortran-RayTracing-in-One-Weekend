# Fortran RayTracing in One Weekend

Fortran implmentation of Peter Shirley's [RayTracing in One Weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html).

Uses [fortran TEV bindings](https://github.com/lewisfish/fortran_tev_bindings) to [TEV](https://github.com/Tom94/tev/) to show live image render updates.

Parallelised using openMP.

## Example Output

For image size 1200x800 with 500 samples.

![Example Output](https://github.com/lewisfish/Fortran-RayTracing-in-One-Weekend/raw/main/image.png)

## Building

Can be bulit using the [Fortran Package Manager](https://fpm.fortran-lang.org/en/index.html):

```
fpm run --profile release --flag -fopemp
```
