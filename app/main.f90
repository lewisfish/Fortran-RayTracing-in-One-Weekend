module main
    implicit none

    contains

    real function hit_sphere(centre, radius, ray)
        
        use vec3_class
        use ray_mod

        implicit none

        type(vec3),  intent(IN) :: centre
        real,        intent(IN) :: radius
        type(ray_t), intent(IN) :: ray
        
        real :: a, c, discrim, half_b
        type(vec3) :: oc, a_tmp

        oc = ray%origin() - centre
        a_tmp = ray%direction()
        a = a_tmp%length_squared()
        half_b = oc .dot. ray%direction()
        c = oc%length_squared() - radius**2
        discrim = half_b**2 - a*c

        if(discrim < 0.)then
            hit_sphere = -1.
        else
            hit_sphere = (-half_b - sqrt(discrim)) / a
        end if

    end function hit_sphere

    subroutine write_colour(colour, samples_per_pixel)

        use iso_fortran_env, only : output_unit
        use vec3_class
        use utils, only : clamp

        implicit none
    
        type(vec3), intent(IN) :: colour
        integer, intent(IN) :: samples_per_pixel

        real :: r, g, b, scale

        r = colour%x
        g = colour%y
        b = colour%z
    
        scale = 1. / real(samples_per_pixel)
        r = sqrt(r * scale)
        g = sqrt(g * scale)
        b = sqrt(b * scale)

        write(output_unit,"(3(I3.1, 1x),A)") int(256*clamp(r, 0., 0.999)),int(256*clamp(g, 0., 0.999)),int(256*clamp(b, 0., 0.999))
    
    end subroutine write_colour


    recursive function ray_color(ray, world, depth) result(res)

        use vec3_class
        use ray_mod
        use material_mod
        use hittable_list_mod
        use utils, only : infinity

        implicit none

        type(ray_t), intent(IN) :: ray
        integer :: depth
        type(hittable_list), intent(IN) :: world

        type(hit_record) :: rec
        type(vec3) :: unit_direction, res, attenuation
        type(ray_t) :: scattered
        real :: t

        if(depth <= 0)then
            res = vec3(0., 0., 0.)
            return
        end if

        if(world%hit(ray, 0.001, infinity, rec))then
            if(rec%mat_ptr%scatter(ray, rec, attenuation, scattered))then
                res = attenuation * ray_color(scattered, world, depth-1)
                return
            end if
        end if

        unit_direction = ray%direction()
        unit_direction = unit_direction%magnitude()
        t = 0.5 * (unit_direction%y + 1.)
        res = (1. - t)*vec3(1., 1., 1.) + t*vec3(.5, .7, 1.)

    end function ray_color

    function random_scene() result(world)
        
        use utils, only : random
        use hittable_list_mod
        use sphere_mod
        use material_mod

        implicit none

        type(hittable_list) :: world
        class(material), pointer :: sphere_material
        type(sphere), target, save, allocatable :: sph(:)
        type(dielectric), target, save :: material1(530)
        type(lambertian), target, save :: ground_material, material2(530)
        type(metal), target, save :: material3(530)
        type(vec3) :: centre, albedo, tmp_vec
        integer :: cnter, a, b, i, mat1, mat2, mat3
        real :: choose_mat, fuzz

        ground_material = lambertian(vec3(0.5, 0.5, 0.5))
        allocate(world%list(530), sph(530))
        sph(1) = sphere(vec3(0., -1000., 0.), 1000., ground_material)
        do i = 1, 530
            allocate(world%list(i)%p, source=sph(i))
        end do
        world%list(1)%p => sph(1)

        cnter = 2
        mat1 = 1
        mat2 = 1
        mat3 = 1

        do a = -11, 11
            do b = -11, 11
                choose_mat = random()
                centre = vec3(a + .9*random(), .2, b + 0.9*random())

                tmp_vec = centre - vec3(4., .2, .0)
                if(tmp_vec%length() > 0.9)then
                    if(choose_mat < 0.8)then
                        !diffuse
                        albedo = random_vec3() * random_vec3()
                        material2(mat2) = lambertian(albedo)
                        sphere_material => material2(mat2)
                        mat2 = mat2 + 1
                        sph(cnter) = sphere(centre, 0.2, sphere_material)
                        world%list(cnter)%p => sph(cnter)
                        cnter = cnter + 1
                    elseif(choose_mat < 0.95)then
                        ! metal
                        albedo = random_vec3(0.5, 1.)
                        fuzz = random(0., 0.5)
                        material3(mat3) = metal(albedo, fuzz)
                        sphere_material => material3(mat3)
                        mat3 = mat3 + 1
                        sph(cnter) = sphere(centre, 0.2, sphere_material)
                        world%list(cnter)%p => sph(cnter)
                        cnter = cnter + 1
                    else
                        !glass
                        material1(mat1) = dielectric(1.5)
                        sphere_material => material1(mat1)
                        mat1 = mat1 + 1
                        sph(cnter) = sphere(centre, 0.2, sphere_material)
                        world%list(cnter)%p => sph(cnter)
                        cnter = cnter + 1
                    end if
                end if
                if(cnter >= 528)goto 10
            end do
        end do
10 continue
    material1(mat1+1) = dielectric(1.5)
    sph(cnter) = sphere(vec3(.0, 1., 0.), 1., material1(mat1+1))
    world%list(cnter)%p => sph(cnter)
    cnter = cnter + 1

    material2(mat2+1) = lambertian(vec3(0.4, 0.2, 0.1))
    sph(cnter) = sphere(vec3(-4., 1., 0.), 1., material2(mat2+1))
    world%list(cnter)%p => sph(cnter)
    cnter = cnter + 1

    material3(mat3+1) = metal(vec3(.7, .6, .5), 0.0)
    sph(cnter) = sphere(vec3(4., 1., 0.), 1., material3(mat3+1))
    world%list(cnter)%p => sph(cnter)

    end function random_scene

end module main
program p
    
    use iso_fortran_env, only : output_unit, error_unit
    use main
    use vec3_class
    use ray_mod
    use sphere_mod
    use hittable_list_mod
    use camera_mod
    use material_mod
    use utils, only : random, init_rng
    use omp_lib
    use tev_mod, only : tevipc

    implicit none
    
    type(tevipc) :: tev
    integer     :: image_width, image_height, i, j, s, samples_per_pixel, max_depth, progress, id
    real        :: aspect_ratio, u, v, dist_to_focus, aperture,x, y, z
    type(vec3)  :: pixel_color, lookat, lookfrom, vup
    type(ray_t) :: r

    real, allocatable :: rgb(:,:,:)

    ! type(sphere), allocatable, target :: sph(:)
    type(hittable_list) :: world
    type(camera) :: cam
    ! type(lambertian) :: material_ground, material_centre
    ! type(dielectric) :: material_left
    ! type(metal)      :: material_right

    tev = tevipc()

    !image
    aspect_ratio = 3. / 2.
    image_width = 500
    image_height = int(image_width / aspect_ratio)
    progress = image_height-1
    allocate(rgb(0:image_height-1, 0:image_width-1, 3))
    samples_per_pixel = 50
    max_depth = 10

    !world
    call init_rng()
    world = random_scene()
    ! material_ground = lambertian(vec3(.8, .8, 0.))
    ! material_centre = lambertian(vec3(.1, .2, 0.5))
    ! material_left = dielectric(1.5)
    ! material_right = metal(vec3(.8, .6, 0.2), .0)

    ! allocate(world%list(5), sph(5))
    ! do i = 1, 5
    !     allocate(world%list(i)%p, source=sph(i))
    ! end do

    ! sph(1) = sphere(vec3( 0.0, -100.5, -1.0), 100.0, material_ground)
    ! sph(2) = sphere(vec3( 0.0,    0.0, -1.0),   0.5, material_centre)
    ! sph(3) = sphere(vec3(-1.0,    0.0, -1.0),   0.5, material_left)
    ! sph(4) = sphere(vec3(-1.0,    0.0, -1.0),  -0.4, material_left)
    ! sph(5) = sphere(vec3( 1.0,    0.0, -1.0),   0.5, material_right)

    ! do i = 1, 5
    !     world%list(i)%p => sph(i)
    ! end do

    !camera
    lookfrom = vec3(13., 2., 3.)
    lookat = vec3(0., 0., 0.)
    vup = vec3(0., 1., 0.)
    aperture = 0.1
    ! dist_to_focus_vec = lookfrom - lookat
    dist_to_focus = 10.0!dist_to_focus_vec%length()
    cam = camera(lookfrom, lookat, vup, 20., aspect_ratio, aperture, dist_to_focus)
    call tev%create_image("RT", image_width, image_height, ["R", "G", "B"], .true.)
    !render
    ! write(output_unit,"(A,A,I4,1x,I3,A)")"P3",new_line("a"), image_width,image_height, new_line("a"),"255",new_line("a")
!$omp parallel default(none) shared(tev,progress,cam, rgb, samples_per_pixel, image_height, image_width, max_depth, world)&
!$omp& private(i, j, s, u, v, pixel_color, r, x, y, z, id)
    id = 0!omp_get_thread_num()

!$omp do
    do j = 0, image_width-1
        ! if(id == 0)then
        !     write(error_unit,"(A,I3.1)", advance="No")achar(27)//"[1000D"//"Scanlines remaing:",progress
        !     flush(error_unit)
        ! end if
        !$omp atomic
            progress = progress - 1
        do i = image_height-1,0,-1

            pixel_color = vec3(0., 0., 0.)
            do s = 1, samples_per_pixel
                u = real(j + random()) / (image_height-1.)
                v = real(i + random()) / (image_width-1.)
                r = cam%get_ray(u, v)
                pixel_color = pixel_color + ray_color(r, world, max_depth)
            end do
            x = pixel_color%x
            y = pixel_color%y
            z = pixel_color%z

            !$omp atomic
                rgb(image_height-1-i,j,1) = rgb(image_height-1-i,j,1) + x
            !$omp atomic
                rgb(image_height-1-i,j,2) = rgb(image_height-1-i,j,2) + y
            !$omp atomic
                rgb(image_height-1-i,j,3) = rgb(image_height-1-i,j,3) + z

if(mod(progress,10) == 0)then
!$omp critical
if(id == 0)call tev%update_image("RT", (rgb(:,:,1:1)), ["R"], 0, 0, .true., .true.)
if(id == 0)call tev%update_image("RT", (rgb(:,:,2:2)), ["G"], 0, 0, .true., .true.)
if(id == 0)call tev%update_image("RT", (rgb(:,:,3:3)), ["B"], 0, 0, .true., .true.)
!$omp end critical
end if
            ! call write_colour(pixel_color, samples_per_pixel)
        end do
    end do
    !$OMP end  do
    !$OMP end parallel
    ! do i = image_height-1, 0, -1
    !     do j = 0, image_width-1
    !         pixel_color = vec3(rgb(j,i,1), rgb(j,i,2), rgb(j,i,3))
    !         call write_colour(pixel_color, samples_per_pixel)
    !     end do
    ! end do
    ! write(error_unit,*)new_line("a")
    call tev%update_image("RT", (rgb(:,:,1:1)), ["R"], 0, 0, .true., .true.)
    call tev%update_image("RT", (rgb(:,:,2:2)), ["G"], 0, 0, .true., .true.)
    call tev%update_image("RT", (rgb(:,:,3:3)), ["B"], 0, 0, .true., .true.)


end program p