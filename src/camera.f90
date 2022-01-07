module camera_mod
    
    use vec3_class
    use ray_mod

    implicit none

    type :: camera
        real       :: vfov, aspect_ratio, aperture, focus_dist
        real       :: theta, h, viewport_height, viewport_width
        type(vec3) :: lookfrom, lookat, vup
        real,       private :: lens_radius
        type(vec3), private :: origin, horizontal, vertical, lower_left_corner, u, v, w
        contains
        procedure :: get_ray
    end type camera

    interface camera
        module procedure init_camera
    end interface camera

    contains
    
    function init_camera(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, focus_dist) result(res)
        
        use utils, only : degrees_to_radians

        implicit none
        
        type(camera) :: res

        type(vec3), intent(IN) :: lookfrom, lookat, vup
        real,       intent(IN) :: vfov, aspect_ratio, aperture, focus_dist
                
        res%vfov = vfov
        res%aspect_ratio = aspect_ratio
        res%focus_dist = focus_dist

        res%theta = degrees_to_radians(vfov)
        res%h = tan(res%theta/2.)
        res%viewport_height = 2. * res%h
        res%viewport_width = res%aspect_ratio * res%viewport_height

        res%w = lookfrom - lookat
        res%w = res%w%magnitude()
        res%u = vup .cross. res%w
        res%u = res%u%magnitude()
        res%v = res%w .cross. res%u

        res%origin = lookfrom
        res%horizontal = res%focus_dist * res%viewport_width * res%u
        res%vertical = res%focus_dist * res%viewport_height * res%v
        res%lower_left_corner = res%origin - res%horizontal/2. - res%vertical/2. - res%focus_dist*res%w

        res%lens_radius = aperture / 2.
    end function init_camera

    type(ray_t) function get_ray(this, s, t)

        implicit none

        class(camera) :: this
        real, intent(IN) :: s, t
        
        type(vec3) :: rd, offset

        rd = this%lens_radius * random_in_unit_disk()
        offset = this%u * rd%x + this%v * rd%y

        get_ray = ray_t(offset + this%origin, this%lower_left_corner + s*this%horizontal + t*this%vertical - this%origin - offset)

    end function get_ray

end module camera_mod