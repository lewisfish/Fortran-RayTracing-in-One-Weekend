module sphere_mod

    use material_mod
    use vec3_class

    implicit none
    
    type, extends(hittable) :: sphere
        real       :: radius
        type(vec3) :: centre
        class(material), pointer :: mat_ptr => null()
        contains
        procedure :: hit => hit_sphere
    end type sphere

    interface sphere
        module procedure init_sphere
    end interface sphere

    contains
    

    type(sphere) function init_sphere(centre, radius, mat)
        implicit none
        type(vec3), intent(IN) :: centre
        real,       intent(IN) :: radius
        class(material), target :: mat
        
        init_sphere%centre = centre
        init_sphere%radius = radius
        init_sphere%mat_ptr => mat

    end function init_sphere

    logical function hit_sphere(this, ray, t_min, t_max, rec)

        use vec3_class
        use ray_mod
        
        implicit none

        class(sphere) :: this

        real,             intent(IN)  :: t_min, t_max
        type(ray_t),      intent(IN)  :: ray
        type(hit_record), intent(OUT) :: rec
    
        type(vec3) :: oc, a_tmp, outward_normal
        real       :: a, c, discrim, half_b, root, sqrtd

        oc = ray%orig - this%centre
        a_tmp = ray%dir
        a = a_tmp%length_squared()
        half_b = oc .dot. ray%dir
        c = oc%length_squared() - this%radius**2
        discrim = half_b**2 - a*c

        if(discrim < 0.)then
            hit_sphere = .false.
            return
        end if

        sqrtd = sqrt(discrim)
        root = (-half_b - sqrtd) / a
        if(root < t_min .or. t_max < root)then
            root = (-half_b + sqrtd) / a
            if(root < t_min .or. t_max < root)then
                hit_sphere = .false.
                return
            end if
        end if

        rec%t = root
        rec%p = ray%at(rec%t)
        outward_normal = (rec%p - this%centre) / this%radius
        call rec%set_face_normal(ray, outward_normal)
        rec%mat_ptr => this%mat_ptr

        hit_sphere = .true.

    end function hit_sphere
end module sphere_mod