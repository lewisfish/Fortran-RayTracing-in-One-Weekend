module material_mod

    use vec3_class
    use ray_mod

    implicit none
    
    type, abstract :: material
        contains
        procedure(matInterface), deferred :: scatter
    end type material
    
    type :: hit_record
        type(vec3) :: p, normal
        real       :: t
        logical    :: front_face
        class(material), pointer :: mat_ptr => null()
        contains
        procedure :: set_face_normal
    end type hit_record

    abstract interface
        logical function matInterface(this, ray_in, rec, attenuation, scattered)
            use vec3_class
            use ray_mod
            import material, hit_record
            implicit none
            class(material)  :: this
            type(hit_record), intent(IN)  :: rec
            type(ray_t),      intent(IN)  :: ray_in
            type(vec3),       intent(OUT) :: attenuation
            type(ray_t),      intent(OUT) :: scattered
        end function matInterface
    end interface

    type, abstract :: hittable
        contains
        procedure(hitInterface), deferred :: hit
    end type hittable
    
    abstract interface
        logical function hitInterface(this, ray, t_min, t_max, rec)
            use vec3_class
            use ray_mod
            import hittable, hit_record
            implicit none
            class(hittable) :: this
            real,             intent(IN)  :: t_min, t_max
            type(ray_t),      intent(IN)  :: ray
            type(hit_record), intent(OUT) :: rec
        end function hitInterface
    end interface


!===============================================================

        type, extends(material) :: lambertian
            type(vec3) :: albedo
            contains
            procedure :: scatter => scatter_lamb
        end type lambertian

        type, extends(material) :: metal
            type(vec3) :: albedo
            real :: fuzz
            contains
            procedure :: scatter => scatter_metal
        end type metal

        interface metal
            module procedure init_metal
        end interface metal

        type, extends(material) :: dielectric
            real :: ir
            contains
            procedure :: scatter => scatter_dielectric
        end type dielectric


    contains

        logical function scatter_lamb(this, ray_in, rec, attenuation, scattered)
            
            implicit none

            class(lambertian) :: this

            type(ray_t), intent(IN)      :: ray_in
            type(ray_t), intent(OUT)     :: scattered
            type(vec3), intent(OUT)      :: attenuation
            type(hit_record), intent(IN) :: rec

            type(vec3) :: scatter_direction

            scatter_direction = rec%normal + random_unit_vector()
            !catch degenerate scatter direction
            if(scatter_direction%near_zero())scatter_direction = rec%normal

            scattered = ray_t(rec%p, scatter_direction)
            attenuation = this%albedo
            scatter_lamb = .true.

        end function scatter_lamb


        type(metal) function init_metal(albedo, f)
            
            implicit none
            
            type(vec3), intent(IN) :: albedo
            real,       intent(IN) :: f
            
            init_metal%albedo = albedo
            if(f < 1.)then
                init_metal%fuzz = f
            else
                init_metal%fuzz = 1.
            end if

        end function init_metal


        logical function scatter_metal(this, ray_in, rec, attenuation, scattered)
            
            implicit none

            class(metal) :: this

            type(ray_t),      intent(IN)  :: ray_in
            type(ray_t),      intent(OUT) :: scattered
            type(vec3),       intent(OUT) :: attenuation
            type(hit_record), intent(IN)  :: rec

            type(vec3) :: reflected, dir

            dir = ray_in%direction()
            dir = dir%magnitude()
            reflected = reflect(dir, rec%normal)
            scattered = ray_t(rec%p, reflected + this%fuzz * random_in_unit_sphere())
            attenuation = this%albedo
            scatter_metal = (scattered%direction() .dot. rec%normal) > 0.

        end function scatter_metal

        logical function scatter_dielectric(this, ray_in, rec, attenuation, scattered)
            
            use vec3_class, only : refract
            use utils,      only : random

            implicit none

            class(dielectric) :: this

            type(ray_t),      intent(IN)  :: ray_in
            type(ray_t),      intent(OUT) :: scattered
            type(vec3),       intent(OUT) :: attenuation
            type(hit_record), intent(IN)  :: rec

            type(vec3) :: unit_direction, direction
            real       :: refraction_ratio, cos_theta, sin_theta
            logical    :: cannot_refract

            attenuation = vec3(1., 1., 1.)
            if(rec%front_face)then
                refraction_ratio = 1./this%ir
            else
                refraction_ratio = this%ir
            end if

            unit_direction = ray_in%direction()
            unit_direction = unit_direction%magnitude()

            cos_theta = min((-1.)*unit_direction .dot. rec%normal, 1.0)
            sin_theta = sqrt(1.0 - cos_theta*cos_theta)

            cannot_refract = refraction_ratio * sin_theta > 1.0

            if(cannot_refract .or. reflectance(cos_theta, refraction_ratio) > random())then
                direction = reflect(unit_direction, rec%normal)
            else
                direction = refract(unit_direction, rec%normal, refraction_ratio)
            end if
            scattered = ray_t(rec%p, direction)
            scatter_dielectric = .true.

        end function scatter_dielectric

        real function reflectance(cosine, ref_idx)
            
            implicit none
            
            real, intent(IN) :: cosine, ref_idx
            
            real :: r0

            r0 = (1. - ref_idx) / (1. + ref_idx)
            r0 = r0**2
            reflectance = r0 + (1-r0) * (1.-cosine)**5

        end function reflectance

        subroutine set_face_normal(this, ray, outward_normal)
        
        implicit none
        
        class(hit_record) :: this
        
        type(ray_t), intent(IN) :: ray
        type(vec3),  intent(IN) :: outward_normal
    
        this%front_face = (ray%direction() .dot. outward_normal) < 0.

        if(this%front_face)then
            this%normal = outward_normal
        else
            this%normal = (-1.)*outward_normal
        end if

    end subroutine set_face_normal


end module material_mod