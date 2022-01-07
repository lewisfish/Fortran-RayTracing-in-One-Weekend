module ray_mod
    
    use vec3_class

    implicit none
    
    type :: ray_t
        type(vec3) :: orig, dir
        contains
        procedure :: origin => orig
        procedure :: direction
        procedure :: at
    end type ray_t

    interface ray_t
        module procedure init_ray
    end interface ray_t

    contains
    
    type(ray_t) function init_ray(origin, dir)

        implicit none

        type(vec3), intent(IN) :: origin, dir

        init_ray%orig = origin
        init_ray%dir = dir

        end function init_ray

    type(vec3) function orig(this)
        
        implicit none

        class(ray_t) :: this
        
        orig = this%orig

    end function orig

    type(vec3) function direction(this)
        
        implicit none

        class(ray_t) :: this
        
        direction = this%dir

    end function direction

    type(vec3) function at(this, t)
        
        implicit none

        class(ray_t)     :: this
        real, intent(IN) :: t
        
        at = this%orig + t * this%dir

    end function at

end module ray_mod