module hittable_list_mod
    
    use material_mod

    implicit none
    
    type, extends(hittable) :: hittable_list
        type(container), allocatable :: list(:)
        contains
        procedure :: hit => hit_list
    end type hittable_list

    type :: container
        class(hittable), pointer :: p => null()
    end type container

    contains
    
    logical function hit_list(this, ray, t_min, t_max, rec)
        
        implicit none
        
        class(hittable_list) :: this
        
        type(ray_t), intent(IN) :: ray
        real,        intent(IN) :: t_min, t_max
        type(hit_record), intent(OUT)        :: rec
        
        type(hit_record) :: temp_rec
        logical :: hit_anything
        real :: closest_so_far
        integer :: i

        hit_anything = .false.
        closest_so_far = t_max

        do i = 1, size(this%list)
            if(this%list(i)%p%hit(ray, t_min, closest_so_far, temp_rec))then
                hit_anything = .true.
                closest_so_far = temp_rec%t
                rec = temp_rec
            end if
        end do

        hit_list = hit_anything

    end function hit_list

end module hittable_list_mod