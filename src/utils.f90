module utils
    
    implicit none

    real, parameter :: pi=4.*atan(1.), infinity=huge(1.)

    interface random
        module procedure random_uni
        module procedure random_uni_range
    end interface random

    contains
    
    real function degrees_to_radians(degrees)
        
        implicit none
        
        real, intent(IN):: degrees
        
        degrees_to_radians = degrees * pi / 180.

    end function degrees_to_radians



    subroutine init_rng(input_seed, fwd)
    ! initiate RNG state with reproducable state

        implicit none
        
        integer, optional, intent(IN) :: input_seed(:)
        logical, optional, intent(IN) :: fwd

        integer, allocatable :: seed(:)
        integer              :: n, i
        logical              :: ffwd
        real                 :: a


        call random_seed(size=n)
        allocate(seed(n))

        if(present(input_seed))then
            seed = 0
            seed = input_seed
        else
            seed = 1234567
        end if

        if(present(fwd))then
            ffwd = fwd
        else
            ffwd = .true.
        end if

        call random_seed(put=seed)

        !fast forward rng state 100 times to avoid any potential bad seeds
        if(ffwd)then
            call random_seed(get=seed)
            do i = 1, 100
                a = random()
                call random_seed(get=seed)
            end do
        end if
    end subroutine init_rng

    real function random_uni()

        implicit none

        call random_number(random_uni)
        
    end function random_uni

    real function random_uni_range(min, max)

        implicit none

        real, intent(IN) :: min, max

        random_uni_range = min + (max-min)*random_uni()

    end function random_uni_range

    real function clamp(x, min, max)

        implicit none

        real, intent(IN) :: x, min, max
        
        if(x < min)then
            clamp = min
            return
        end if
        if(x > max)then
            clamp = max
            return
        end if
        clamp = x

    end function clamp

end module utils