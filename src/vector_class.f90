Module vec3_class

    type :: vec3
        real :: x, y, z
        contains

        procedure :: magnitude       => magnitude_fn
        procedure :: length          => length
        procedure :: length_squared  => length_squared
        procedure :: near_zero       => near_zero

        generic   :: operator(.dot.)   => vec_dot_vec, vec_dot_mat
        generic   :: operator(.cross.) => vec_cross_vec
        generic   :: operator(/)       => vec_div_scal
        generic   :: operator(*)       => vec_mult_vec, vec_mult_scal, scal_mult_vec
        generic   :: operator(+)       => vec_add_vec, vec_add_scal, scal_add_vec
        generic   :: operator(-)       => vec_minus_vec, vec_minus_scal, scal_minus_vec

        procedure, pass(a), private :: vec_dot_vec
        procedure, pass(a), private :: vec_dot_mat

        procedure, pass(a), private :: vec_cross_vec

        procedure, pass(a), private :: vec_div_scal

        procedure, pass(a), private :: vec_mult_vec
        procedure, pass(a), private :: vec_mult_scal
        procedure, pass(b), private :: scal_mult_vec

        procedure, pass(a), private :: vec_add_vec
        procedure, pass(a), private :: vec_add_scal
        procedure, pass(b), private :: scal_add_vec

        procedure, pass(a), private :: vec_minus_vec
        procedure, pass(a), private :: vec_minus_scal
        procedure, pass(b), private :: scal_minus_vec

    end type vec3

    private
    public :: magnitude, vec3, abs, length, max, invert, near_zero, reflect, refract
    public :: random_in_unit_sphere, random_vec3, random_unit_vector,random_in_unit_hemisphere, random_in_unit_disk


    interface abs
        module procedure abs_vec
    end interface abs

    interface max
        module procedure max_vec
    end interface max

    interface random_vec3
        module procedure random_vec3_uni
        module procedure random_vec3_uni_range
    end interface random_vec3

    contains

        type(vec3) function refract(uv, n, etai_over_etat)
            
            implicit none
            
            type(vec3), intent(IN) :: uv, n
            real,       intent(IN) :: etai_over_etat
            
            type(vec3) :: r_out_perp, r_out_parallel
            real       :: cos_theta

            cos_theta = min((-1.)*uv .dot. n, 1.)
            r_out_perp = etai_over_etat * (uv + cos_theta*n)
            r_out_parallel = (-1.)*sqrt(abs(1. - r_out_perp%length_squared())) * n

            refract = r_out_perp + r_out_parallel

        end function refract


        type(vec3) function reflect(v, n)
            
            implicit none

            type(vec3), intent(IN) :: v, n
                
            reflect = v - 2.* (v .dot. n) * n

        end function reflect

        logical function near_zero(this)

            implicit none

            class(vec3) :: this
            real :: s

            s = 1e-8
            near_zero = .false.
            if((abs(this%x) < s) .and. (abs(this%y) < s) .and. (abs(this%z) < s))near_zero=.true.
            
        end function near_zero

        type(vec3) function random_vec3_uni()
                
            use utils, only : random

            implicit none
            
            random_vec3_uni = vec3(random(), random(), random())

        end function random_vec3_uni


        type(vec3) function random_vec3_uni_range(min, max)
                
            use utils, only : random

            implicit none

            real, intent(IN) :: min, max
            
            random_vec3_uni_range = vec3(random(min, max), random(min, max), random(min, max))

        end function random_vec3_uni_range

        type(vec3) function random_in_unit_sphere()

            implicit none

            type(vec3) :: p

            do
                p = random_vec3(-1., 1.)
                if(p%length_squared() >= 1)continue
                exit
            end do
            
            random_in_unit_sphere = p

        end function random_in_unit_sphere

        type(vec3) function random_unit_vector()

            implicit none
            
            random_unit_vector = random_in_unit_sphere()
            random_unit_vector = random_unit_vector%magnitude()

        end function random_unit_vector


        type(vec3) function random_in_unit_hemisphere(normal)

            implicit none
            
            type(vec3), intent(IN) :: normal

            type(vec3) :: in_unit_sphere

            in_unit_sphere = random_in_unit_sphere()
            if((in_unit_sphere .dot. normal) > 0.)then
                random_in_unit_hemisphere = in_unit_sphere
            else
                random_in_unit_hemisphere = (-1.) * in_unit_sphere
            end if

        end function random_in_unit_hemisphere

        type(vec3) function random_in_unit_disk()

            use utils, only : random

            implicit none

            type(vec3) :: p

            do
                p = vec3(random(-1., 1.), random(-1., 1.), 0.)
                if(p%length_squared() >= 1.)continue
                exit
            end do  

            random_in_unit_disk = p          
        end function random_in_unit_disk


        type(vec3) function abs_vec(this)

            implicit none

            type(vec3), intent(IN) :: this

            abs_vec = vec3(abs(this%x), abs(this%y), abs(this%z))

        end function abs_vec

        type(vec3) function max_vec(this, val)

            implicit none

            type(vec3), intent(IN) :: this
            real, intent(IN) :: val

            max_vec = vec3(max(this%x, val), max(this%y, val), max(this%z, val))

        end function max_vec

        type(vec3) function min_fn(this, val)

            implicit none

            type(vec3), intent(IN) :: this
            real, intent(IN) :: val

            min_fn = vec3(min(this%x, val), min(this%y, val), min(this%z, val))

        end function min_fn

        type(vec3) function vec_minus_vec(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            type(vec3),  intent(IN) :: b

            vec_minus_vec = vec3(a%x - b%x, a%y - b%y, a%z - b%z)

        end function vec_minus_vec


        type(vec3) function vec_add_scal(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            real,        intent(IN) :: b

            vec_add_scal = vec3(a%x + b, a%y + b, a%z + b)

        end function vec_add_scal


        type(vec3) function scal_add_vec(a, b)

            implicit none

            class(vec3), intent(IN) :: b
            real,          intent(IN) :: a

            scal_add_vec = vec3(b%x + a, b%y + a, b%z + a)

        end function scal_add_vec


        type(vec3) function vec_minus_scal(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            real,          intent(IN) :: b

            vec_minus_scal = vec3(a%x - b, a%y - b, a%z - b)

        end function vec_minus_scal


        type(vec3) function scal_minus_vec(a, b)

            implicit none

            class(vec3), intent(IN) :: b
            real,          intent(IN) :: a

            scal_minus_vec = vec3(b%x - a, b%y - a, b%z - a)

        end function scal_minus_vec


        type(vec3) function vec_add_vec(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            type(vec3),  intent(IN) :: b

            vec_add_vec = vec3(a%x + b%x, a%y + b%y, a%z + b%z)

        end function vec_add_vec


        elemental function vec_dot_vec(a, b) result (dot)

            implicit none

            class(vec3), intent(IN) :: a
            type(vec3),  intent(IN) :: b
            real :: dot

            dot = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)

        end function vec_dot_vec

        elemental function vec_cross_vec(a, b) result (cross)

            implicit none

            class(vec3), intent(IN) :: a
            type(vec3),  intent(IN) :: b
            type(vec3) :: cross

            cross = vec3(a%y * b%z - a%z * b%y, &
                         a%z * b%x - a%x * b%z, &
                         a%x * b%y - a%y * b%x)

        end function vec_cross_vec

        function vec_dot_mat(a, b) result (dot)

            implicit none

            class(vec3), intent(IN) :: a
            real,          intent(IN) :: b(4, 4)
            type(vec3) :: dot

            dot%x = b(1, 1)*a%x + b(2, 1)*a%y + b(3, 1)*a%z + b(4, 1)*1.
            dot%y = b(1, 2)*a%x + b(2, 2)*a%y + b(3, 2)*a%z + b(4, 2)*1.
            dot%z = b(1, 3)*a%x + b(2, 3)*a%y + b(3, 3)*a%z + b(4, 3)*1.

        end function vec_dot_mat

        type(vec3) function vec_mult_vec(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            type(vec3),  intent(IN) :: b

            vec_mult_vec = vec3(a%x * b%x, a%y * b%y, a%z * b%z)

        end function vec_mult_vec


        type(vec3) function vec_mult_scal(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            real,          intent(IN) :: b

            vec_mult_scal = vec3(a%x * b, a%y * b, a%z * b)

        end function vec_mult_scal


        type(vec3) function scal_mult_vec(a, b)

            implicit none

            class(vec3), intent(IN) :: b
            real,          intent(IN) :: a

            scal_mult_vec = vec3(a * b%x, a * b%y, a * b%z)

        end function scal_mult_vec


        type(vec3) function vec_div_scal(a, b)

            implicit none

            class(vec3), intent(IN) :: a
            real,         intent(IN) :: b

            vec_div_scal = vec3(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal


        type(vec3) function magnitude_fn(this)

            implicit none

            class(vec3) :: this

            real :: tmp

            tmp = sqrt(this%x**2 + this%y**2 + this%z**2)
            magnitude_fn = this / tmp

        end function magnitude_fn


        real function length(this)

            implicit none

            class(vec3) :: this

            length = sqrt(this%x**2 + this%y**2 + this%z**2)

        end function length

        real function length_squared(this)

            implicit none

            class(vec3) :: this

            length_squared = this%x**2 + this%y**2 + this%z**2

        end function length_squared

    pure function invert(A) result(B)
        !! from http://fortranwiki.org/fortran/show/Matrix+inversion
        !! Performs a direct calculation of the inverse of a 4×4 matrix.
        real, intent(in) :: A(4,4)   !! Matrix

        real             :: B(4,4)   !! Inverse matrix
        real             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1./(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
        - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
        + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
        - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end function invert
end Module vec3_class