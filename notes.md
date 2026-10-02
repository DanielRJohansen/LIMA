PS C:\Users\Daniel\git_repo\LIMA\distribution> .\release.bat
Current version in CMakeLists.txt: 1.3.0
Version to release (Enter keeps 1.3.0):
### LIMA 1.3.0 from commit c486a3240804c2a9149f70cd34c648461aabe265

### Building for Windows, CUDA architectures 89-real;90-real;100-real;120
--
 Resource files: C:/Users/Daniel/git_repo/LIMA/resources/logo/logoicon.rc
-- Configuring done (1.3s)
-- Generating done (0.1s)
-- Build files have been written to: C:/Users/Daniel/git_repo/LIMA/build/release-windows
[13/56] Building CXX object code\LIMA_BASE\CMakeFiles\LIMA_BASE.dir\src\Forcefield.cpp.obj
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\utility(274): warning C4267: 'initializing': conversion from 'size_t' to '_Ty2', possible loss of data
        with
        [
            _Ty2=int
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\utility(274): note: the template instantiation context (the oldest one first) is
C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\src\Forcefield.cpp(71): note: see reference to function template instantiation 'std::pair<const std::string,int>::pair<const std::string&,unsigned __int64,0>(_Other1,_Other2 &&) noexcept(false)' being compiled
        with
        [
            _Other1=const std::string &,
            _Other2=unsigned __int64
        ]
C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\src\Forcefield.cpp(71): note: see the first reference to 'std::pair<const std::string,int>::pair' in 'AtomtypeDatabase::GetActiveIndex'
[33/56] Building CXX object code\LIMA_MD\CMakeFiles\LIMA_MD.dir\src\Display.cpp.obj
C:\Users\Daniel\git_repo\LIMA\code\LIMA_MD\src\Display.cpp(546): warning C4305: 'return': truncation from 'int' to 'bool'
[37/56] Building CXX object code\LIMA_MD\CMakeFiles\LIMA_MD.dir\src\Environment.cpp.obj
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\optional(166): warning C4244: '=': conversion from 'const double' to '_Ty', possible loss of data
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\optional(166): note: the template instantiation context (the oldest one first) is
C:\Users\Daniel\git_repo\LIMA\code\LIMA_MD\src\Environment.cpp(772): note: see reference to function template instantiation 'std::optional<float> &std::optional<float>::operator =<const double&,0>(_Ty2) noexcept' being compiled
        with
        [
            _Ty2=const double &
        ]
C:\Users\Daniel\git_repo\LIMA\code\LIMA_MD\src\Environment.cpp(772): note: see the first reference to 'std::optional<float>::operator =' in 'Environment::UpdateSimstatus'
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\optional(312): note: see reference to function template instantiation 'void std::_Optional_construct_base<_Ty>::_Assign<const double&>(_Ty2) noexcept' being compiled
        with
        [
            _Ty=float,
            _Ty2=const double &
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(462): warning C4244: 'initializing': conversion from 'const double' to '_Ty', possible loss of data
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(462): note: the template instantiation context (the oldest one first) is
C:\Users\Daniel\git_repo\LIMA\code\LIMA_MD\src\Environment.cpp(751): note: see reference to function template instantiation 'float &std::vector<float,std::allocator<float>>::emplace_back<const double&>(const double &)' being compiled
C:\Users\Daniel\git_repo\LIMA\code\LIMA_MD\src\Environment.cpp(751): note: see the first reference to 'std::vector<float,std::allocator<float>>::emplace_back' in 'Environment::UpdateSimstatus'
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\vector(924): note: see reference to function template instantiation '_Ty &std::vector<_Ty,std::allocator<_Ty>>::_Emplace_one_at_back<const double&>(const double &)' being compiled
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\vector(845): note: see reference to function template instantiation '_Ty &std::vector<_Ty,std::allocator<_Ty>>::_Emplace_back_with_unused_capacity<const double&>(const double &)' being compiled
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\vector(860): note: see reference to function template instantiation 'void std::_Construct_in_place<float,const double&>(_Ty &,const double &) noexcept' being compiled
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(472): note: see reference to function template instantiation '_Ty *std::construct_at<_Ty,const double&>(_Ty *const ,const double &) noexcept(<expr>)' being compiled
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(476): warning C4244: 'initializing': conversion from 'const double' to 'float', possible loss of data
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(463): warning C4244: 'initializing': conversion from 'const double' to '_Ty', possible loss of data
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(463): note: the template instantiation context (the oldest one first) is
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\vector(860): note: see reference to function template instantiation 'void std::_Construct_in_place<float,const double&>(_Ty &,const double &) noexcept' being compiled
        with
        [
            _Ty=float
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\xutility(472): note: see reference to function template instantiation '_Ty *std::construct_at<_Ty,const double&>(_Ty *const ,const double &) noexcept' being compiled
        with
        [
            _Ty=float
        ]
[49/56] Building CUDA object code\LIMA_CONVEXHULLENGINE\CMakeFiles\LIMA_CONVEXHULLENGINE.dir\src\ConvexHullEngine.cu.ob
C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_CONVEXHULLENGINE\src\ConvexHullEngine.cu(197): warning #20054-D: dynamic initialization is not supported for a function-scope static __shared__ variable within a __device__/__global__ function
        __declspec(__shared__) TriList clippedFacets;
                                       ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_CONVEXHULLENGINE\src\ConvexHullEngine.cu(197): warning #20054-D: dynamic initialization is not supported for a function-scope static __shared__ variable within a __device__/__global__ function
        __declspec(__shared__) TriList clippedFacets;
                                       ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_CONVEXHULLENGINE\src\ConvexHullEngine.cu(197): warning #20054-D: dynamic initialization is not supported for a function-scope static __shared__ variable within a __device__/__global__ function
        __declspec(__shared__) TriList clippedFacets;
                                       ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec2.hpp(101): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec3.hpp(107): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_vec4.hpp(105): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x2.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat2x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x3.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __device__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\detail\type_vec1.hpp(95): warning #20012-D: __host__ annotation is ignored on a non-virtual function("vec") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr vec() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat3x4.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x2.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x3.hpp(36): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __device__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                           ^

C:\Users\Daniel\git_repo\LIMA\code\dependencies\glm\ext\../detail/type_mat4x4.hpp(35): warning #20012-D: __host__ annotation is ignored on a non-virtual function("mat") that is explicitly defaulted on its first declaration
                __declspec(__device__) __declspec(__host__) constexpr mat() = default ;
                                                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_CONVEXHULLENGINE\src\ConvexHullEngine.cu(197): warning #20054-D: dynamic initialization is not supported for a function-scope static __shared__ variable within a __device__/__global__ function
        __declspec(__shared__) TriList clippedFacets;
                                       ^

ConvexHullEngine.cu
tmpxft_000020dc_00000000-7_ConvexHullEngine.compute_120.cudafe1.cpp
[50/56] Building CUDA object code\LIMA_ENGINE\CMakeFiles\LIMA_ENGINE.dir\src\Engine.cu.obj
C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\include\BoxGrid.cuh(72): warning #69-D: integer conversion resulted in truncation
                        uint16_t blockId = 0xFFFFFFFF;
                                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                       ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                          ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(398): warning #68-D: integer conversion resulted in a change of sign
        const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                      ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EnergyMinimization.cuh(214): warning #68-D: integer conversion resulted in a change of sign
                const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\include\BoxGrid.cuh(72): warning #69-D: integer conversion resulted in truncation
                        uint16_t blockId = 0xFFFFFFFF;
                                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(22): warning #177-D: variable "byteSize" was declared but never referenced
                const size_t byteSize = sizeof(SuperCluster) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(SuperClusterMeta) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(int) * (nBlocks + 1) * 2;
                             ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                       ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                          ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(398): warning #68-D: integer conversion resulted in a change of sign
        const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                      ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EnergyMinimization.cuh(214): warning #68-D: integer conversion resulted in a change of sign
                const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\PME.cuh(500): warning #177-D: variable "X" was declared but never referenced
                                int X = ix - 1 + dx;
                                    ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\SuperclusterTaskBuilder.cuh(335): warning #177-D: variable "scIndexInSelf" was declared but never referenced
        const int scIndexInSelf = blockIdx.y;
                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\include\BoxGrid.cuh(72): warning #69-D: integer conversion resulted in truncation
                        uint16_t blockId = 0xFFFFFFFF;
                                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(22): warning #177-D: variable "byteSize" was declared but never referenced
                const size_t byteSize = sizeof(SuperCluster) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(SuperClusterMeta) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(int) * (nBlocks + 1) * 2;
                             ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                       ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                          ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(398): warning #68-D: integer conversion resulted in a change of sign
        const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                      ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EnergyMinimization.cuh(214): warning #68-D: integer conversion resulted in a change of sign
                const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\PME.cuh(500): warning #177-D: variable "X" was declared but never referenced
                                int X = ix - 1 + dx;
                                    ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\SuperclusterTaskBuilder.cuh(335): warning #177-D: variable "scIndexInSelf" was declared but never referenced
        const int scIndexInSelf = blockIdx.y;
                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\include\BoxGrid.cuh(72): warning #69-D: integer conversion resulted in truncation
                        uint16_t blockId = 0xFFFFFFFF;
                                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(22): warning #177-D: variable "byteSize" was declared but never referenced
                const size_t byteSize = sizeof(SuperCluster) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(SuperClusterMeta) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(int) * (nBlocks + 1) * 2;
                             ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                       ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                          ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(398): warning #68-D: integer conversion resulted in a change of sign
        const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                      ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EnergyMinimization.cuh(214): warning #68-D: integer conversion resulted in a change of sign
                const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\PME.cuh(500): warning #177-D: variable "X" was declared but never referenced
                                int X = ix - 1 + dx;
                                    ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\SuperclusterTaskBuilder.cuh(335): warning #177-D: variable "scIndexInSelf" was declared but never referenced
        const int scIndexInSelf = blockIdx.y;
                  ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_BASE\include\BoxGrid.cuh(72): warning #69-D: integer conversion resulted in truncation
                        uint16_t blockId = 0xFFFFFFFF;
                                           ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(22): warning #177-D: variable "byteSize" was declared but never referenced
                const size_t byteSize = sizeof(SuperCluster) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(SuperClusterMeta) * nBlocks * SuperClustersControl::maxClustersPerBlock + sizeof(int) * (nBlocks + 1) * 2;
                             ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                       ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                          ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\ParticleClusters.cuh(331): warning #221-D: floating-point value does not fit in required floating-point type
                meanPositionsOfPClusters[i] = Float3{ ((float)(1e+300)) ,((float)(1e+300)) , ((float)(1e+300)) };
                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(398): warning #68-D: integer conversion resulted in a change of sign
        const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                      ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EnergyMinimization.cuh(214): warning #68-D: integer conversion resulted in a change of sign
                const int scIdGlobal = (blockIdx.x * nScsPerBlock + threadIdx.y) < nSuperclusters ? blockIdx.x * nScsPerBlock + threadIdx.y : -1;
                                                                                                                                              ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\PME.cuh(500): warning #177-D: variable "X" was declared but never referenced
                                int X = ix - 1 + dx;
                                    ^

C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\SuperclusterTaskBuilder.cuh(335): warning #177-D: variable "scIndexInSelf" was declared but never referenced
        const int scIndexInSelf = blockIdx.y;
                  ^

Engine.cu
tmpxft_00006c70_00000000-7_Engine.compute_120.cudafe1.cpp
C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(17): warning C4083: expected ')'; found identifier 'E0020'
C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\EngineKernels.cuh(21): warning C4068: unknown pragma 'diag_suppress'
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\utility(296): warning C4244: 'initializing': conversion from 'const _Ty1' to '_Ty1', possible loss of data
        with
        [
            _Ty1=float
        ]
        and
        [
            _Ty1=int
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\utility(296): note: the template instantiation context (the oldest one first) is
C:\Users\Daniel\git_repo\LIMA\code\LIMA_ENGINE\src\Engine.cu(694): note: see reference to function template instantiation 'void std::sort<std::_Vector_iterator<std::_Vector_val<std::_Simple_types<_Ty>>>,Engine::TestAlgorithms::<lambda_1>::()::<lambda_1>>(const _RanIt,const _RanIt,_Pr)' being compiled
        with
        [
            _Ty=std::pair<float,int>,
            _RanIt=std::_Vector_iterator<std::_Vector_val<std::_Simple_types<std::pair<float,int>>>>,
            _Pr=Engine::TestAlgorithms::<lambda_1>::()::<lambda_1>
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\algorithm(8404): note: see reference to function template instantiation 'void std::_Sort_unchecked<std::pair<float,int>*,_Fn>(_RanIt,_RanIt,__int64,_Pr)' being compiled
        with
        [
            _Fn=Engine::TestAlgorithms::<lambda_1>::()::<lambda_1>,
            _RanIt=std::pair<float,int> *,
            _Pr=Engine::TestAlgorithms::<lambda_1>::()::<lambda_1>
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\algorithm(8374): note: see reference to function template instantiation '_BidIt std::_Insertion_sort_unchecked<_RanIt,_Pr>(const _BidIt,const _BidIt,_Pr)' being compiled
        with
        [
            _BidIt=std::pair<float,int> *,
            _RanIt=std::pair<float,int> *,
            _Pr=Engine::TestAlgorithms::<lambda_1>::()::<lambda_1>
        ]
C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.44.35207\include\algorithm(8249): note: see reference to function template instantiation 'std::pair<int,int>::pair<float,int,0>(const std::pair<float,int> &) noexcept' being compiled
[56/56] Linking CXX executable code\LIMA_TESTS\limaclitest.exe
Packaged C:\Users\Daniel\git_repo\LIMA\distribution\out\1.3.0\lima-1.3.0-windows-x64.zip (214 MB)

### Testing the packaged Windows build
Testing C:\Users\Daniel\git_repo\LIMA\distribution\out\1.3.0\lima-1.3.0-windows-x64\lima.exe
Output in C:\Users\Daniel\git_repo\LIMA\tests\clitests\_runs\20261002-204800
Leave the windows alone, they close by themselves
Test Every command has a test                                Success (1.90 [us])

=== lima dispatcher ===
Test lima dispatcher                                         12 commands respond to help, bad input exits 2 (2.31 [s])

=== lima makesimparams ===
> lima makesimparams

Test lima makesimparams                                      16 parameters (73.55 [ms])

=== lima makebox ===
> lima makebox --box-size 7 --name box

Test lima makebox                                            Empty 7 nm box (76.21 [ms])

=== lima editconf ===
> lima editconf -c met.gro -t met.top --conf-out out.gro --set-center 1.5 1.5 1.5 --rotate 0 0 1.5708

> lima render -f out.gro -t met.top

Test lima editconf                                           Rigid, centered, rotated (distance error 0.0009 nm) (4.82 [s])

=== lima togmx ===
> lima togmx -f 6lzm.pdb --name lzm

> lima render -f lzm.gro -t lzm.top

Test lima togmx                                              162 residues, 2610 atoms (4.88 [s])

=== lima solvate ===
> lima solvate -c met_box4.gro -t met_box4.top

> lima render -f met_box4_solvated.gro -t met_box4_solvated.top

Test lima solvate                                            2180 waters added (4.85 [s])

=== lima insertmolecule ===
> lima insertmolecule --conf-source met.gro --top-source met.top --conf-target box5.gro --top-target box5.top --position 2 2 2

> lima render -f box5.gro -t box5.top

Test lima insertmolecule                                     Inserted at the requested position (4.83 [s])

=== lima insertmolecules ===
> lima insertmolecules --conf-source met.gro --top-source met.top --conf-target met_box6.gro --top-target met_box6.top -n 10 --rotate-randomly -d
Task "FindIntersectIteration" - Calls: 1, Total: 1.336 ms, Average: 1.336 ms

Test lima insertmolecules                                    10 molecules, closest contact 0.49 nm (748.43 [ms])

=== lima em ===
> lima em -c metsol_clash.gro -t metsol.top --conf-out em.gro -d

Test lima em                                                 Clash resolved: 0.08 nm -> 0.28 nm (810.02 [ms])

=== lima mdrun ===
> lima mdrun -c metsol.gro -t metsol.top -s md_params.txt --conf-out out.gro --trajectory traj.trr --uff -d
Step #019999    Avg. time: 0.18msEngine time 3.586279

               Wall t (s)
       Time:         3.586
                 (ns/day)    (hour/ns)
Performance:       963.673          0.025

Test lima mdrun                                              2689 atoms, trajectory 3.2 MB (4.42 [s])

=== lima buildmembrane ===
> lima buildmembrane --lipids POPC 70 cholesterol 30 --box-size 8 --seed 1 -d
Step #001699    Avg. time: 0.39msbuildmembrane finished with a min max-force of 99.775

Test lima buildmembrane                                      216 lipids, 65% POPC, 0.59 nm^2/lipid (1.59 [s])

=== lima render ===
> lima render -f metsol.gro -t metsol.top --highlight 0 1 2

Test lima render
        Window after 0.4s, closed cleanly, screenshot render.bmp (5.57 [s])

=== lima selftest ===
> lima selftest
Selftest successful
Test lima selftest                                           Success (1.02 [s])


#--- Unittesting finished with 14 successes of 14 tests ---#


### Building for Linux in WSL
### Checking out c486a3240804c2a9149f70cd34c648461aabe265
### Building for CUDA architectures 89-real;90-real;100-real;120
-- The CXX compiler identification is GNU 14.2.0
-- The CUDA compiler identification is NVIDIA 13.2.86
-- The C compiler identification is GNU 14.2.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /usr/bin/g++-14 - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Detecting CUDA compiler ABI info
-- Detecting CUDA compiler ABI info - done
-- Check for working CUDA compiler: /usr/local/cuda-13.2/bin/nvcc - skipped
-- Detecting CUDA compile features
-- Detecting CUDA compile features - done
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /usr/bin/gcc-14 - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Found OpenGL: /usr/lib/x86_64-linux-gnu/libOpenGL.so
-- Found CUDAToolkit: /usr/local/cuda-13.2/targets/x86_64-linux/include (found version "13.2.86")
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE
--
 Resource files:
-- Configuring done (22.4s)
-- Generating done (0.0s)
-- Build files have been written to: /home/lima/lima-release/build
[104/104] Linking CXX executable code/LIMA/lima
### Packaging
Packaged /mnt/c/Users/Daniel/git_repo/LIMA/distribution/out/1.3.0/lima-1.3.0-linux-x86_64.tar.gz (90M)
### Smoke testing on the GPU
+ /home/lima/lima-release/lima-1.3.0-linux-x86_64/bin/lima --help
+ /home/lima/lima-release/lima-1.3.0-linux-x86_64/bin/lima makebox --box-size 5
+ /home/lima/lima-release/lima-1.3.0-linux-x86_64/bin/lima solvate -c met_box4.gro -t met_box4.top
+ /home/lima/lima-release/lima-1.3.0-linux-x86_64/bin/lima em -c metsol_clash.gro -t metsol.top --conf-out em.gro
+ printf 'n_steps=500\ndt=2\ndata_logging_interval=50\n'
+ /home/lima/lima-release/lima-1.3.0-linux-x86_64/bin/lima mdrun -c metsol.gro -t metsol.top -s params.txt --conf-out md.gro --trajectory md.trr
Engine time 0.114259

               Wall t (s)
       Time:         0.114
                 (ns/day)    (hour/ns)
Performance:       756.180          0.032
Smoke tests passed

### Debian package
Packaged /mnt/c/Users/Daniel/git_repo/LIMA/distribution/out/1.3.0/lima_1.3.0_amd64.deb (47M), depends: libc6 (>= 2.38), libglfw3 (>= 3.3), libglx0, libopengl0, libtbb12 (>= 2021.4.0)
### Test-installing in a clean Ubuntu container
debconf: delaying package configuration, since apt-utils is not installed
The .deb installs and runs

### Arch package
Wrote /mnt/c/Users/Daniel/git_repo/LIMA/distribution/out/1.3.0/PKGBUILD, depends: nvidia-utils glfw libglvnd onetbb glibc
### Test-building and installing in a clean Arch container
    lima-1.3.0-linux-x86_64.tar.gz ... Passed
The PKGBUILD builds, installs and runs

### Uploading draft release v1.3.0
HTTP 403: Resource not accessible by personal access token (https://api.github.com/repos/DanielRJohansen/LIMA/releases)

RELEASE FAILED: Uploading the release failed
PS C:\Users\Daniel\git_repo\LIMA\distribution>