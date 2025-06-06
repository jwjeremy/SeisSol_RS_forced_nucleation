# SPDX-FileCopyrightText: 2019 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

function(get_arch_flags architecture compiler)
    set(HAS_REDZONE OFF PARENT_SCOPE)

    # Westmere cpu architecture
    if ("${HOST_ARCH}" STREQUAL "wsm")
        set(HAS_REDZONE ON PARENT_SCOPE)
        set(CPU_ARCH_FLAGS "-msse3" PARENT_SCOPE)
    
    # Sandy Bridge cpu architecture
    elseif ("${HOST_ARCH}" STREQUAL "snb")
        set(HAS_REDZONE ON PARENT_SCOPE)
        set(CPU_ARCH_FLAGS "-mavx" PARENT_SCOPE)
    
    # Haswell cpu architecture
    elseif ("${HOST_ARCH}" STREQUAL "hsw")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX2" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-mavx2" "-mfma" PARENT_SCOPE)
        endif()

    # Knights Corner (Xeon Phi)
    elseif ("${HOST_ARCH}" STREQUAL "knc")
        set(HAS_REDZONE ON PARENT_SCOPE)
        set(CPU_ARCH_FLAGS "-mmic" "-fma" PARENT_SCOPE)
    
    # Knight Landing (Xeon Phi)
    elseif ("${HOST_ARCH}" STREQUAL "knl")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xMIC-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-mavx512f" "-mavx512cd" "-mavx512pf" "-mavx512er" "-mfma" PARENT_SCOPE)
        endif()
    
    # Skylake cpu architecture
    elseif ("${HOST_ARCH}" STREQUAL "skx")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=skylake-avx512" PARENT_SCOPE)
        endif()

    elseif ("${HOST_ARCH}" STREQUAL "thunderx2t99")
        if (compiler STREQUAL "Intel")

        elseif(compiler MATCHES "GNU|Clang")
	    # Note: mcpu/march/mtune are weird on arm, see:
	    # https://community.arm.com/developer/tools-software/tools/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
            set(CPU_ARCH_FLAGS "-mcpu=thunderx2t99" PARENT_SCOPE)
        endif()
    # AMD Zen 1 to 4
    elseif ("${HOST_ARCH}" STREQUAL "naples")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-march=core-avx2" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=znver1" "-mtune=znver1" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "rome")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-march=core-avx2" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=znver2" "-mtune=znver2" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "milan")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-march=core-avx2" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=znver3" "-mtune=znver3" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "bergamo")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=znver4" "-mtune=znver4" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "turin")
        set(HAS_REDZONE ON PARENT_SCOPE)
        if (compiler STREQUAL "Intel")
            set(CPU_ARCH_FLAGS "-xCORE-AVX512" "-fma" PARENT_SCOPE)
        elseif(compiler MATCHES "GNU|Clang|IntelLLVM")
            set(CPU_ARCH_FLAGS "-march=znver5" "-mtune=znver5" PARENT_SCOPE)
        endif()

    # IBM power 9
    elseif ("${HOST_ARCH}" STREQUAL "power9")
        if (compiler MATCHES "GNU|Clang")
            set(CPU_ARCH_FLAGS "-mtune=power9" PARENT_SCOPE)
        endif()

    elseif ("${HOST_ARCH}" STREQUAL "a64fx")
        if (compiler STREQUAL "FujitsuClang")
            set(CPU_ARCH_FLAGS "-KSVE " PARENT_SCOPE)
        else()
            set(CPU_ARCH_FLAGS "-mcpu=a64fx" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "neon")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "sve128")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "sve256")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "sve512")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "sve1024")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "sve2048")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "rvv128")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "rvv256")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "rvv512")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "rvv1024")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "rvv2048")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    elseif ("${HOST_ARCH}" STREQUAL "rvv4096")
        set(CPU_ARCH_FLAGS "" PARENT_SCOPE)
    
    elseif ("${HOST_ARCH}" STREQUAL "apple-m1")
        if (compiler MATCHES "GNU")
            set(CPU_ARCH_FLAGS "-march=armv8.4-a -mtune=generic" PARENT_SCOPE)
        endif()

    elseif ("${HOST_ARCH}" STREQUAL "apple-m2")
        if (compiler MATCHES "GNU")
            set(CPU_ARCH_FLAGS "-march=armv8.5-a" "-mtune=generic" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "apple-m3")
        if (compiler MATCHES "GNU")
            set(CPU_ARCH_FLAGS "-march=armv8.6-a" "-mtune=generic" PARENT_SCOPE)
        endif()
    elseif ("${HOST_ARCH}" STREQUAL "apple-m4")
        if (compiler MATCHES "GNU")
            set(CPU_ARCH_FLAGS "-march=armv9.2-a" "-mtune=generic" PARENT_SCOPE)
        endif()
    endif()

    if (compiler MATCHES "NVHPC|PGI")
        #NOTE: PGI-based compiler does not have `-mno-red-zone` flag
        set(HAS_REDZONE OFF PARENT_SCOPE)
    endif()

endfunction()
