cmake_minimum_required(VERSION 3.16)

project(
    affine_gap_aligner
    VERSION 1.0.0
    DESCRIPTION "Exact pairwise global alignment with affine gap penalties"
    LANGUAGES C
)

option(AGA_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)

add_library(affine_align_lib src/affine_align.c)
add_library(affine_align::lib ALIAS affine_align_lib)

target_include_directories(
    affine_align_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(affine_align_lib PUBLIC c_std_11)

add_executable(affine-align src/main.c)
target_link_libraries(affine-align PRIVATE affine_align::lib)

if(MSVC)
    target_compile_options(affine_align_lib PRIVATE /W4)
    target_compile_options(affine-align PRIVATE /W4)
else()
    set(AGA_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    target_compile_options(affine_align_lib PRIVATE ${AGA_WARNING_FLAGS})
    target_compile_options(affine-align PRIVATE ${AGA_WARNING_FLAGS})

    if(AGA_ENABLE_SANITIZERS)
        set(AGA_SANITIZER_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_compile_options(affine_align_lib PRIVATE ${AGA_SANITIZER_FLAGS})
        target_compile_options(affine-align PRIVATE ${AGA_SANITIZER_FLAGS})
        target_link_options(affine-align PRIVATE ${AGA_SANITIZER_FLAGS})
    endif()
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(affine-align-tests tests/test_affine_align.c)
    target_link_libraries(affine-align-tests PRIVATE affine_align::lib)

    if(MSVC)
        target_compile_options(affine-align-tests PRIVATE /W4)
    else()
        target_compile_options(affine-align-tests PRIVATE ${AGA_WARNING_FLAGS})
        if(AGA_ENABLE_SANITIZERS)
            target_compile_options(affine-align-tests PRIVATE ${AGA_SANITIZER_FLAGS})
            target_link_options(affine-align-tests PRIVATE ${AGA_SANITIZER_FLAGS})
        endif()
    endif()

    add_test(NAME library-tests COMMAND affine-align-tests)
    add_test(
        NAME cli-direct-input
        COMMAND affine-align
                --seq-a ACGT
                --seq-b AGT
                --match 2
                --mismatch 1
                --gap-open 3
                --gap-extend 1
                --format json
    )
    set_tests_properties(
        cli-direct-input
        PROPERTIES PASS_REGULAR_EXPRESSION "\"score\": 3"
    )

    add_test(
        NAME cli-fasta-input
        COMMAND affine-align --fasta ${CMAKE_CURRENT_SOURCE_DIR}/examples/pair.fasta
    )
    set_tests_properties(
        cli-fasta-input
        PROPERTIES PASS_REGULAR_EXPRESSION "Affine Global Alignment"
    )
endif()

install(TARGETS affine_align_lib affine-align)
install(FILES include/affine_align.h DESTINATION include)
