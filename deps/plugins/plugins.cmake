# open-plugins includes some private wxWidgets headers (wx/private/json.h, wx/gtk/private/*.h, ...)
# that a regular wxWidgets install does not ship: let it see the wxWidgets source include tree.
set(_op_wx_src_include "${CMAKE_CURRENT_BINARY_DIR}/dep_wxWidgets-prefix/src/dep_wxWidgets/include")
if (MSVC)
    set(_op_wx_private_include "/I${_op_wx_src_include}")
else ()
    set(_op_wx_private_include "-idirafter ${_op_wx_src_include}")
endif ()

AnycubicSlicer_add_cmake_project(
    OpenPlugins
    GIT_REPOSITORY "https://github.com/ANYCUBIC-3D/open-plugins.git"
    # Pinned to the open-plugins revision matching this source tree (early Feb 2026).
    # Later revisions require wxWidgets 3.3.
    GIT_TAG "f8199f36b2d4f749521bbe1b32d1c645d3b1cfa3"
    GIT_SHALLOW OFF
    DEPENDS ${WXWIDGETS_PKG} ${OPENSSL_PKG} dep_Boost ${CURL_PKG}
    CMAKE_ARGS
        -DENABLE_I18N=OFF
        "-DCMAKE_CXX_FLAGS_INIT:STRING=${_op_wx_private_include}"
)