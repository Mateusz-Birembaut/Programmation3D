HEADERS += \
    src/Camera.h \
    src/Exception.h \
    src/GLError.h \
    src/GLProgram.h \
    src/GLShader.h \
    src/Mesh.h \
    src/Scene.h \
    src/Trackball.h \
    src/Vec3.h \
    src/imageLoader.h

SOURCES += \
    gmini.cpp \
    src/Camera.cpp \
    src/GLError.cpp \
    src/GLProgram.cpp \
    src/GLShader.cpp \
    src/Mesh.cpp \
    src/Trackball.cpp

DISTFILES += \
    src/shader.frag \
    src/shader.vert
