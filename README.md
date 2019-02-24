# MathPlotter
A function plotter written in C++ with only OpenGl / ImGui / mathExp libraries

### Guide to Use
Add to **compiler** options
-DGLEW_STATIC -DNO_FREETYPE

Add to **linker** options:
-static-libgcc -lglew32 -lglfw3 -lopengl32 -lglu32 -lgdi32 -lwinmm

Or Download the "Math.exe" file
