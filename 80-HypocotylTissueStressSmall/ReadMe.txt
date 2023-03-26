Simulation of cracking in the hypocotyl.

OutsideFaces.txt - outside faces, used to change growth and threshold
AllFaces.txt - used to set FEM parameters
CrossWalls.txt - used to set FEM parameters for crosswalls
HorizontalGrowth.txt - two vertices used to measure transverse strain/growth
VerticalGrowth.txt - two vertices used to measure longitudinal strain/growth
CenterVertex.txt - center vertex for Dirichlet (1,1,1)
Dirichlet.txt - vertices in XZ plane for Dirichlet (0,1,0)

Crack2.txt - vertices to fracture for 2 cell crack
Crack4.txt - vertices to fracture for 4 cell crack
Crack8.txt - vertices to fracture for 8 cell crack

CrackedCells2.txt - 3D cells for 2 cell crack
CrackedCells4.txt - 3D cells for 4 cell crack
CrackedCells8.txt - 3D cells for 8 cell crack

Shared parameters
Young's modulus longitudinal 150
Young's modulus transverse 1050
Crosswalls 1050:1050
Poisson's ratio .3
Growth Time: 25

WT/Qua
Inner wall .2um
Outer wall 2um
Inner growth 1
Outer growth 10
Inner strain thresh .02
Outer strain thresh .01

Dwarf
Inner wall .2um
Outer wall 2um
Inner growth 1
Outer growth 1 
Inner strain thresh .02
Outer strain thresh .02

Uniform
Inner wall .2um
Outer wall .2um
Inner growth 1
Outer growth 1 
Inner strain thresh .02
Outer strain thresh .02

Lengths pressurized
approx 65.5 in height
approx 113 in width

Single Cell Simulation for Tissue Stress/Specified growth:

Uniform - Cell:315 has cross section area:339.684, wall area:13.2203  
2 Walls - Cell:315 has cross section area:334.018, wall area:51.8848
2.5 Walls - Cell:315 has cross section area:332.874, wall area:61.5537
3 Walls  - Cell:315 has cross section area:331.799, wall area:71.231   
All walls stiff - Cell:315 has cross section area:324.623, wall area:129.24

Used an in-between amount for inner/outer factors (2.5 walls), inner 25.69, outer 5.4, Cell Type: 1 inner, 2 outer


To make video with text:
ffmpeg -i Frame-%5d.jpg -vf "crop=648:1226,drawtext=fontfile=../ArialItalic.ttf:text='qua2-1':fontcolor=white:fontsize=48:box=1:boxcolor=black@0.5:boxborderw=5:x=(w-text_w)/2:y=(7.1*h)/8" qua2-1.mp4

Frame size: 648x1226

Videos
After simulation run 45 steps
Turn clipping off run 30 steps
Do 400 degree rotation over 120 steps
Hold for 60 steps
