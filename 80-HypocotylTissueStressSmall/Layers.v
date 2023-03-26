// Make a block of 3D plant tissue
// Coordinate system (x,y,z) or ...
//                   (a,r,l) (angle, radius, length) for cylindrical

[Main]
Debug: 1
Epsilon: .00001
SegSz: 4 4 4 // size of segments

// The next section defines the cell layers, Go Wild!
// These are the cell sizes in segments
LayerSegs: 5 5 21
LayerSegs: 5 5 21
LayerSegs: 5 5 21
LayerSegs: 5 5 21
LayerSegs: 5 5 21

// Starting points of cells in segments
LayerStart: 0 0 0
LayerStart: 0 5 0
LayerStart: 0 10 0
LayerStart: 0 15 0
LayerStart: 0 20 0

// Ending points of cells in segments
LayerEnd: 5 5 126
LayerEnd: 5 10 126
LayerEnd: 5 15 126
LayerEnd: 5 20 126
LayerEnd: 5 25 126

// Stagger cells
Stagger: 0
Stagger: 6
Stagger: 0
Stagger: 6
Stagger: 0

// Make layers or not (good for debugging)
LayerMake: true
LayerMake: true
LayerMake: true
LayerMake: true
LayerMake: true
