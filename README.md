# Frame2DClass
2 dimensional rigid jointed frame analysis class

## Frame2DClass
 
 ## Class Inheritance
 * [FEMSolver]
 *  -> [Frame2D]
 * [FEMSolver] is base class and includes several matrix operation for the standard FEM solutions.
 * [Frame2D] is a class for the FEM solution process and include data structure of 2 dimensional rigid jointed Frame.
 * [Frame2D] use [CSV] class for handling csv file to store and retrieve associated arrays;
 * Solution process run for the loaded truss to analyse deformations, reactions and element forces.
 * Multiple load cases will be solved simultaneously.
 * Html result tables are generated during the process.
 * Model of 2d Frame can be generated within class by assigning values to variables.
 * Or model can be created by loading CSV file.
 * Model csv file is very simple comma separated text file in-wich the properties of FEM element, boundary conditions and loads are written.

# Technical Reference

## Assumptions

1. Structure behave as a linear system.
2. Stress and strain inside of the nembers are small enough to be in the range of elastic portion.
3. Displacements of joints/nodes are small enough sothat secondary effects will be negalected.
4. Members are large enough to prevent bucklings.
5. Member end rotations are the same at the same node.

## FEM model

Each element/member connected to 2 joints/nodes.
Joint has 3 degrees of freedom, ux, uy and θz.
where
ux is translation along x axis.
uy is translation along y axis.
θz is rotation about z axis.

Element deform as a beam-bending and quadratic interpolation function between 2 nodes is used.
Pin-jointed truss can be idealized as a Pintruss3D.

## Applicable field
It can solve structural axial and bending elements problems such as moment resisting plane frame, continuous beams and gable frame.
Most of the applicable structures are made of reinforced concrete framed buildings.
