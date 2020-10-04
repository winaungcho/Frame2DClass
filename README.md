# Frame2DClass
2 dimensional rigid jointed frame analysis class

 * Frame2DClass
 *
 * Class Inheritance
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
