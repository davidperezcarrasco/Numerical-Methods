######################################################
# Constant values for the simulations
ArraySize = 129
ConstantD = 0.00000325
ConstantK = 0.05
WorldSize = 0.05
delta_t = 1.0 / 24.0
delta_x = WorldSize / ArraySize
T1 = 0.115
T2 = 0.0085
##
######################################################
######################################################
# Question PW1
##Inicializamos los State con los atributos a 0
StateI_J = [0.0 , [0.0]* ArraySize , [0]* ArraySize ]
StateI_GS =  [0.0 , [0.0]* ArraySize , [0]* ArraySize ]
StateCN_J =  [0.0 , [0.0]* ArraySize , [0]* ArraySize ]
StateCN_GS =  [0.0 , [0.0]* ArraySize , [0]* ArraySize ]
StateTime =  0
StateConcentration =  1
StateCellType =  2
##
######################################################
######################################################
# Visualization Parameters
# Amount of pixels per cell on the X axis
PixelsPerCell = 6
# Dimensions of the visualization window
WindowWidth = PixelsPerCell * ArraySize
WindowHeight = WindowWidth / 2
##
######################################################
######################################################
 
def CopyArray(inArray, outArray):
    assert(len(inArray) == len(outArray))
    for i in range(len(inArray)):
        outArray[i] = inArray[i]
        
def SetInitialState():
    ######################################################
    # Question P1
    #Inicializamos el tiempo a 0 en todos los state
    InitTime = 0.0
    StateI_J[StateTime] = InitTime
    EnforceBoundaryConditions(StateI_J[StateConcentration])
    
    StateI_GS[StateTime] = InitTime
    EnforceBoundaryConditions(StateI_GS[StateConcentration])
    
    StateCN_J[StateTime] = InitTime
    EnforceBoundaryConditions(StateCN_J[StateConcentration])
    
    StateCN_GS[StateTime] = InitTime
    EnforceBoundaryConditions(StateCN_GS[StateConcentration])
    

    ## End QUestion P1
    ######################################################

def EnforceBoundaryConditions(array):
    ######################################################
    # Question P2
    ##Aplicamos Neumann Boundary Conditions
    array[0] = array[1]
    array[-1] = array[-2]
    array[64] = 1
    ## End question P2
    ######################################################

def ComputeError(array_previous, array_next):
    ######################################################
    #### Question P5
    ## Calculate the error using the L2 distance
    sum = 0.0
    ux = [0.0]*ArraySize
    for i in range(1,ArraySize-1):
        ux[i] = array_next[i] - array_previous[i] #u_j^m+1 - u_j^m
        sum += sq(ux[i])
    L2 = sqrt(sum) #L2 norm
    return L2
    ## End Question P5
    ######################################################

def ComputeCellType(array):
    ##################################################
    #### Question P4  
    # Compute the cell type
    for i in range(ArraySize):
        array[StateCellType][i] = 0 #CellType inicial
        if array[StateConcentration][i] > T2:
            array[StateCellType][i] = 1 #tipo 1
        if array[StateConcentration][i] > T1:
            array[StateCellType][i] = 2 #tipo 2
    ## End Question P4
    ##################################################
       
def TimeStep(array, discretization, solver):
    # Set maximum iterations and tolerance threshold
    maxIterations = 20
    tolerance = 1e-10 # For Question P5
    
    ##################################################
    #### Question P3.1
    ## Initialize mu, aii (diagonal) and aij (values next to the diagonal)
    #describimos los distintos valores tal y como indicamos en la teoría
    if discretization == 'Implicit':
        mu = ConstantD*(delta_t/sq(delta_x))
        eta = ConstantK*delta_t
        aii = 1+2*mu+eta
        aij = -mu
    elif discretization == 'Crank-Nicolson':
        mu = ConstantD*(delta_t/(2*sq(delta_x)))
        eta = ConstantK*delta_t
        aii = 1+2*mu+eta
        aij = -mu
    else:
        raise ValueError("Discretization not implemented")
    #### End of question P3.1
    ##################################################
    
    
    ##################################################
    #### Question 3.2
    # Create an auxiliary vector u^{(m+1),(n)} of ArraySize, and initialize 
    # it by copying the array's contents into it##
    uCurrent = [0.0]*ArraySize
    CopyArray(array[StateConcentration],uCurrent) #copiamos el State actual al nuevo
    #### End of question P3.2
    ##################################################
    
    # Iterative Solver implementation:
    for iteration in range(maxIterations):
        # Empty the current values of u^{(m+1),(n+1)}
        uNext = [0.0]*ArraySize 

        ##############################################
        #### Question P3.3
        ## Update the boundary of u^{(m+1),(n+1)} with the values of u^{m}
        #Atribuímos los extremos del State anterior a los del nuevo
        uNext[0] = array[StateConcentration][0]
        uNext[-1] = array[StateConcentration][-1]
     
        #### End of question P3.3
        ##############################################
        
        # Iterative loop: apply Jacobi or Gauss-Seidel formula to calculate each cell of next x
        # Exclude the boundaries in the loop!
        for i in range(1,ArraySize-1):
            ##########################################
            #### Question P3.4
            ## Obtain the value of the summatory (T4) for the different discretizations.
            ## This value will be the same irrespective of the selected solver.   
            #Aplicamos las fórmulas de la teoría         
            if solver == 'Jacobi':
                summatory = aij*(uCurrent[i-1]+uCurrent[i+1])
            elif solver == 'Gauss-Seidel':
                summatory = aij*(uCurrent[i+1]+uNext[i-1])
            else:
                raise ValueError("Solver not implemented")
            #### End of question P3.4
            ##########################################

            ##########################################
            #### Question P3.5
            ## Obtain the value of the b term for the different solvers. This value will
            ## be the same irrespective of the discretization method       
            #Aplicamos las fórmulas de la teoría      
            if discretization == 'Implicit':
                bi = array[StateConcentration][i]
            elif discretization == 'Crank-Nicolson':
                bi = (1-2*mu)*array[StateConcentration][i]+mu*array[StateConcentration][i-1]+mu*array[StateConcentration][i+1]
            else:
                raise ValueError("Discretization not implemented")
            #### End of question P3.5
            ##########################################
                      
            ##########################################
            #### Question P3.6
            ## Implement Jacobi/Gauss-Seidel forumla. Given the above computations of 
            ## P3.4 and P3.5, the implementation should be the same for both cases     
            
            #Aplicamos las fórmulas de la teoría        
            uNext[i] = (bi-summatory)/aii

            #### End of question P3.6
            ##########################################
            
         ##############################################
        #### Question P3.7
        # Update the boundaries for u^{(m+1),(n+1)} using Neumann
        #Aplicamos las Boundary Conditions al nuevo State
        EnforceBoundaryConditions(uNext)
    
        #### Question P3.7
        #############################################
        # Calculate the error of u^{(m+1),(n+1)} with respect to u^{(m+1),(n)} 
        error = ComputeError(uCurrent,uNext)
        # If the error is less than the tolerance, then break the loop
        ## (For Q 2.5)
        if error <= tolerance: 
            print("Error < tolerance. Stopping at {} iterations...".format(iteration+1))
            break
        # Update u^{(m+1),(n)} with the obtained u^{(m+1),(n+1)}
        #Atribuímos el nuevo State al auxiliar
        CopyArray(uNext,uCurrent)


        
    ##################################################
    #### Question P3.8
    # Copy the result of the solver into the current state by:
    #Actualizamos el State
    CopyArray(uCurrent,array[StateConcentration])
    # Then, update the time
    array[StateTime] += delta_t
    ## End Question P3.8
    ##################################################
    
    ## Enforce boundary conditions
    EnforceBoundaryConditions(array[StateConcentration])
    ## Compute cell type
    ComputeCellType(array)
        
def DrawState():
    OffsetX = 50
    OffsetY = 0.8 * WindowHeight
    TextBoxSize = 40

    for i in range(ArraySize):
        PixelsX = i * (PixelsPerCell - 1)
        
        # Question PW6   
        ##################################################
        # Implicit_Jacobi: Fill with the name of the index
        SimY = StateI_J[StateConcentration][i]
        fill(SimY, 0.0, (1 - SimY))
        rect(OffsetX + PixelsX, 220, PixelsPerCell - 1, -150 * SimY)
        
        if StateI_J[StateCellType][i] == 2:
            fill(1, 0, 0)
        elif StateI_J[StateCellType][i] == 1:
            fill(1, 1, 1)
        elif StateI_J[StateCellType][i] == 0:
            fill(0, 0, 1)      
        for a in range(1,5):
            circle(OffsetX+PixelsX+PixelsPerCell/2, 250-a*PixelsPerCell-1, PixelsPerCell-1)
        
        # Implict-GS: Fill with the name of the index
        SimY = StateI_GS[StateConcentration][i]
        fill(SimY, 0.0, (1 - SimY))
        rect(OffsetX + PixelsX, 420, PixelsPerCell - 1, -150 * SimY)

        if StateI_GS[StateCellType][i] == 2:
            fill(1, 0, 0)
        elif StateI_GS[StateCellType][i] == 1:
            fill(1, 1, 1)
        elif StateI_GS[StateCellType][i] == 0:
            fill(0, 0, 1)      
        for a in range(1,5):
            circle(OffsetX + PixelsX+PixelsPerCell/2, 450-a*PixelsPerCell-1, PixelsPerCell-1)
        
        # CN-Jacobi: Fill with the name of the index
        SimY = StateCN_J[StateConcentration][i]
        fill(SimY, 0.0, (1 - SimY))
        rect(OffsetX + PixelsX, 620, PixelsPerCell - 1, -150 * SimY)

        if StateCN_J[StateCellType][i] == 2:
            fill(1, 0, 0)
        elif StateCN_J[StateCellType][i] == 1:
            fill(1, 1, 1)
        elif StateCN_J[StateCellType][i] == 0:
            fill(0, 0, 1)      
        for a in range(1,5):
            circle(OffsetX + PixelsX+PixelsPerCell/2, 650-a*PixelsPerCell-1, PixelsPerCell-1)

        # CN-Implicit: Fill with the name of the index
        SimY = StateCN_GS[StateConcentration][i]
        fill(SimY, 0.0, (1 - SimY))
        rect(OffsetX + PixelsX, 820, PixelsPerCell - 1, -150 * SimY)

        if StateCN_GS[StateCellType][i] == 2:
            fill(1, 0, 0)
        elif StateCN_GS[StateCellType][i] == 1:
            fill(1, 1, 1)
        elif StateCN_GS[StateCellType][i] == 0:
            fill(0, 0, 1)      
        for a in range(1,5):
            circle(OffsetX + PixelsX+PixelsPerCell/2, 850-a*PixelsPerCell-1, PixelsPerCell-1)

    # Protect the figure's name
    fill(1.0,1.0,1.0)
    rect(3.0,3.0,800-6,TextBoxSize)
    
def setup():
    SetInitialState()
    size(750, 860)
    
    WindowWidth = 750
    WindowHeight = 860

    colorMode(RGB, 1.0)
    noStroke()
    textSize(24)
    frameRate(1/delta_t)

def draw():
    background(0.9)

    TimeStep(StateI_J,'Implicit','Jacobi')
    TimeStep(StateI_GS,'Implicit','Gauss-Seidel')
    TimeStep(StateCN_J,'Crank-Nicolson','Jacobi')
    TimeStep(StateCN_GS,'Crank-Nicolson','Gauss-Seidel')
    DrawState()

    print("")

    # Label.
    fill(0.0)
    text("Morphogen diffusion model", 220, 32)
    text("Discretization:\nImplicit", 60, 100)
    text("Discretization:\nImplicit", 60, 300)
    text("Discretization:\nCrank-Nicolson", 60, 500)
    text("Discretization:\nCrank-Nicolson", 60, 700)
    text("Solver: Jacobi", 540, 100)
    text("Solver:\nGauss-Seidel", 540, 300)
    text("Solver: Jacobi", 540, 500)
    text("Solver:\nGauss-Seidel", 540, 700) 
