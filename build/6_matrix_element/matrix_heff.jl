using Random
using LinearAlgebra
using SparseArrays
using PyCall

L = 8;
Number_Of_Noise = 4*L^2-6*L+13-3;
#SEED = parse(Int64,ARGS[1])
Random.seed!(6000)
NOISE = 2*rand(Float64,Number_Of_Noise).-1;


#Rx(theta) = exp(-1im*(theta*[-1 1;1 -1]/2));
#Rx(theta) = [cos(theta/2) -1im*sin(theta/2) ; -1im*sin(theta/2)  cos(theta/2)];
Rx(theta) = exp(-1im*theta*[1 1;1 1]/2);
Rx(pi)

int(x) = floor(Int,x)

#Rx(theta) = exp(-1im*theta*[1 1;1 1]/2);
#Rx(theta) = [cos(theta/2) -1im*sin(theta/2) ; -1im*sin(theta/2)  cos(theta/2)];#

#round.(-exp(-1im*pi*([1 1;1 1]/2)); digits = 3)

Ry(theta) = [cos(theta/2) -sin(theta/2) ; sin(theta/2) cos(theta/2)];

Pauli_X = [0 1;1 0];
Pauli_Y = [1 -1im;1im 0];
Pauli_Z = [1 0;0 -1];

Hadamard(noise) = Ry(pi/2+noise)*Pauli_Z;

X = [0 1;1 0];

"""

Following function takes a 2x2 matrix (Gate) and qubit position (Qubit) and
returns the resultant matrix.

For example, the matrix for the gate U acting on the 3-rd qubit for N=5
qubit system is given by   I (x) I (x) U (x) I (x) I; where (x) is the
tensor product.

"""

function Matrix_Gate(Gate, Qubit) # Previously known as multi qubit gate.
    
    ## The case Qubit=1 is treated differently because we need to
    # initialize the matrix as U before starting the kronecker product.
    
    if Qubit == 1
        
        M = sparse(Gate)
        for i=2:L
            M = kron(M, sparse([1 0;0 1]))
        end
        
    else
        
        M = sparse([1 0;0 1])
        for i=2:L
            if i == Qubit
                M = kron(M, Gate)
            else
                M = kron(M, sparse([1 0;0 1]))
            end
        end
    end
    
    return M
end;

Identity(dimension) = 1* Matrix(I, dimension, dimension);

"""

The following function returns a controlled U gate matrix.

Input  : c (integer), t(integer), U (unitary operator).
Output : Matrix of the multicontrolled U gate with control qubit c and target qubit t.

"""

function CU(U,c,t)
    
    I2 = sparse([1 0;0 1])
    Z = sparse([1 0;0 -1])

    PI_0 = (I2+Z)/2
    PI_1 = (I2-Z)/2
     
    #function Rx(Noise)
        #A = cos((pi+Noise)/2)
        #B = -1im*sin((pi+Noise)/2)
        #return 1im*[A B;B A]
    #end
    
    Matrices = Dict("I" => I2,"PI_0" => PI_0,"U" => U, "PI_1" => PI_1)
    
    p0 = fill("I", L)
    p1 = fill("I", L)
    
    p0[c] = "PI_0"
    p1[c] = "PI_1"
    p1[t] = "U"

    
    PI_0_matrix = Matrices[p0[1]]
    for i = 2:L
        PI_0_matrix = kron(PI_0_matrix,Matrices[p0[i]])
    end        
        
    PI_1_matrix = Matrices[p1[1]]   
    for i = 2:L
        PI_1_matrix = kron(PI_1_matrix,Matrices[p1[i]])        
    end
           
    #return p0,p1
    return PI_0_matrix + PI_1_matrix     
end;

"""

The following returns a multicontrolled U gate matrix.

Input  : c (list), t(integer), U (unitary operator).
Output : Matrix of the multicontrolled U gate with control qubits c and target qubit t.

"""

function MCU(c,t,U)
    
    p0 = fill("I", L)
    p1 = fill("I", L)

    
    if typeof(c) == Int64
        p0[c] = "PI_1"
        p1[t] = "PI_1"
        
    else
        for i in c
            p0[i] = "PI_1"
            p1[i] = "PI_1"
        end
    end
    
    p0[t] = "I"
    p1[t] = "U"

    
    I = sparse([1 0;0 1])
    Z = sparse([1 0;0 -1])
    X = sparse([0 1;1 0])
    PI_0 = (I+Z)/2
    PI_1 = (I-Z)/2
     
    Matrices = Dict("I" => I,"PI_0" => PI_0,"U" => U, "PI_1" => PI_1)
    
    PI_0_matrix = Matrices[p0[1]]
    for i = 2:L
        PI_0_matrix = kron(PI_0_matrix,Matrices[p0[i]])
    end        
        
    PI_1_matrix = Matrices[p1[1]]   
    for i = 2:L
        PI_1_matrix = kron(PI_1_matrix,Matrices[p1[i]])        
    end
             
    # The identity in the following line needs to be replaced.
    return Identity(2^L) - PI_0_matrix + PI_1_matrix     
end;

function MCX_Reconstructed(DELTA)
    
    # Following iterates over the noise list.
    Noise_Counter = 1
    
    #C_1 = [];
    #C_2 = [];
    #C_3 = [];
    #C_4 = [];
    #C_5 = [];
    #C_6 = [];
    
    # Creating an empty matrix to store the MCX matrix.
    MCX = sparse(Identity(2^L));
    
    #=
    The following loops generates all the controlled Rx gates as
    described in arXiv:1303.3557. It generates six layers of gates
    as mentioned in the paper. The loops can be checked by running
    each of them manually.
    
    =#
    # C_1.
    for i = 1:L-2
        for j = 1:i
            #push!(C_1,[j,L-i,L-i+j])
            
            epsilon = NOISE[Noise_Counter]
            MCX = CU(Rx((pi/2^j)+DELTA*epsilon), L-i, L-i+j)*MCX
            Noise_Counter += 1
            
        end
    end

    # C_2.
    for i = 2:L
        #push!(C_2,[i-2,1,i])
        
        epsilon = NOISE[Noise_Counter]
        MCX = CU(Rx((pi/2^(i-2))+DELTA*epsilon), 1, i)*MCX
        Noise_Counter += 1
        
    end

    # C3 = - C1.
    for i = L-2:-1:1
        for j = i:-1:1
            #push!(C_3,[j,L-i,L-i+j])        
            epsilon = NOISE[Noise_Counter]
            MCX = CU(Rx((-pi/2^j)+DELTA*epsilon), L-i, L-i+j)*MCX
            Noise_Counter += 1
        end
    end

    # C_4.
    for i = 1:L-3
        for j = 1:i
            #push!(C_4,[j,L-i-1,L-i+j-1])          
            epsilon = NOISE[Noise_Counter]
            MCX = CU(Rx((pi/2^j)+DELTA*epsilon), L-i-1, L-i-1+j)*MCX
            Noise_Counter += 1
        end    
    end

    # C_5.
    for i = 2:L-1
        #push!(C_5,[i-2,1,i])       
        epsilon = NOISE[Noise_Counter]
        MCX = CU(Rx((-pi/2^(i-2))+DELTA*epsilon), 1, i)*MCX
        Noise_Counter += 1
        
    end

    # C6 = - C4.
    for i = L-3:-1:1
        for j = i:-1:1
            #push!(C_6,[j,L-i-1,L-i-1+j])         
            epsilon = NOISE[Noise_Counter]
            MCX = CU(Rx((-pi/2^j)+DELTA*epsilon), L-i-1, L-i-1+j)*MCX
            Noise_Counter += 1
            
        end    
    end
    return sparse(MCX)
end    ;

List_of_H = [];

#= U0 = MCZ = X^(1) X^(2)...X^(L-1) H^(t) MCX X^(1) X^(2)...X^(L-1) H^(t). =#

U0_XHL_Gates = []
for i = 1:L-1
    push!(U0_XHL_Gates,["X",i])
end    
push!(U0_XHL_Gates,["H",L])

#=
Creating a list for the gates on the right hand side of MCX gate.
=#
U0_XHR_Gates = []
for i = 1:L-1
    push!(U0_XHR_Gates,["X",i])
end
push!(U0_XHR_Gates,["H",L])

function U0_reconstructed(DELTA)
    
    # Following iterates over the noise list.
    Noise_Counter = 1
    
    #C_1 = [];
    #C_2 = [];
    #C_3 = [];
    #C_4 = [];
    #C_5 = [];
    #C_6 = [];
    
    # Creating an empty matrix to store the MCX matrix.
    MCX = sparse(Identity(2^L));
    
    #=
    The following loops generates all the controlled Rx gates as
    described in arXiv:1303.3557. It generates six layers of gates
    as mentioned in the paper. The loops can be checked by running
    each of them manually.
    
    =#
    # C_1.
    for i = 1:L-2
        for j = 1:i
            #push!(C_1,[j,L-i,L-i+j])
            
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((pi/2^j)+DELTA*epsilon), L-i, L-i+j)
            MCX = CRX*MCX
            push!(List_of_H,[pi/2^j+DELTA*epsilon, L-i, L-i+j])
            Noise_Counter += 1
            
        end
    end

    # C_2.
    for i = 2:L
        #push!(C_2,[i-2,1,i])
        
        epsilon = NOISE[Noise_Counter]
        CRX = CU(Rx((pi/2^(i-2))+DELTA*epsilon), 1, i)
        MCX = CRX*MCX
        push!(List_of_H, [pi/2^(i-2)+DELTA*epsilon, 1, i])
        Noise_Counter += 1
        
    end

    # C3 = - C1.
    for i = L-2:-1:1
        for j = i:-1:1
            #push!(C_3,[j,L-i,L-i+j])        
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((-pi/2^j)+DELTA*epsilon), L-i, L-i+j)
            MCX = CRX*MCX
            push!(List_of_H,[-pi/2^j+DELTA*epsilon, L-i, L-i+j])
            Noise_Counter += 1
        end
    end

    # C_4.
    for i = 1:L-3
        for j = 1:i
            #push!(C_4,[j,L-i-1,L-i+j-1])          
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((pi/2^j)+DELTA*epsilon), L-i-1, L-i-1+j)
            MCX = CRX*MCX
            push!(List_of_H,[pi/2^j+DELTA*epsilon,L-i-1, L-i-1+j])
            Noise_Counter += 1
        end    
    end

    # C_5.
    for i = 2:L-1
        #push!(C_5,[i-2,1,i])       
        epsilon = NOISE[Noise_Counter]
        CRX = CU(Rx((-pi/2^(i-2))+DELTA*epsilon), 1, i)
        MCX = CRX*MCX
        push!(List_of_H,[-pi/2^(i-2)+DELTA*epsilon,1, i])
        Noise_Counter += 1
        
    end

    # C6 = - C4.
    for i = L-3:-1:1
        for j = i:-1:1
            #push!(C_6,[j,L-i-1,L-i-1+j])         
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((-pi/2^j)+DELTA*epsilon), L-i-1, L-i-1+j)
            MCX = CRX*MCX
            push!(List_of_H,[-pi/2^j+DELTA*epsilon,L-i-1, L-i-1+j])
            Noise_Counter += 1
            
        end    
    end

    #=
    Noise counter starts at the total number of gates required for
    the construction of the MCX value. 
    
    Total number noise created = Number of gates for MCX + Number of gate on the left +
                                Number of gate on right.
    =#
    
    XHL_Matrix = sparse(Identity(2^L))
    for i in U0_XHL_Gates
        
        if i[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Hadamard(DELTA*epsilon), i[2]) 
            XHL_Matrix = XHL_Matrix*CRX
            push!(List_of_H,["H", DELTA*epsilon,i[2]])
            Noise_Counter += 1 
            
        elseif i[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Rx(pi+DELTA*epsilon),i[2])
            XHL_Matrix = XHL_Matrix*CRX
            push!(List_of_H,["X", DELTA*epsilon,i[2]])
            Noise_Counter += 1 
            
        end
    end
    

    XHR_Matrix = sparse(Identity(2^L))
    for j in U0_XHR_Gates
        if j[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Hadamard(DELTA*epsilon), j[2])
            XHR_Matrix = XHR_Matrix*CRX
            push!(List_of_H,["H", DELTA*epsilon,j[2]])
            Noise_Counter += 1 
            
        elseif j[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Rx(pi+DELTA*epsilon),j[2])
            XHR_Matrix = XHR_Matrix*CRX
            push!(List_of_H,["X", DELTA*epsilon,j[2]])
            Noise_Counter += 1 
        end
    end
    #= MCZ = X^(1) X^(2)...X^(L-1) H^(t) MCX X^(1) X^(2)...X^(L-1) H^(t) = MCZ. =#
    return sparse(XHL_Matrix*MCX*XHR_Matrix)
end;


#= Ux = H^(1) H^(2) ... H^(L-1) X^(1) X^(2)...X^(L-1) MCX X^(1) X^(2)...X^(L-1) H^(1) H^(2) ... H^(L-1). =#
Ux_XHL_Gates = []
for i = 1:L-1
    push!(Ux_XHL_Gates,["H",i])
end    
for i = 1:L-1
    push!(Ux_XHL_Gates,["X",i])
end  

Ux_XHR_Gates = []
for i = 1:L-1
    push!(Ux_XHR_Gates,["X",i])
end    
for i = 1:L-1
    push!(Ux_XHR_Gates,["H",i])
end

function Ux_reconstructed(DELTA)

    Noise_Counter = 2*L^2-4*L+7

    # Creating an empty matrix to store the MCX matrix.
    MCX = sparse(Identity(2^L));

    # C_1.
    for i = 1:L-2
        for j = 1:i
            #push!(C_1,[j,L-i,L-i+j])
            
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((pi/2^j)+DELTA*epsilon), L-i, L-i+j)
            MCX = CRX*MCX
            push!(List_of_H,[pi/2^j+DELTA*epsilon,L-i, L-i+j])
            Noise_Counter += 1
            
        end
    end

    # C_2.
    for i = 2:L
        #push!(C_2,[i-2,1,i])
        
        epsilon = NOISE[Noise_Counter]
        CRX = CU(Rx((pi/2^(i-2))+DELTA*epsilon), 1, i)
        MCX = CRX*MCX
        push!(List_of_H,[pi/2^(i-2)+DELTA*epsilon, 1, i])
        Noise_Counter += 1
        
    end

    # C3 = - C1.
    for i = L-2:-1:1
        for j = i:-1:1
            #push!(C_3,[j,L-i,L-i+j])        
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((-pi/2^j)+DELTA*epsilon), L-i, L-i+j)
            MCX = CRX*MCX
            push!(List_of_H,[-pi/2^j+DELTA*epsilon,L-i, L-i+j])
            Noise_Counter += 1
        end
    end

    # C_4.
    for i = 1:L-3
        for j = 1:i
            #push!(C_4,[j,L-i-1,L-i+j-1])          
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((pi/2^j)+DELTA*epsilon), L-i-1, L-i-1+j)
            MCX = CRX*MCX
            push!(List_of_H,[pi/2^j+DELTA*epsilon, L-i-1, L-i-1+j])
            Noise_Counter += 1
        end    
    end

    # C_5.
    for i = 2:L-1
        #push!(C_5,[i-2,1,i])       
        epsilon = NOISE[Noise_Counter]
        CRX = CU(Rx((-pi/2^(i-2))+DELTA*epsilon), 1, i)
        MCX = CRX*MCX
        push!(List_of_H,[-pi/2^(i-2)+DELTA*epsilon,1, i])
        Noise_Counter += 1
        
    end

    # C6 = - C4.
    for i = L-3:-1:1
        for j = i:-1:1
            #push!(C_6,[j,L-i-1,L-i-1+j])         
            epsilon = NOISE[Noise_Counter]
            CRX = CU(Rx((-pi/2^j)+DELTA*epsilon), L-i-1, L-i-1+j)
            MCX = CRX*MCX
            push!(List_of_H,[-pi/2^j+DELTA*epsilon,L-i-1, L-i-1+j])
            Noise_Counter += 1
            
        end    
    end


    XHL_Matrix = sparse(Identity(2^L))
    for i in Ux_XHL_Gates
        
        if i[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Hadamard(DELTA*epsilon), i[2])
            XHL_Matrix = XHL_Matrix*CRX
            push!(List_of_H,["H", DELTA*epsilon,i[2]])
            Noise_Counter += 1 
            
        elseif i[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Rx(pi+DELTA*epsilon),i[2])
            XHL_Matrix = XHL_Matrix*CRX
            push!(List_of_H,["X", DELTA*epsilon,i[2]])
            Noise_Counter += 1 
            
        end
    end
    
    XHR_Matrix = sparse(Identity(2^L))
    for j in Ux_XHR_Gates
        if j[1] == "H"          
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Hadamard(DELTA*epsilon), j[2]) 
            XHR_Matrix = XHR_Matrix*CRX
            push!(List_of_H,["H", DELTA*epsilon,j[2]])
            Noise_Counter += 1          
        elseif j[1] == "X"         
            epsilon = NOISE[Noise_Counter]
            CRX = Matrix_Gate(Rx(pi+DELTA*epsilon),j[2])
            XHR_Matrix = XHR_Matrix*CRX
            push!(List_of_H,["X", DELTA*epsilon,j[2]])
            Noise_Counter += 1 
        end
    end
    return sparse(XHL_Matrix*MCX*XHR_Matrix)
end;

Grover(DELTA) = collect(Ux_reconstructed(DELTA) * U0_reconstructed(DELTA));

U0 = collect(U0_reconstructed(0.0));
Ux = collect(Ux_reconstructed(0.0));
A = ones(2^L,2^L);
U_x = (2/2^L)*A-Identity(2^L);

function CRX_to_H(List)
    
    I2 = [[1,0] [0,1]]
    X = [[0,1] [1,0]]
    Z = [[1,0] [0,-1]]
    if List[1] == "H"
        Noise = List[2]
        Qubit = List[3]
        Unitary = Matrix_Gate(I2-Hadamard(Noise),Qubit)
        return Unitary
    elseif List[1] == "X"
        Noise = List[2]
        Qubit = List[3]
        Unitary = Matrix_Gate(Rx(pi/2+Noise)-I2,Qubit)
        return Unitary
    else #= The matrix is a CRX(theta+delta*epsilon).=#
        Angle = List[1]
        Control_Qubit = int(List[2])
        Target_Qubit = int(List[3])
        #= H = ((I-Z)/2)_c \otimes ((I+X)/2)_t.=#
        Matrices = Dict("I" => I2,"U" => (I2+X)/2, "PI_1" => (I2-Z)/2)
        p = fill("I", L)
        p[Control_Qubit] = "PI_1"
        p[Target_Qubit] = "U"
        Unitary = Matrices[p[1]]
        for i = 2:L
            Unitary = kron(Unitary,Matrices[p[i]])
        end         
        return Unitary  
    end
end;


HL = List_of_H;

function k_term(k)
        p = Identity(2^L);
        q = Identity(2^L);
        for i = k:length(HL)
            #= exp(-iH_{k+1} theta_{k+1})...exp(-iH_{N-1} theta_{N-1})exp(-iH_{N} theta_{N}). =#
            x1 = HL[length(HL)-i+k][1]#= theta/"H"/"X".=#
            x2 = HL[length(HL)-i+k][2]#= control_qubit/delta*epsilon/pi+delta*epsilon.=#
            x3 = HL[length(HL)-i+k][3]#= target_qubit/qubit/qubit.=#
            if x1 == "H"
                CRX = Matrix_Gate(Hadamard(x2),int(x3))
            elseif x1 == "X"
                CRX = Matrix_Gate(Rx(pi+x2),int(x3))
            else
                CRX = CU(Rx(x1), int(x2), int(x3))
            end
            p *= CRX # Counting in reverse.
              
            #= exp(iH_{k+1} theta_{k+1})...exp(iH_{N-1} theta_{N-1})exp(iH_{N} theta_{N}). =#
            y1 = HL[i][1]#= theta/"H"/"X".=#
            y2 = HL[i][2]#= control_qubit/delta*epsilon/pi+delta*epsilon.=#
            y3 = HL[i][3]#= target_qubit/qubit/qubit.=#
            #= The exponential factors has positive sign. We have to do some extra work. =#
            M1 = collect(CRX_to_H(HL[i]))
            
            if y1 == "H"
                CRX = exp(1im*y2*M1) #Matrix_Gate(Hadamard(x2),int(x3))
            elseif y1 == "X"
                CRX = exp(1im*y2*M1)#Matrix_Gate(Rx(pi+x2),int(x3))
            else
                CRX = exp(1im*y1*M1)#CU(Rx(x1), int(x2), int(x3))
            end            
            q *= CRX # Counting normally.
        end
    return p*CRX_to_H(HL[k])*q
end


function H_eff(DELTA)
    
    #= The following loop sums over all epsilon to get H_eff. =#
    h_eff = zeros(2^L,2^L);
    for i = 1:length(HL)
        h_eff += NOISE[i].*k_term(i)
    end
    return h_eff
end     

difference(DELTA) = Grover(0.0) - H_eff(DELTA);

my_delta = parse(Float64,ARGS[1])
D = H_eff(my_delta);

py"""
f = open('matrix_data'+'.txt', 'w')
def Write_file(row_index, column_index, element):
    f = open('matrix_data'+'.txt', 'a')
    f.write(str(row_index) +'\t'+ str(column_index)+ '\t' + str(element) +'\n')
"""

for i=1:2^L
    for j=1:2^L
        py"Write_file"(i,j,D[i,j])
    end
end
