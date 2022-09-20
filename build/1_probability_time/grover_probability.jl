using SparseArrays
using LinearAlgebra
using Random

L = 12;
Number_Of_Noise = 4*L^2-6*L+13;
Random.seed!(7000)
NOISE = 2*rand(Float64,Number_Of_Noise).-1;


Rx(theta) = exp(-1im*theta*[1 1;1 1]/2);
#Rx(theta) = [cos(theta/2) -1im*sin(theta/2) ; -1im*sin(theta/2)  cos(theta/2)];#

round.(-exp(-1im*pi*([1 1;1 1]/2)); digits = 3)

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
#Identity(3)

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


#=
The number of noise is total number of gates in the linear decomposition plus
the number of gates required to convert the MCX into a MCZ gate.
=#


A = ones(2^L,2^L);
U_x = (2/2^L)*A-Identity(2^L);

#=
The following function creates a multicontrolled X gate.

Input: DELTA (noise).
Output: Matrix of MCX.
=#

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
    MCX = Identity(2^L);
    
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
    return MCX
end    ;

# Total number of gates.
#print("Total number of gates = ",2*L^2-6*L+5)

#MCX[2^L,2^L-1],MCX[2^L-1,2^L]
#length(C_1)+length(C_2)+length(C_3)+length(C_4)+length(C_5)+length(C_6)

#=

MCZ = X^(1) X^(2)...X^(L-1) H^(t) MCX X^(1) X^(2)...X^(L-1) H^(t) = MCZ.
Creating a list for the gates on the left hand side of MCX gate.
=#
XHL_Gates = []
for i = 1:L-1
    push!(XHL_Gates,["X",i])
end    
push!(XHL_Gates,["H",L])

#=
Creating a list for the gates on the right hand side of MCX gate.
=#
XHR_Gates = [["H",L]]
for i = 1:L-1
    push!(XHR_Gates,["X",i])
end    

#=
The following function returns the matrix of U_0.
Input: Noise control parameter DELTA.
Output: Matrix of U_0.
=#

function U0_reconstructed(DELTA)
    
    Noise_Counter = 1

    # Creating an empty matrix to store the MCX matrix.
    MCX = Identity(2^L);
    
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


    #=
    Noise counter starts at the total number of gates required for
    the construction of the MCX value. 
    
    Total number noise created = Number of gates for MCX + Number of gate on the left +
                                Number of gate on right.
    =#
    
    XHL_Matrix = Identity(2^L)
    for i in XHL_Gates
        
        if i[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            XHL_Matrix = XHL_Matrix*Matrix_Gate(Hadamard(DELTA*epsilon), i[2]) 
            Noise_Counter += 1 
            
        elseif i[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            XHL_Matrix = XHL_Matrix*Matrix_Gate(1im*Rx(pi+DELTA*epsilon),i[2])
            Noise_Counter += 1 
            
        end
    end
    

    XHR_Matrix = Identity(2^L)
    for j in XHR_Gates
        if j[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            XHR_Matrix = XHR_Matrix*Matrix_Gate(Hadamard(DELTA*epsilon), j[2]) 
            Noise_Counter += 1 
            
        elseif j[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            XHR_Matrix = XHR_Matrix*Matrix_Gate(1im*Rx(pi+DELTA*epsilon),j[2])
            Noise_Counter += 1 
        end
    end
    #= MCZ = X^(1) X^(2)...X^(L-1) H^(t) MCX X^(1) X^(2)...X^(L-1) H^(t) = MCZ. =#
    return XHL_Matrix*MCX*XHR_Matrix
end;



function Ux_reconstructed(DELTA)

    Noise_Counter = 2*L^2-4*L+7

    # Creating an empty matrix to store the MCX matrix.
    MCX = Identity(2^L);
    
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




    
    #=
    Noise counter starts at the total number of gates required for
    the construction of the MCX value. 
    
    Total number noise created = Number of gates for MCX + Number of gate on the left +
                                Number of gate on right.
    =#
    
    
    HL_Matrix = Identity(2^L)
    for i in 1:L
        epsilon = NOISE[Noise_Counter]
        HL_Matrix = HL_Matrix*Matrix_Gate(Hadamard(DELTA*epsilon), i) 
        Noise_Counter += 1         
    end
    
    XHL_Matrix = Identity(2^L)
    for i in XHL_Gates
        
        if i[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            XHL_Matrix = XHL_Matrix*Matrix_Gate(Hadamard(DELTA*epsilon), i[2]) 
            Noise_Counter += 1 
            
        elseif i[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            XHL_Matrix = XHL_Matrix*Matrix_Gate(1im*Rx(pi+DELTA*epsilon),i[2])
            Noise_Counter += 1 
            
        end
    end
    
    XHR_Matrix = Identity(2^L)
    for j in XHR_Gates
        if j[1] == "H"
            
            epsilon = NOISE[Noise_Counter]
            XHR_Matrix = XHR_Matrix*Matrix_Gate(Hadamard(DELTA*epsilon), j[2]) 
            Noise_Counter += 1 
            
        elseif j[1] == "X"
            
            epsilon = NOISE[Noise_Counter]
            XHR_Matrix = XHR_Matrix*Matrix_Gate(1im*Rx(pi+DELTA*epsilon),j[2])
            Noise_Counter += 1 
        end
    end
    
    HR_Matrix = Identity(2^L)
    for i in 1:L
        epsilon = NOISE[Noise_Counter]
        HR_Matrix = HR_Matrix*Matrix_Gate(Hadamard(DELTA*epsilon), i) 
        Noise_Counter += 1         
    end
    
    #= MCZ = X^(1) X^(2)...X^(L-1) H^(t) MCX X^(1) X^(2)...X^(L-1) H^(t) = MCZ. =#
    return HL_Matrix*XHL_Matrix*MCX*XHR_Matrix*HR_Matrix
end;

Grover(DELTA) = Ux_reconstructed(DELTA) * U0_reconstructed(DELTA);



Normalized(Wavefunction) = Wavefunction/norm(Wavefunction);

function Pxbar(Wavefunction)
    s = sum(Wavefunction[2:length(Wavefunction)])
    return s*(conj.(s))/(2^L-1)
end

Psi_0(L) = (1/sqrt(2^L))*ones(ComplexF64,2^L);


py"""
f = open('probability_data'+'.txt', 'w')
def Write_file(p1, p2, i):
    f = open('probability_data'+'.txt', 'a')
    f.write(str(p1) +'\t'+ str(p2)+ '\t' + str(i) +'\n')
"""


psi = Psi_0(L)


Delta = parse(Float64,ARGS[1])

U =Grover(Delta)

for i=0:1000
    if i == 0
        p1 = psi[1]*conj.(psi[1])
        p2 = Pxbar(psi)
        py"Write_file"(real(p1),real(p2),i)
    else
        global psi = U*psi
        #global psi = Normalized(psi)
        p1 = psi[1]*conj.(psi[1])
        p2 = Pxbar(psi)
        py"Write_file"(real(p1),real(p2),i)
    end
end
        




