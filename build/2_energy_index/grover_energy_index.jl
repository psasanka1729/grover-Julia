using SparseArrays
using LinearAlgebra
using Random
using PyCall
L = 12;
Number_Of_Noise = 4*L^2-6*L+13;
SEED = parse(Int64,ARGS[1]);
Random.seed!(SEED)
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


#Grover(0.1);

sigma_x = [0 1;1 0];
sigma_y = [0 -1im;1im 0];
sigma_z = [1 0;0 -1];
sigma_plus = (sigma_x+1im*sigma_y)/2;
sigma_minus = (sigma_x-1im*sigma_y)/2;

adjoint(psi) = psi'; # conjugate transpose.

psi_to_rho(psi) = psi'*psi;

exp_val(psi, op) = adjoint(psi)*real(op*psi);

norm_sq(psi) = real(adjoint(psi)*psi);

function normalize(psi, tol=1.e-9)
    ns = norm_sq(psi)^0.5 
    return psi/ns
end;
 
function is_herm(M,tol=1.e-9)
    
    if size(M)[1] != size(M)[2]
        return "False"
    end
    diff = M-adjoint(M)
    
    return max(abs(collect(Iterators.flatten([diff])))) < tol
end;

function is_unitary(M,tol=1.e-9)
    
    if size(M)[1] != size(M)[2]
        return "False"
    end
    diff = M*adjoint(M)-Identity(size(M)[1])
    return max(abs(collect(Iterators.flatten([diff])))) < tol
end

function eigu(U,tol = 1.e-9)
    
    (E_1,V_1) = (eigvals(U+adjoint(U)), eigvecs(U+adjoint(U)))
    U_1 = adjoint(V_1)*U*V_1
    H_1 = adjoint(V_1)*(U+adjoint(U))*V_1
    non_diag_lst = []
    j = 1
    while j < size(U_1)[1]
        k = 1
        while k < size(U_1)[1]
            if j!=k && abs(U_1[j,k]) > tol
                if j ∉ non_diag_lst
                    push!(non_diag_lst,j)
                end
                if k ∉ non_diag_lst
                    push!(non_diag_lst,k)
                end
                  k += 1
            end
            j += 1
        end
    end
    
    if length(non_diag_lst) > 0
        non_diag_lst = sort(non_diag_lst)
        U_1_cut = U_1[non_diag_lst,:][:,non_diag_lst]
        (E_2_cut, V_2_cut) = (eigvals(1im*(U_1_cut-adjoint(U_1_cut))), eigvecs(1im*(U_1_cut-adjoint(U_1_cut))))
        V_2 = Identity(size(U)[1])
        for j=1:length(non_diag_lst)
            V_2[non_diag_lst[j],non_diag_lst] = V_2_cut[j,:]
        end
        V_1 = V_1*V_2
        U_1 = adjoint(V_2)*adjoint(U_1)*adjoint(V_2)
    end
    
    U_1 = diag(U_1)
    inds= sortperm(imag(log(U_1)))
    return U_1[inds],V_1[:inds]
end         
#eigu(Grover(0.0))

using PyCall
py"""
import numpy
import numpy.linalg
def adjoint(psi):
    return psi.conjugate().transpose()
def psi_to_rho(psi):
    return numpy.outer(psi,psi.conjugate())
def exp_val(psi, op):
    return numpy.real(numpy.dot(adjoint(psi),op.dot(psi)))
def norm_sq(psi):
    return numpy.real(numpy.dot(adjoint(psi),psi))
def normalize(psi,tol=1e-9):
    ns=norm_sq(psi)**0.5
    if ns < tol:
        raise ValueError
    return psi/ns
def is_herm(M,tol=1e-9):
    if M.shape[0]!=M.shape[1]:
        return False
    diff=M-adjoint(M)
    return max(numpy.abs(diff.flatten())) < tol
def is_unitary(M,tol=1e-9):
    if M.shape[0]!=M.shape[1]:
        return False
    diff=M.dot(adjoint(M))-numpy.identity((M.shape[0]))
    return max(numpy.abs(diff.flatten())) < tol
def eigu(U,tol=1e-9):
    (E_1,V_1)=numpy.linalg.eigh(U+adjoint(U))
    U_1=adjoint(V_1).dot(U).dot(V_1)
    H_1=adjoint(V_1).dot(U+adjoint(U)).dot(V_1)
    non_diag_lst=[]
    j=0
    while j < U_1.shape[0]:
        k=0
        while k < U_1.shape[0]:
            if j!=k and abs(U_1[j,k]) > tol:
                if j not in non_diag_lst:
                    non_diag_lst.append(j)
                if k not in non_diag_lst:
                    non_diag_lst.append(k)
            k+=1
        j+=1
    if len(non_diag_lst) > 0:
        non_diag_lst=numpy.sort(numpy.array(non_diag_lst))
        U_1_cut=U_1[non_diag_lst,:][:,non_diag_lst]
        (E_2_cut,V_2_cut)=numpy.linalg.eigh(1.j*(U_1_cut-adjoint(U_1_cut)))
        V_2=numpy.identity((U.shape[0]),dtype=V_2_cut.dtype)
        for j in range(len(non_diag_lst)):
            V_2[non_diag_lst[j],non_diag_lst]=V_2_cut[j,:]
        V_1=V_1.dot(V_2)
        U_1=adjoint(V_2).dot(U_1).dot(V_2)
    # Sort by phase
    U_1=numpy.diag(U_1)
    inds=numpy.argsort(numpy.imag(numpy.log(U_1)))
    return (U_1[inds],V_1[:,inds]) # = (U_d,V) s.t. U=V*U_d*V^\dagger
"""
#G = py"eigu"(Grover(0.0));

#real(1im*log.(G[1]))

#Gr = py"eigu"(Grover(0.0))[1]
#real(1im*log.(Gr[1]))

function Entropy(Psi)   
    
    LS = Int64(L/2)

    #Psi = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1im];
    
    # Normalizing Psi.
    Psi = Psi/norm(Psi) 
    
    psi(s) = Psi[(2^LS)*(s-1)+1:(2^LS)*s]
    
    #=
        psi(s_p) is a row matrix/vector. psi(s) is a column matrix/vector.      
        Dimension of rhoA is N/2 x N/2. 
        The element <s|rhoA|sp> is given by psi_sp^\dagger * psi_s.
    =#
    
    
    # psi(s_p)^\dagger * psi(s) is the element of (s,s_p) of rho_AB. 
    rhoA(s,s_p) = psi(s_p)' * psi(s)
    
    
    # Following function returns the reduced density matrix rho_A.
    function rhoA_Matrix()
        
        LS = Int64(L/2)
            
        # Creates a zero matrix to store the density matrix.
        M = zeros(ComplexF64,2^LS,2^LS)
        
        #=
        rho is Hermitian, it is sufficient to calculate the elements above the diagonal.
        The the elements below the diagonal can be replace by the complex cpnjugate of the
        elements above the diagonal.
        =#
    
        for i=1:2^LS
            for j=1:2^LS
                if i<=j
                    M[i,j] = rhoA(i,j)
                else
                    # Elements below diagonals are replaced by the elements above the diagonal.
                    M[i,j] = M[j,i]' 
                end
            end
        end
        return M
    end;
    
    w = eigvals(rhoA_Matrix()) # Eigenvalues of the reduced density matrix.
    #=
    The following loop calculates S = - sum \lamba_i * log(\lambda_i).
    =#
    
    # Array to store the log of the eigenvalues.
    DL = zeros(ComplexF64,2^LS)
    for i=1:length(w)
        if abs(w[i]) < 1.e-8 # Avoid zeros.
            continue
        else
            DL[i] = log(w[i])
        end
    end
    return real(-sum(w.*DL)) # S = -tr(rho *log(rho)).
end;

Bin2Dec(BinaryNumber) = parse(Int, string(BinaryNumber); base=2);
Dec2Bin(DecimalNumber) = string(DecimalNumber, base=2);

List = [i for i=0:2^L-1]; # List with numbers from 0 to 2^L-1.

#=
The following function converts all numbers in decimals in the above list 
 from 0 to 2^L -1 to binary.
=#

function List_Bin(Lst)
    
    l = []
    
    for i in Lst
        
        i_Bin = Dec2Bin(i)
        
        #=
        While converting numbers from decimal to binary, for example, 1
        is mapped to 1, to make sure that
        every numbers have N qubits in them, the following loop adds leading 
        zeros to make the
        length of the binary string equal to N. Now, 1 is mapped to 000.....1
        (string of length N).
        =#
        
        while length(i_Bin) < L
            i_Bin = "0"*i_Bin
        end
            
        # Puts the binary number in the list l after its length is L.
        push!(l,i_Bin)
    end
    return l
end;

#=
    The following function takes a binary string as input and rolls the qubits by one and
    returns the rolled string.
=#

Roll_String(Binary_String) = last(Binary_String)*Binary_String[1:L-1];

#=
    The following function takes a wavefunction as input and performs one roll
    on the qubits and returns the resultant wavefunction.
=#

function Psi_Roll(Initial_Psi)
    
    #=
        The following list contains all possible 2^N qubits after one roll 
        is performed on them.
        For example, the first position 0001 is changed to 1000.
    =#
    
    # Rolls every string in the list List by one qubit.
    Rl = [ Roll_String(i) for i in List_Bin(List)]
    
    #=
        The following list contains the previous list but in decimal numbers. For example,
        for N =4, the first position 1 is changed to 8.
    =#
    
    Rl_d = [Bin2Dec(i) for i in Rl]
    
    #=
        The following loop rearranges the coefficients of Psi after rolling. 
        For example, for N = 4, if the first coefficient 0001 is mapped to the
        eighth coefficient 1000 after one rotation of the qubits. 
        The coefficient of the rolled Psi in the i ^ th position is in the
        Rl_d[i] ^ th positon of the initial Psi.
    =#
    
    Psi_Rolled = []
    
    for i=1:2^L
        
        # Rearranging the coefficients according to the list l_d.
        push!(Psi_Rolled,Initial_Psi[Rl_d[i]+1])
        
        #= The addition of 1 to the index is necessary because Julia counts from 1,
           rather than 0. But we need to represent all the numbers from 1 to 16 using 
           four binary digits. So we started with the List = [0 to 15], then after
           rolling we just add 1 to each index to make it work for Julia.
        =#
    end
    return Psi_Rolled
end

#=
The following function performs specified number of rolls Num on the qubits.
=#

function N_Rolled(Num, Initial_Psi)
    
    if Num == 0 
        return Initial_Psi
    else
        
        s = Psi_Roll(Initial_Psi)
        for i=1:Num-1
            s = Psi_Roll(s)
        end
        return s
    end
end


function Average_Entropy(Initial_Psi)
    
    list_of_entropies = []
    #=
    The loop calculates all the entropies and returns a list containing them.
    =#
    for i=1:L
        S = Entropy(N_Rolled(i,Initial_Psi))
        push!(list_of_entropies,S)
    end
    return sum(list_of_entropies)/length(list_of_entropies)
end;

#Average_Entropy([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
#List_Bin(List)
#Psi_Roll([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

py"""
f = open('plot_data'+'.txt', 'w')
def Write_file(Noise, Energy, Entropy):
    f = open('plot_data'+'.txt', 'a')
    f.write(str(Noise) +'\t'+ str(Energy)+ '\t' + str(Entropy) +'\n')
"""

delta = 0.23
Op = Grover(delta)
EIGU = py"eigu"(Op)
X = string(delta)
Y = real(1im*log.(EIGU[1]))
V = EIGU[2]
    
for j=1:2^L
    py"Write_file"(delta, real(Y[j]), Average_Entropy(V[1:2^L,j:j]))
end
