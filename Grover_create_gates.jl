using Random
using LinearAlgebra
using SparseArrays
using PyCall
using Statistics

L = 10;
Random.seed!(3000);

py"""
f = open('Grover_gates_data'+'.txt', 'w')
def Write_file1(X, Y, Z):
    f = open('Grover_gates_data'+'.txt', 'a')
    f.write(str(X) + '\t' + str(Y)+ '\t' + str(Z) +'\n')
"""

Identity(dimension) = 1* Matrix(I, dimension, dimension);

function Write_Gates_to_File(L)
    Gate_count = 0
    U0_XHL_Gates = []
    for i = 1:L-1
        push!(U0_XHL_Gates,["X",i])
    end    
    push!(U0_XHL_Gates,["H",L])

    U0_XHR_Gates = []
    for i = 1:L-1
        push!(U0_XHR_Gates,["X",i])
    end
    push!(U0_XHR_Gates,["H",L])
    
    MCX = sparse(Identity(2^L));
    
    XHL_Matrix = sparse(Identity(2^L))
    for i in U0_XHL_Gates
        if i[1] == "H"   
            
            py"Write_file1"("H",0.0,i[2])
            Gate_count +=1
            
        elseif i[1] == "X"

            py"Write_file1"("X",0.0,i[2])
            Gate_count +=1
 
        end
    end
    
    #= Constructing the multicontrolled Toffoli gate. =# 
    # C_1.
    for i = 1:L-2
        for j = 1:i

            py"Write_file1"(pi/2^j, L-i, L-i+j)
            Gate_count +=1
            
        end
    end

    # C_2.
    for i = 2:L
        
        py"Write_file1"(pi/2^(i-2), 1, i)
        Gate_count +=1

    end

    # C3 = - C1.
    for i = L-2:-1:1
        for j = i:-1:1

            py"Write_file1"(-pi/2^j, L-i, L-i+j)
            Gate_count +=1

            
        end
    end

    # C_4.
    for i = 1:L-3
        for j = 1:i

            py"Write_file1"(pi/2^j, L-i-1, L-i-1+j)
            Gate_count +=1
   
        end    
    end

    # C_5.
    for i = 2:L-1

        py"Write_file1"(-pi/2^(i-2), 1, i)
        Gate_count +=1
  
    end

    # C6 = - C4.
    for i = L-3:-1:1
        for j = i:-1:1

            py"Write_file1"(-pi/2^j, L-i-1, L-i-1+j)
            Gate_count +=1
                
        end    
    end

    XHR_Matrix = sparse(Identity(2^L))
    for j in U0_XHR_Gates
        if j[1] == "H"

            py"Write_file1"("H", 0.0,j[2])
            Gate_count +=1

                
        elseif j[1] == "X"
            
            py"Write_file1"("X",0.0,j[2])
            Gate_count +=1
  
        end
    end

    #U0_matrix = sparse(XHL_Matrix*MCX*XHR_Matrix)    

    
    #= Ux matrix. =#
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
    

    # Creating an empty matrix to store the Ux matrix.
    MCX = sparse(Identity(2^L));
    
    XHL_Matrix = sparse(Identity(2^L))
    for i in Ux_XHL_Gates
        
        if i[1] == "H"

            py"Write_file1"("H", 0.0,i[2])
            Gate_count +=1
            
        elseif i[1] == "X"

            py"Write_file1"("X", 0.0,i[2])
            Gate_count +=1
     
        end
    end
    
    #= Contructing the multicontrolled Toffoli gate. =#
    # C_1.
    for i = 1:L-2
        for j = 1:i
            
            py"Write_file1"(pi/2^j, L-i, L-i+j)
            Gate_count +=1
                
        end
    end

    # C_2.
    for i = 2:L
        
        py"Write_file1"(pi/2^(i-2), 1, i)
        Gate_count +=1
  
    end

    # C3 = - C1.
    for i = L-2:-1:1
        for j = i:-1:1

            py"Write_file1"(-pi/2^j, L-i, L-i+j)
            Gate_count +=1
   
        end
    end

    # C_4.
    for i = 1:L-3
        for j = 1:i

            py"Write_file1"(pi/2^j, L-i-1, L-i-1+j)
            Gate_count +=1
   
        end    
    end

    # C_5.
    for i = 2:L-1

        py"Write_file1"(-pi/2^(i-2), 1, i)
        Gate_count +=1
    
    end

    # C6 = - C4.
    for i = L-3:-1:1
        for j = i:-1:1

            py"Write_file1"(-pi/2^j, L-i-1, L-i-1+j)
            Gate_count +=1
   
        end    
    end


    XHR_Matrix = sparse(Identity(2^L))
    for j in Ux_XHR_Gates
        if j[1] == "H"          
            
            py"Write_file1"("H", 0.0,j[2])
            Gate_count +=1
   
        elseif j[1] == "X"         
            
            py"Write_file1"("X", 0.0,j[2])
            Gate_count +=1
   
        end
    end
    return Gate_count
end; 

NOISE = 2*rand(Float64,Write_Gates_to_File(L)).-1;

py"""
f = open('Noise_data'+'.txt', 'w')
def Write_file2(Z):
    f = open('Noise_data'+'.txt', 'a')
    f.write(str(Z) +'\n')
"""

for i in NOISE
    py"Write_file2"(i)
end

#real(round.(collect(Ux_matrix*U0_matrix),digits=2))

#List_of_H[8]

#=
H3 = List_of_H[8] # CX
Noise = H3[1]
Control_Qubit = int(H3[2])
Target_Qubit = int(H3[3])
Matrices = Dict("I" => I,"U" => (I-CX(0.0))/2, "PI_1" => (I-Z)/2)
p1 = fill("I", L)
p1[Control_Qubit] = "PI_1"
p1[Target_Qubit] = "U"
H_k3 = Matrices[p1[1]]
for i = 2:L
    H_k3 = kron(H_k3,Matrices[p1[i]])
end 
exp(-1im*(H3[1])*H_k3)=#

#List_of_V[8]

U_0 = [-1 0 0 0; 0 1 0 0; 0 0 1 0;0 0 0 1];
A = ones(2^L,2^L);
U_x = (2/2^L)*A-Identity(2^L); # 2\s><s|-I
G_exact = U_x*U_0
#U_x

#G_from_U = -round.(List_of_V[10]*List_of_V[9]*List_of_V[8]*List_of_V[7]*List_of_V[6]*List_of_V[5]*List_of_V[4]*List_of_V[3]*List_of_V[2]*List_of_V[1],digits=2)

#=
G_from_H = Identity(2^L);
G_from_H_list = []
for U_k in reverse(List_of_H) # reverse is necessary.
        if U_k[1] == "H"
            Noise = U_k[2] # delta*epsilon.
            Qubit = U_k[3] # qubit.
            H_k = Matrix_Gate([1 0;0 1]-1/sqrt(2)*[1 1;1 -1],Qubit) #= H_had = I2-Had. =#
            G_from_H = G_from_H*exp(-1im*(pi/2+Noise)*collect(H_k))
            push!(G_from_H_list,exp(-1im*(pi/2+Noise)*collect(H_k)))
        elseif U_k[1] == "X"
            Noise = U_k[2] # delta*epsilon.
            Qubit = U_k[3] # qubit.
            H_k = Matrix_Gate([1 1;1 1],Qubit) #= H_X = X+I2. =#
            G_from_H = G_from_H*exp(-1im*(pi/2+Noise)*collect(H_k))
            push!(G_from_H_list,exp(-1im*(pi/2+Noise)*collect(H_k)))
        else
            Angle = U_k[1]
            Control_Qubit = int(U_k[2])
            Target_Qubit = int(U_k[3])
            #= H = ((I-Z)/2)_c \otimes ((I+X)/2)_t.=#
            Matrices = Dict("I" => I,"U" => (I-CX(0.0))/2, "PI_1" => (I-Z)/2)
            p1 = fill("I", L)
            p1[Control_Qubit] = "PI_1"
            p1[Target_Qubit] = "U"
            H_k = Matrices[p1[1]]
            for i = 2:L
                H_k = kron(H_k,Matrices[p1[i]]);
            end   
            G_from_H = G_from_H*exp(-1im*(Angle)*collect(H_k))
            push!(G_from_H_list,exp(-1im*(Angle)*collect(H_k)))
        end
end
real(round.(-G_from_H,digits=8))=#

#=G_from_V = Identity(2^L)
for k in List_of_V
    G_from_V = k*G_from_V
end
G_from_V=#

#-real(Matrix_Gate(Hadamard(0.0),1)*Matrix_Gate(CX(0.0),1)*CU(CX(0.0),1,2)*(Matrix_Gate(CX(0.0),1))*(Matrix_Gate(Hadamard(0.0),1)))

#G_exact

#=
k = 1
    f_k = Identity(2^L);
        for i = k:length(List_of_U)-1
            f_k = f_k * collect(List_of_U[length(List_of_U)-i+k])
        end 
        
        #= Corresponding H for the kth term. =#
        U_k = List_of_H[k]
        if U_k[1] == "H"
            Noise = U_k[2] # delta*epsilon.
            Qubit = U_k[3] # qubit.
            H_k = Matrix_Gate(I2-H,Qubit) #= H_had = I2-Had. =#
        elseif U_k[1] == "X"
            Noise = U_k[2] # delta*epsilon.
            Qubit = U_k[3] # qubit.
            H_k = Matrix_Gate(I+CX(0.0),Qubit) #= H_X = X+I2. =#
 
        else
            Noise = U_k[1]
            Control_Qubit = int(U_k[2])
            Target_Qubit = int(U_k[3])
            #= H = ((I-Z)/2)_c \otimes ((I+X)/2)_t.=#
            Matrices = Dict("I" => I,"U" => (I+CX(0.0))/2, "PI_1" => (I-Z)/2)
            p1 = fill("I", L)
            p1[Control_Qubit] = "PI_1"
            p1[Target_Qubit] = "U"
            H_k = Matrices[p1[1]]
            for i = 2:L
                H_k = kron(H_k,Matrices[p1[i]])
            end                                 
        end
K1 = real(round.(f_k*H_k*(f_k'),digits=2))
#real(round.(H_k,digits=2))
#real(round.(H_k,digits=2))=#

#=
# 10 th term.
H10 = List_of_H[10] # H
Noise = H10[2] # delta*epsilon.
Qubit10 = H10[3] # qubit.
H_k10 = collect(Matrix_Gate(I-H,Qubit10))
# 9 th term.
H9 = List_of_H[9] # X
Noise = H9[2] # delta*epsilon.
Qubit9 = H9[3] # qubit.
H_k9 = List_of_U[10]*Matrix_Gate(-[0 1;1 0]+I,Qubit9)*(List_of_U[10])'
real(round.(H_k9,digits=2))
# 8 th term.
H8 = List_of_H[8] # CX
Noise = H8[1]
Control_Qubit = int(H8[2])
Target_Qubit = int(H8[3])
Matrices = Dict("I" => I,"U" => (I+CX(0.0))/2, "PI_1" => (I-Z)/2)
p1 = fill("I", L)
p1[Control_Qubit] = "PI_1"
p1[Target_Qubit] = "U"
H_k8 = Matrices[p1[1]]
for i = 2:L
    H_k8 = kron(H_k8,Matrices[p1[i]])
end 
H_k8 = List_of_U[10]*List_of_U[9]*H_k8*List_of_U[9]'*List_of_U[10]'
real(round.(H_k8,digits=2))
# 7 th term
H7 = List_of_H[7] # X
Noise = H7[2] # delta*epsilon.
Qubit7 = H7[3] # qubit.
H_k7 = List_of_U[10]*List_of_U[9]*List_of_U[8]*Matrix_Gate(CX(0.0)+I,Qubit7)*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k7,digits=2))
# 6 th term.
H6 = List_of_H[6] # H
Noise = H6[2] # delta*epsilon.
Qubit6 = H6[3] # qubit.
H_k6 = List_of_U[10]*List_of_U[9]*List_of_U[8]*List_of_U[7]*Matrix_Gate(I-H,Qubit6)*List_of_U[7]'*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k6,digits=2))
# 5 th term.
H5 = List_of_H[5] # H
Noise = H5[2] # delta*epsilon.
Qubit5 = H5[3] # qubit.
H_k5 = List_of_U[10]*List_of_U[9]*List_of_U[8]*List_of_U[7]*List_of_U[6]*Matrix_Gate(I-H,Qubit5)*List_of_U[6]'*List_of_U[7]'*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k5,digits=2))
# 4 th term
H4 = List_of_H[4] # X
Noise = H4[2] # delta*epsilon.
Qubit4 = H4[3] # qubit.
H_k4 = List_of_U[10]*List_of_U[9]*List_of_U[8]*List_of_U[7]*List_of_U[6]*List_of_U[5]*
Matrix_Gate(CX(0.0)+I,Qubit4)*List_of_U[5]*List_of_U[6]*List_of_U[7]'*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k4,digits=2))
# 3 rd term
H3 = List_of_H[3] # CX
Noise = H3[1]
Control_Qubit = int(H3[2])
Target_Qubit = int(H3[3])
Matrices = Dict("I" => I,"U" => (I+CX(0.0))/2, "PI_1" => (I-Z)/2)
p1 = fill("I", L)
p1[Control_Qubit] = "PI_1"
p1[Target_Qubit] = "U"
H_k3 = Matrices[p1[1]]
for i = 2:L
    H_k3 = kron(H_k3,Matrices[p1[i]])
end 
H_k3 = List_of_U[10]*List_of_U[9]*List_of_U[8]*List_of_U[7]*List_of_U[6]*List_of_U[5]*List_of_U[4]*H_k3*List_of_U[4]'*List_of_U[5]'*List_of_U[6]'*List_of_U[7]'*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k3,digits=2))
# 2 nd term
H2 = List_of_H[2] # H
Noise = H2[2] # delta*epsilon.
Qubit2 = H2[3] # qubit.
H_k2 = List_of_U[10]*List_of_U[9]*List_of_U[8]*List_of_U[7]*List_of_U[6]*List_of_U[5]*List_of_U[4]*List_of_U[3]*Matrix_Gate(I-H,Qubit2)*List_of_U[3]'*List_of_U[4]'*List_of_U[5]'*List_of_U[6]'*List_of_U[7]'*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k2,digits=2))
# 1 st term
H1 = List_of_H[1] # X
Noise = H1[2] # delta*epsilon.
Qubit1 = H1[3] # qubit.
H_k1 = List_of_U[10]*List_of_U[9]*List_of_U[8]*List_of_U[7]*List_of_U[6]*List_of_U[5]*List_of_U[4]*List_of_U[3]*List_of_U[2]*Matrix_Gate(CX(0.0)+I,Qubit1)*List_of_U[2]'*List_of_U[3]'*List_of_U[4]'*List_of_U[5]'*List_of_U[6]'*List_of_U[7]'*List_of_U[8]'*List_of_U[9]'*List_of_U[10]'
real(round.(H_k1,digits=2))=#

#=
Hks = [H_k10,H_k9,H_k8,H_k7,H_k6,H_k5,H_k4,H_k3,H_k2,H_k1];
Zz = zeros(2^L,2^L)
h_effective = Noise_Used[10]*H_k10+Noise_Used[9]*H_k9+Noise_Used[8]*H_k8+Noise_Used[7]*H_k7+Noise_Used[6]*H_k6+Noise_Used[5]*H_k5+Noise_Used[4]*H_k4+Noise_Used[3]*H_k3+Noise_Used[2]*H_k2+Noise_Used[1]*H_k1
h_effective = DELTA*h_effective
E_effective = eigen(h_effective[3:2^L,3:2^L]).values
E_effective  = sort(real(E_effective),rev = true)=#

#=py"""
f = open('new_eigenvalues_data'+'.txt', 'w')
def Write_file2(delta, effective, exact):
    f = open('new_eigenvalues_data'+'.txt', 'a')
    f.write(str(delta) + '\t' + str(effective)+ '\t' + str(exact) +'\n')
"""
Num = 100;
for i = 1:Num
    delta = 0.1*(i/Num)
    EE = MyEigenvalues(delta);
    effective = EE[1]
    exact = EE[2]
    for j = 1:2^L-2
        py"Write_file2"(delta,effective[j],exact[j])
        #println(effective[j])
        #println(exact[j])
    end
end=#
