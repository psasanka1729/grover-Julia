using Random
using LinearAlgebra
using SparseArrays
using PyCall
using Statistics

L = 14;

py"""
f = open('Grover_gates_data'+'.txt', 'w')
def Write_file1(X, Y, Z):
    f = open('Grover_gates_data'+'.txt', 'a')
    f.write(str(X) + '\t' + str(Y)+ '\t' + str(Z) +'\n')
"""

Identity(dimension) = 1* Matrix(I, dimension, dimension);

function Write_Gates_to_File(L)
NOISE = 2*rand(Float64,Write_Gates_to_File(L)).-1;
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

Write_Gates_to_File(14)
