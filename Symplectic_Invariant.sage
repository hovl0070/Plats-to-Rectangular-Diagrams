#This Computes Birman's Symplectic invariant for 2 and 3 bridge links.  It has not been generalized to higher bridge numbers yet.  I have not discovered a clever pattern for the lifts yet.  

def two_bridge_symplectic_invariant(braid1,braid2):
    strands=4
    ######### Initialization ###########
    #The idea is to realize some pattern in the matrix representations so that this function can be called for any bridge numbers.  That is   why there are some variables here. 
    n=strands-2  #the square matrices representing the braids are strands-2.
    k=n//2     #these are the 4 k x k submatrices
    k_identity_matrix=matrix.identity(k)
    k_zero_matrix=matrix.zero(k)
    k_ones_matrix=matrix.ones(k)

    braidlifts=list(range(strands-1))
    
    ######## Creating matrices that are lifts of the braids ### 
    ### This should be a for loop  going through range(strands-1) once I work out the pattern for higher bridge numbers as well
    B1=matrix([[1,1],[0,1]])
    braidlifts[0]=B1    

    B2=matrix([[1,0],[1,1]])
    braidlifts[1]=B2    

    B3=matrix([[1,1],[0,1]])
    braidlifts[2]=B3

   #### Now we multiply our matrices according to the braid words given. 
    B1_matrix=matrix.identity(n) #initializing our matrix
    for i in range(len(braid1)):
        if braid1[i]>=0:         #its a positive multiplication
            B1_matrix=B1_matrix*braidlifts[braid1[i]-1] #multiply the new word on the right
        else: 
            temp=-braid1[i]-1
            B1_matrix=B1_matrix*braidlifts[temp].inverse() #multiply the new word on the right
     
    B2_matrix=matrix.identity(n) #initializing our matrix
    for i in range(len(braid2)):
        if braid2[i]>=0:         #its a positive multiplication
            B2_matrix=B2_matrix*braidlifts[braid2[i]-1] #multiply the new word on the right
        else: 
            temp=-braid2[i]-1
            B2_matrix=B2_matrix*braidlifts[temp].inverse() #multiply the new word on the right
            
    P1=B1_matrix.submatrix(k,0,k,n-k)
    Q1=B1_matrix.submatrix(k,n-k,k,n-k)
    R1=B1_matrix.submatrix(0,0,n-k,n-k)

    P2=B2_matrix.submatrix(k,0,k,n-k)
    Q2=B2_matrix.submatrix(k,n-k,k,n-k)
    R2=B2_matrix.submatrix(0,0,n-k,n-k)

    
    D1, U1, V1 = P1.smith_form(integral=True)
    D2, U2, V2 = P2.smith_form(integral=True)
    p1=P1[k-1,k-1]
    p2=P2[k-1,k-1]
    if p1 !=p2:
        print('The p-values are not the same! Your braids are not the same link type.  They are not even mutants of each other!')
        print('p-value of the first braid is:',p1)
        print('p-value of the second braid is:',p2)
    else:    
        print('braid_1: p1=',p1, "Det(R1):",det(R1), "Det(Q1):,",det(Q1))
        print('braid_2: p2=',p2, "Det(R2):",det(R2), "Det(Q2):,",det(Q2))
    return [B1_matrix,B2_matrix]



def three_bridge_symplectic_invariant(braid1,braid2):
    strands=6
    ######### Initialization ###########
    #The idea is to realize some pattern in the matrix representations so that this function can be called for any bridge numbers.  That is   why there are some variables here. 
    n=strands-2  #the square matrices representing the braids are strands-2.
    k=n//2     #these are the 4 k x k submatrices
    k_identity_matrix=matrix.identity(k)
    k_zero_matrix=matrix.zero(k)
    k_ones_matrix=matrix.ones(k)

    braidlifts=list(range(strands-1))
    
    ######## Creating matrices that are lifts of the braids ### 
    ### This should be a for loop  going through range(strands-1) once I work out the pattern for higher bridge numbers as well
    A=matrix.zero(k)
    A[[0],[0]]=-1   
    B1=k_identity_matrix.augment(A)
    temp=k_zero_matrix.augment(k_identity_matrix)
    B1=B1.stack(temp)
    braidlifts[0]=B1    

    A=matrix.zero(k)
    A[[0],[0]]=1
    temp=A.augment(k_identity_matrix)
    B2=k_identity_matrix.augment(k_zero_matrix)
    B2=B2.stack(temp)
    braidlifts[1]=B2    

    A=matrix.ones(k)-2*matrix.identity(k)
    B3=k_identity_matrix.augment(A)
    temp=k_zero_matrix.augment(k_identity_matrix)
    B3=B3.stack(temp)
    braidlifts[2]=B3

    A=matrix.zero(k)
    A[[k-1],[k-1]]=1   
    temp=A.augment(k_identity_matrix)
    B4=k_identity_matrix.augment(k_zero_matrix)
    B4=B4.stack(temp)
    braidlifts[3]=B4    

    A=matrix.zero(k)
    A[[k-1],[k-1]]=-1   
    temp=k_identity_matrix.augment(A)
    B5=k_zero_matrix.augment(k_identity_matrix)
    B5=temp.stack(B5)
    braidlifts[4]=B5    

   #### Now we multiply our matrices according to the braid words given. 
    B1_matrix=matrix.identity(n) #initializing our matrix
    for i in range(len(braid1)):
        if braid1[i]>=0:         #its a positive multiplication
            B1_matrix=B1_matrix*braidlifts[braid1[i]-1] #multiply the new word on the right
        else: 
            temp=-braid1[i]-1
            B1_matrix=B1_matrix*braidlifts[temp].inverse() #multiply the new word on the right
     
    B2_matrix=matrix.identity(n) #initializing our matrix
    for i in range(len(braid2)):
        if braid2[i]>=0:         #its a positive multiplication
            B2_matrix=B2_matrix*braidlifts[braid2[i]-1] #multiply the new word on the right
        else: 
            temp=-braid2[i]-1
            B2_matrix=B2_matrix*braidlifts[temp].inverse() #multiply the new word on the right
            
    P1=B1_matrix.submatrix(k,0,k,n-k)
    Q1=B1_matrix.submatrix(k,n-k,k,n-k)
    R1=B1_matrix.submatrix(0,0,n-k,n-k)

    P2=B2_matrix.submatrix(k,0,k,n-k)
    Q2=B2_matrix.submatrix(k,n-k,k,n-k)
    R2=B2_matrix.submatrix(0,0,n-k,n-k)

    
    D1, U1, V1 = P1.smith_form(integral=True)
    D2, U2, V2 = P2.smith_form(integral=True)
    p1=P1[k-1,k-1]
    p2=P2[k-1,k-1]
    if p1 !=p2:
        print('The p-values are not the same! Your braids are not the same link type.  They are not even mutants of each other!')
        print('p-value of the first braid is:',p1)
        print('p-value of the second braid is:',p2)
    else:    
        print('braid_1: p1=',p1, "Det(R1):",det(R1), "Det(Q1):,",det(Q1))
        print('braid_2: p2=',p2, "Det(R2):",det(R2), "Det(Q2):,",det(Q2))
    return [B1_matrix,B2_matrix]