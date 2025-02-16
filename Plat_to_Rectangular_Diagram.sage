def Plat_to_Rectangular_Diagram(W):
    #INPUT: the braid word whose plat closure is some link.   
    #OUTPUT: A grid G which can then be plotted using draw_grid(G) from  GridPyM
    
    positive_reduced_sawtooth_array=braid_to_reduced_positive_braid(W)
    syllablelist=syllable_string_partitioning(positive_reduced_sawtooth_array)
    Vertices=plat2rectdiag(syllablelist)
    G=vertexset_to_XOGrid(Vertices)

    return G


def braid_to_positive_braid(W):
    #Input is an array [1 -2 3 5 -1 -2] etc.
    #Output is a positive braid.  This could be adjusted to output just the positive array!
    n=max(abs(i)+1 for i in W) #This is the braid group
    L=list(range(1,n)) #This will be used to generate the positive replacement words
    positivewordreplacement=[]
    for j in range(1,n):
        positivewordreplacement.append(L[j-(n+1)::-1]+L+L[:j-1:-1]) #This array holds in each index the replacement positive word
        PositiveW=[];
    for i in W:
        if i>0:
            PositiveW.append(i)
        else:
             PositiveW.append(positivewordreplacement[i])
    PositiveW=flatten(PositiveW) #This is the positive braid word replacing the input W
    Braidgroup=BraidGroup(n,'s'); #The generators are s_0-s_(n-2) 
    PositiveBraidW=Braidgroup(PositiveW)
    return PositiveBraidW

def braid_to_reduced_positive_braid(W):
    #Input is an array [1 -2 3 5 -1 -2] etc.
    #Output is a reduced positive braid. This is a positive braid with no garside elements.      The link type is unchanged.
    n=max(abs(i)+1 for i in W) #This is the braid group
    if n%2==0:  #if n is even then we are fine, but if the last strand of the plat is free then we need to make sure we are still in an even braid group.
        n=n
    else:
        n=n+1
        
    B = BraidGroup(n) #this is braid group
    garside=[j for i in range(n, -1, -1) for j in range(1,n) if i > j] 
    garside=B(garside) #this generates the garside element in the braid group Bn
    braid_word=B(W)    #this takes our array and makes it a braid
    L1=braid_word.left_normal_form() #this is the left normal form.
    if len(L1)!=1:   #there is more than one factor involved in this braid
        first_factor_of_L1=L1[0]
        if first_factor_of_L1.lcm(garside)==garside: #this means the first factor is a power of the garside element
            L1_new = L1[1:]   #I have removed the power of the garside element
        else:
            L1_new=L1
    else:
        L1_new=L1
    L1_new=list(L1_new) #this returns a list with each entry a factor of the braid and no garside elements
    postive_reduced_braid_list=[]
    for factors in L1_new:
        factors=list(factors.Tietze())  #this turns the braid back into an list
        postive_reduced_braid_list.append(factors)
    postive_reduced_braid_list=flatten(postive_reduced_braid_list)    
    #return postive_reduced_braid_list
     #This will take the array of integers and put it into sawtooth form. 
    sawtooth=0
    while sawtooth==0:
        for i in range(len(postive_reduced_braid_list)-1):
            if postive_reduced_braid_list[i+1]>postive_reduced_braid_list[i]+1:
                temp=postive_reduced_braid_list[i]
                postive_reduced_braid_list[i]=postive_reduced_braid_list[i+1]
                postive_reduced_braid_list[i+1]=temp
                sawtooth=0
                break
            else:
                sawtooth=1
    return postive_reduced_braid_list    #the output is a positive reduced braid in sawtooth form

#This takes a postive reduced braid array in sawtooth form and returns just the increasing syllables.  

def syllable_string_partitioning(array):
    syllable=[]
    syllables=[]
    for i in range(len(array)-1):
        if (array[i]+1 == array[(i+1)]):
            syllable.append(array[i])
        else:
            syllable.append(array[i])
            syllables.append(syllable)
            syllable=[]
    if array[len(array)-1]-1==array[len(array)-1]:
        syllable.append( array[len(array)])
    else: syllables.append([array[len(array)-1]])
    W_input=[[syllables[i][0],syllables[i][-1]] for i in range(len(syllables))]
    return W_input #the output is the increasing syllables.  This is input for plat2rectdiag function



def plat2rectdiag(syllablelist):
    # This is figuring out how many bridges there are. This will tell us how many strands there are. 
    syllableends = [j for i, j in syllablelist]
    maxcrossing: int = max(syllableends)
    if maxcrossing % 2 == 1:
        n = (maxcrossing + 1) / 2
    else:
        n = (maxcrossing + 2) / 2
    n = int(n)

    # This is initializing the rectangular diagram: strands are initially in the first 2n positions, and each bridge starts on the next column.
    strands = list(range(1, (2 * n) + 1))
    vertexset = []
    i = 1
    while i <= n:
        vertexset.append((i, 2 * i - 1))
        vertexset.append((i, 2 * i))
        i += 1


    # For each syllable, we add two new vertices to the vertexset and update where each strand is
    # print(strands)
    k=1
    while syllablelist:
        # print('step ', k)
        currentsyllable = syllablelist[0]
        
        #Adds two new vertices
        newvertextop = (n + k, strands[currentsyllable[0] - 1])
        newvertexbot = (n + k, strands[currentsyllable[1]] + 1)
        
        # x and y values in vertex set without the new vertices
        vsx = [i for i, j in vertexset]
        vsy = [j for i, j in vertexset]
        
        shiftedindices = list(filter(lambda x: vsy[x] >= strands[currentsyllable[1]] + 1, range(len(vsy))))
        for i in shiftedindices:
            vsy[i] = vsy[i] + 1
        vertexset = list(zip(vsx, vsy))
        vertexset = vertexset + [newvertextop, newvertexbot]
        # print(vertexset)
        

        # Now we need to update where the strands are. It's an annoying list manipulation, but nothing complicated.
        strands.insert(currentsyllable[1], strands[currentsyllable[1]])
        strandchange = list(range(currentsyllable[1] + 1, len(strands)))
        for i in strandchange:
            strands[i] += 1
        strands.pop(currentsyllable[0] - 1)
        # print(strands)
        syllablelist.pop(0)
        k += 1
    # print('end of for loop')
    # print('k = ',1)
    # print(strands)

    # Now we need to close off the plat where strands 1&2 are, 3&4 are, etc.
    for i in range(1, n + 1):
        vertexset.append((n+ k + i - 1, strands[2 * i - 2]))
        vertexset.append((n + k + i - 1, strands[2 * i - 1]))

    # Lastly, we'll make our vertex set go in the negative y direction, it'll make our pictures look how we want them to.
    vsx = [i for i, j in vertexset]
    vsy = [j for i, j in vertexset]
    negvsy = [-x for x in vsy]
    vertexset = list(zip(vsx, negvsy))

   

    return vertexset

def rectangulardiagram_to_orientedrectangulardiagram(V):
    #Input is the vertex set of the rectangular diagram
    #The goal is to add an additional entry to the vertex set, either an 0 or 1 this will serve as the orientation. We need on each horizontal and vertical a 0 and a 1.  
    #Add a 0 to the first vertex in the vertex set
    oriented_vertex_set=[V[0]+(0,)]
    #Delete that element from the vertex set
    V.pop(0)
    
    while V: #While the vertex set is non-empty we will add a 0 or 1 to each vertex and then add that to the oriented_vertex_set
        
            xindex=[j for j in range(len(V)) if (V[j][0]==oriented_vertex_set[-1][0])] #This will find the index in vertex set of vertex sharing the horizontal with the last element added to the oriented_vertex_set.
            yindex=[j for j in range(len(V)) if (V[j][1]==oriented_vertex_set[-1][1])] #This will find the index in vertex set of vertex sharing the vertical with the last element added to the oriented_vertex_set.
            
            if xindex:  # If there is a vertex on the horizonal strand that is the one we will add to the oriented_vertex_set
                oriented_vertex_set.append(V[xindex[0]]+((oriented_vertex_set[-1][2]+1) % 2,)) #This adds the updated oriented vertex 
                V.pop(xindex[0]) #That shared horizontal vertex is deleted. 
            elif yindex:     #Maybe you do not have a vertex sharing the same horizontal in vertex_set (its already in oriented_vertex_set) so find the vertex sharing a vertical
                oriented_vertex_set.append(V[yindex[0]]+((oriented_vertex_set[-1][2]+1) % 2,)) #add this to the oriented vertex set
                V.pop(yindex[0])
            else:
                first_vertex=V[0]
                oriented_vertex_set.append(first_vertex+(0,)) #Add a 0 to the first vertex in the vertex set.  This will be a new link component
                V.remove(first_vertex) 
                #break
    return oriented_vertex_set


def vertexset_to_XOGrid(vertexset):
    V1=rectangulardiagram_to_orientedrectangulardiagram(vertexset)
    output = sorted(V1, key=lambda x: x[-1])
    X = [(output[i][0], output[i][1]) for i in range(len(output)) if output[i][2] == 0]
    X=sorted(X, key=lambda x: x[-1])
    Xlist=[X[i][0]-1 for i in range(len(X))]
    O =[(output[i][0], output[i][1]) for i in range(len(output)) if output[i][2] == 1]
    O=sorted(O, key=lambda x: x[-1])
    Olist=[O[i][0]-1 for i in range(len(O))]
    G=[Xlist,Olist]
    return G


def drawrectdiagwrite(vertexset):
#### This will save the tikzfile to your computer as a .txt file
    import numpy as np
    import sys
    sys.stdout = open('RectangularDiagramtizkfile.txt','w')
    print("\\begin{figure}")
    print("\\centering")
    print("\\resizebox{\\columnwidth}{!}")
    print('{')
    print("\\begin{tikzpicture}")
    # First, we draw the horizontal lines remembering that our y's are negative
    vsy = [j for i, j in vertexset]
    miny: int = min(vsy)
    vsyarray = np.array(vsy)
    for i in range(miny, 0):
        currentvertices = np.where(vsyarray == i)[0]
        vertex1 = vertexset[currentvertices[0]]
        vertex2 = vertexset[currentvertices[1]]
        print("\\draw[black, very thick] ", vertex1, ' -- ', vertex2, ';')

    print('')
    print('')

    # Next, we need to draw the white rectangles because we don't want to figure out how to use the knots package in tikz.
    vsx = [i for i, j in vertexset]
    maxx: int = max(vsx)
    vsxarray = np.array(vsx)
    for i in range(1, maxx + 1):
        currentvertices = np.where(vsxarray == i)[0]
        vertex1 = vertexset[currentvertices[0]]
        vertex2 = vertexset[currentvertices[1]]
        realv1 = (vertex1[0] - .05, vertex1[1] - .25)
        realv2 = (vertex2[0] + .05, vertex2[1] + .25)
        print("\\filldraw[white, very thick] ", realv1, ' rectangle ', realv2, ';')

    print('')
    print('')

    # Next, we draw a small circle around each vertex, and draw red ones on red edges
    for i in vertexset:
        print("\\filldraw[black] ", i, ' circle (1.5pt) node[anchor=west]{};')

    print('')
    print('')

    # Lastly, we draw the vertical lines, where we need to check for redness.
    for i in range(1, maxx + 1):
        currentvertices = np.where(vsxarray == i)[0]
        vertex1 = vertexset[currentvertices[0]]
        vertex2 = vertexset[currentvertices[1]]
      
            print("\\draw[black, very thick] ", vertex1, ' -- ', vertex2, ';')

    print('\\end{tikzpicture}')
    print('}')
    print('\\end{figure}')

