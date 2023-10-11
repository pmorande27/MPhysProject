import numpy as np
import SUroutines as SU2

class Chiral(object):
    def __init__(self,N, beta, D = 2) -> None:

        self.N = N

        self.U = np.zeros((N,N,4),dtype=float) 
        self.beta = beta
        a = np.arange(N*N).reshape((N,N))

        for i in range(N):

            for j in range(N):
                    self.U[i,j] = [1,0,1,0]

        self.NN_mask = self.make_NN_mask()
    def make_NN_mask(self):
        '''Makes mask to apply to phi or pi which then gives the matrix parameter values of the nearest neighbors (NN) for each lattice site.
        Hence phi[self.NN_mask] is of shape (N,N,#neighbors,#parameters) i.e (N,N,4,4).
        
        Returns:
        NN_mask: tuple
            tuple of two (N,N,4,1) arrays, each giving the row and column coordinate for all nearest neighbors
        '''
     
        # make a (N,N,2) array storing the row and col indices of each lattice sites
        grid = np.indices((self.N,self.N)) 
        lattice_coords = grid.transpose(1,2,0)

        # shift lattice coordinates by 1 such that the coordinates at (i,j) are those of the right, left, top, and bottom neighbor of lattice site (i,j)
        # rolling axis=1 by -1 means all columns are moved one step to the left with periodic bcs. Hence value of resulting array at (i,j) is (i,j+1), i.e the coordinates of the right neighbor.
        # all of shape (N,N,2)
        right_n = np.roll(lattice_coords, -1, axis=1)
        left_n = np.roll(lattice_coords, 1, axis=1)
        top_n = np.roll(lattice_coords, 1, axis=0)
        bottom_n = np.roll(lattice_coords, -1, axis=0)

        # for each lattice site, for each neighbor, store row and column coordinates
        # order of neighbors: right, left, top, bottom
        NN = np.empty((self.N,self.N,4,2), dtype=int)
        NN[:,:,0,:] = right_n # row and col indices of right neighbors
        NN[:,:,1,:] = left_n
        NN[:,:,2,:] = top_n
        NN[:,:,3,:] = bottom_n

        # make mask to index phi
        # separate the row and column neighbor coordinates for each lattice site: (N,N,4,1)
        NN_rows = NN[:,:,:,0]
        NN_cols = NN[:,:,:,1]
        NN_mask = (NN_rows, NN_cols)

        return NN_mask 
       
    def action_2(self, phi):
        '''
        Computes the action for lattice configuration phi
        phi: (N,N,4) array
            parameter values of SU(2) matrices at each lattice site

        Returns
        S: float
            the action
        '''
        #phi_hc = np.matrix(phi).getH()
        phi_NN = phi[self.NN_mask] # (N,N,4,4): containing the 4 paras of each of the 4 NN
        phi_hc = SU2.hc(phi)
        # sum over lattice unit vectors: to the right and up. Hence only need right and top NN, stored at position 0,3 respectively
        G = np.zeros((self.N,self.N))
        print(phi)
        for i in [0,3]:
        
            A = SU2.dot(phi_hc, phi_NN[:,:,i,:])
            G += SU2.tr(A + SU2.hc(A)) # when getting UFuncTypeError, check that dtype of G and SU2.tr is the same (float64 by default)

        # sum over lattice sites    
        S = -1/2 * self.beta * np.sum(G)

        return S

    @staticmethod
    def action(U, beta):
        G = 0
        top_neighbours = np.roll(U,-1,axis = 0)
        right_neighbours = np.roll(U,1,axis = 1)
        for t in range(len(U)):
            for x in range(len(U)):
                G+= np.dot(np.matrix(U[t,x]).getH(),top_neighbours[t,x])  +np.matrix(np.dot(np.matrix(U[t,x]).getH(),top_neighbours[t,x])).getH()
                G+= np.dot(np.matrix(U[t,x]).getH(),right_neighbours[t,x])  + np.matrix(np.dot(np.matrix(U[t,x]).getH(),right_neighbours[t,x])).getH()
        return -beta/2* np.trace(G)
a = Chiral(2,2)
#print(Chiral.action(a.U,a.beta))
print(a.action_2(a.U))
#print(Chiral.action(a.U,a.beta))
