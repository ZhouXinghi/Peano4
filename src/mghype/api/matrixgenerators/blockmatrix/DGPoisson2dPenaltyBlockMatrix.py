# from block_matrix import BlockMatrix
# from nodal_basis import NodalBasisGaussLegendre
# from finiteelement import DGScalarElement, DGFacetElement
# from quadrature import GaussianQuadrature, Quadrature2d
# import numpy as np

from .BlockMatrix import BlockMatrix
from .nodal_basis import NodalBasisGaussLegendre
from .finiteelement import DGScalarElement, DGFacetElement
from .quadrature import GaussianQuadrature, Quadrature2d
import numpy as np

class DGPoisson2dPenaltyBlockMatrix(BlockMatrix):
    """Classical DG discretisation of 2d classical Poisson formulation with spurious facet function spaces

    Order of function spaces (used in the output dictionaries):
    0 |  0   1   2   3   4   5
    u | q^+ q^- q^f u^+ u^- u^f
    """

    def __init__(self, r_press, basis=NodalBasisGaussLegendre(2)):
        """Initialise instance
        :arg r_press: polynomial degree of pressure space
        :arg r_vel: polynomial degree of velocity
        :arg basis: nodal basis
        """
        dim = 2
        super().__init__(dim, basis)
        self.basis = basis
        self._L_C = [r_press]
        self._L_F = 6 * [r_press]
        self._fs_labels_C = ["u"]
        self._fs_labels_F = ["q+", "q-", "qf", "u+", "u-", "uf"] # projections of the normal derivative onto facets and "flux"

        # Penalty method parameters
        self._theta = 1.0
        self._gamma_F = 2 # this gives 2/h when scaled
        print('Warning: penalty parameter gamma_F is set to 2/h and does not depend on the polynomial degree!')
        # self._gamma_F = r_press * (r_press + 1) # this should be degree * (degree + 1) / h. 
        #                                         # We remove /h here because the contribution of this penalty term cannot be scaled separately from the first penalty term at the moment
 
        # All pairs of function spaces that couple
        self._coupling_CC = [(0, 0)] # _u_u
        
        self._coupling_CF = [(0, 2), # _u_qf
                             (0, 5)] # _u_uf
        
        self._coupling_FC = [(0, 0), # _q+_u
                             (1, 0), # _q-_u
                             (3, 0), # _u+_u
                             (4, 0)  # _u-_u
                             ]
        
        self._coupling_FF = [(0, 0), # _q+_q+
                             (1, 1), # _q-_q-
                             (2, 0), # _qf_q+
                             (2, 1), # _qf_q-
                             (2, 2), # _qf_qf
                             (3, 3), # _u+_u+
                             (4, 4), # _u-_u-
                             (5, 3), # _uf_u+
                             (5, 4), # _uf_u-
                             (5, 5)  # _uf_uf
                             ]
        
        # Lists of entity types
        self._ent_types_CF_FC = ["left", "right", "bottom", "top"] # ! All the facet integrals will be now assembled for the left facet only !
        self._ent_types_FF = self._ent_types_CF_FC # ! All the facet integrals will be now assembled for the left facet only !
        #self._ent_types_FF = ["ih", "iv", "bl", "br", "bb", "bt"] # "interior horizontal/vertical" and "boundary left/right/bottom/top"

        # Cell elements
        self._fe_C = [None] * len(self._L_C)
        self._Xi_C = [None] * len(self._L_C)
        for i in range(len(self._L_C)):
            self._fe_C[i] = DGScalarElement(dim, self._L_C[i])
            self._Xi_C[i] = self._fe_C[i].nodes

        # Facet elements
        self._fe_F = [None] * len(self._L_F)
        self._Xi_F = [None] * len(self._L_F)
        for i in range(len(self._L_F)):
            self._fe_F[i] = DGFacetElement(dim, self._L_F[i])
            self._Xi_F[i] = self._fe_F[i].nodes
        

        self._A_CC = {fs_pair: self._assemble_A_CC(fs_pair) for fs_pair in self._coupling_CC}

        self._A_CF = {fs_pair: {ent_type: self._assemble_A_CF(fs_pair, ent_type) for ent_type in self._ent_types_CF_FC} 
                          for fs_pair in self._coupling_CF
                         }

        self._A_FC = {fs_pair: {ent_type: self._assemble_A_FC(fs_pair, ent_type) for ent_type in self._ent_types_CF_FC} 
                          for fs_pair in self._coupling_FC
                         }

        self._A_FF = {fs_pair: {ent_type: self._assemble_A_FF(fs_pair, ent_type) for ent_type in self._ent_types_FF} 
                          for fs_pair in self._coupling_FF
                         }
        
        # Matrix for rhs
        # Remake this according to the note; now only works for the mass matrix for the Poisson equation with interior penalty method
        self._RHS_CC = self._assemble_RHS_CC()
        
        #MG = MatrixGenerator(r_)

    ###################################
    ### Assemble weak form for A_CC ###
    ###################################
    def _assemble_A_CC(self, fs_pair):
        """Assemble weak form for A_{CC}^{(k,ell)}
        :arg fs_pair: pair of cell function space indices that couple
        k = fs_pair[0]: first function space (to-space): a cell function space
        ell = fs_pair[1]: second function space (from-space): a cell function space
        """
        #assert fs_pair in self._coupling_CC

        k = fs_pair[0]
        ell = fs_pair[1]
        ndof_k = self._fe_C[k].ndof
        ndof_ell = self._fe_C[ell].ndof
        A_k_ell = np.empty((ndof_k, ndof_ell))
        gamma = (0, 0) # no sacling, as the contribution h^2 from dxdy is compensated by 1/h^2 from derivatives \nambla u \cdot \nabla v

        # Quadrature rule: based on the max degree of basis in the function spaces pair
        max_degree = max(self._fe_C[k].degree, self._fe_C[ell].degree)
        quad2d = Quadrature2d(GaussianQuadrature(max_degree + 1))

        if fs_pair == (0, 0):
            #func_k = self._fe_C[k].evaluate
            #func_ell = self._fe_C[ell].evaluate
            factor = 1.0
            
            # Assemble local matrix A^(k,ell)
            for i in range(ndof_k):
                for j in range(ndof_ell):
                    A_k_ell[i, j] = 0.0
                    for q, w in zip(quad2d.points, quad2d.weights):
                        A_k_ell[i, j] += factor * w * (  self._fe_C[k].evaluate_derivative(i, 0, q) * self._fe_C[ell].evaluate_derivative(j, 0, q)
                                                         +
                                                         self._fe_C[k].evaluate_derivative(i, 1, q) * self._fe_C[ell].evaluate_derivative(j, 1, q)
                                                      )
        else:
            A_k_ell = np.zeros((k, ell))

        return {gamma: A_k_ell}

    ###################################
    ### Assemble weak form for A_CF ###
    ###################################   
    def _assemble_A_CF(self, fs_pair, relative_position):
        """Assemble weak form for A_{CF}^{(k,ell)}
        :arg fs_pair: pair of function space indices that couple
        k = fs_pair[0]: first function space (to-space): a cell function space
        ell = fs_pair[1]: second function space (from-space): a facet function space
        :arg relative_position: relative position of the facet within a cell given by
        the "type of entity"
        """
        k = fs_pair[0]
        ell = fs_pair[1]
        ndof_k = self._fe_C[k].ndof
        ndof_ell = self._fe_F[ell].ndof
        A_k_ell = np.empty((ndof_k, ndof_ell))
        gamma = (0, 0)

        # Quadrature rule: based on the max degree of basis in the function spaces pair
        max_degree = max(self._fe_C[k].degree, self._fe_F[ell].degree)
        quad = GaussianQuadrature(max_degree + 1)

        # Matrix restr applied to [1, q] gives [-1, q], [1, q], [q, -1] or [q, 1] -- 
        # 2D coordinates on a facet depending on its relative location
        # Is's used to later restrict a cell-function to a facet for integrating over facet
        if relative_position == "left":
            restr = np.array([[-1.0, 0.0],
                              [ 0.0, 1.0]])
            n_sign = -1.0 # nx, ny = -1.0, 0.0 # outward unit normal components
        elif relative_position == "right":
            restr = np.array([[1.0, 0.0],
                              [0.0, 1.0]])
            n_sign = +1.0 # nx, ny = +1.0, 0.0
        elif relative_position == "bottom":
            restr = np.array([[ 0.0, 1.0],
                              [-1.0, 0.0]])
            n_sign = -1.0 # nx, ny = 0.0, -1.0
        elif relative_position == "top":
            restr = np.array([[0.0, 1.0],
                              [1.0, 0.0]])
            n_sign = +1.0 # nx, ny = 0.0, +1.0
        else:
            raise ValueError("Invalid value of relative_position")        
        
        # The following sets factors and integrands for all the pairs (k,ell);
        if fs_pair == (0, 2): # -int_F q^f v^u dF
            
            if relative_position in ["left", "right"]:
                gamma = (0, 1) # scale by hy
            else: # ["top", "bottom"]
                gamma = (1, 0) # scale by hx

            # Define the integrand function for different
            def func_k(i, q):
                facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                return self._fe_C[k].evaluate(i, facet_coord) # restriction of the cell-function to a facet
            func_ell = self._fe_F[ell].evaluate

            factor = -1.0 * 0.5 * n_sign # -1.0: sign of the flux integral in the weak formulation
                                              #  0.5: mesh-scaling is h/2 

            # Assemble local matrix A^(k,ell)
            for i in range(ndof_k):
                for j in range(ndof_ell):
                    A_k_ell[i, j] = 0.0
                    for q, w in zip(quad.points, quad.weights):
                        A_k_ell[i, j] += factor * w * func_k(i, q) * func_ell(j, q)

        elif fs_pair == (0, 5):
            ### Compute first penalty term \theta * <[[u]], {n*\nabla v}>
            ### it transforms to \theta * <u^f, \nabla v> 
            gamma = (0, 0) # factor 0.5 and saling h compensated by \grad(v)
                           # assuming hx = hy; correct later: should be (\pm 1, \pm 1) in general

            # Define the integrand function for different
            def func_k(i, q):
                facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                if relative_position in ["left", "right"]:
                    return self._fe_C[k].evaluate_derivative(i, 0, facet_coord) # restriction of the cell-function to a facet
                else:
                    return self._fe_C[k].evaluate_derivative(i, 1, facet_coord) # restriction of the cell-function to a facet
            func_ell = self._fe_F[ell].evaluate

            factor = 1.0
            
            # Assemble the penalty term matrix
            A_penalty_1 = np.zeros_like(A_k_ell)
            for i in range(ndof_k):
                for j in range(ndof_ell):
                    for q, w in zip(quad.points, quad.weights):
                        A_penalty_1[i, j] += self._theta * factor * w * func_k(i, q) * func_ell(j, q)

            ### Compute second penalty term \gamma_F * <[[u], [[v]]]>
            # # Direction of integration
            # if relative_position in ["bottom", "right"]:
            #     sg = +1.0
            # else: # ["left", "top"]:
            #     sg = -1.0

            # ! SCALING IS REMOVED FOR THE SECOND PENALTY TERMS ! SHOULD BE A SEPARATE KEY->VALUES IN THE OUTPUT DICTIONARY
            # if relative_position in ["left", "right"]:
            #     gamma = (0, 1) # scale by hy
            # else:
            #     gamma = (1, 0) # scale by hx

            # Define the integrand function for different
            def func_k(i, q):
                facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                return self._fe_C[k].evaluate(i, facet_coord) # restriction of the cell-function to a facet
            func_ell = self._fe_F[ell].evaluate

            factor = 0.5 * 2.0 * n_sign # 0.5, because mesh-scaling is h/2
                                        # 2.0 comes from the jump in the weak form. It should be 1.0 on the boundary!

            # Assemble the penalty term matrix
            A_penalty_2 = np.zeros_like(A_k_ell)
            for i in range(ndof_k):
                for j in range(ndof_ell):
                    for q, w in zip(quad.points, quad.weights):
                        A_penalty_2[i, j] += self._gamma_F * factor * w * func_k(i, q) * func_ell(j, q)

            # Assemble local matrix A^(k,ell)
            # !!! Now both components hve the same scaling (0,0)
            # !!! It should ultimately be two different dictionary entries with different scaling:
            # !!! {gamma1: A_penalty_1, gamma2: A_penalty_2}
            A_k_ell = A_penalty_1 + A_penalty_2
            #print(A_penalty_2)
            pass

        else:
            A_k_ell = np.zeros((ndof_k, ndof_ell))
        
        #print(k, ell, relative_position, "\n", A_k_ell, "\n")
        return {gamma: A_k_ell}

    ###################################
    ### Assemble weak form for A_FC ###
    ###################################
    # Matrices that form projections being multiplied by inv(M) from _A_FF
    def _assemble_A_FC(self, fs_pair, relative_position):
        """Assemble weak form for A_{CF}^{(k,ell)}
        :arg fs_pair: pair of function space indices that couple
        k = fs_pair[0]: first function space (to-space): a facet function space
        ell = fs_pair[1]: second function space (from-space): a cell function space
        :arg relative_position: relative position of the facet within a cell
        """
        k = fs_pair[0]
        ell = fs_pair[1]
        ndof_k = self._fe_F[k].ndof
        ndof_ell = self._fe_C[ell].ndof
        A_k_ell = np.empty((ndof_k, ndof_ell))

        # Quadrature rule: based on the max degree of basis in the function spaces pair
        max_degree = max(self._fe_F[k].degree, self._fe_C[ell].degree)
        quad = GaussianQuadrature(max_degree + 1)

        # Matrix restr applied to [1, q] gives [-1, q], [1, q], [q, -1] or [q, 1] -- 
        # 2D coordinates on a facet depending on its relative location
        # Is's used to later restrict a cell-function to a facet for integrating over facet
        if relative_position == "left":
            restr = np.array([[-1.0, 0.0],
                              [ 0.0, 1.0]])
            # nx, ny = -1.0, 0.0 # outward unit normal components
        elif relative_position == "right":
            restr = np.array([[1.0, 0.0],
                              [0.0, 1.0]])
            # nx, ny = +1.0, 0.0
        elif relative_position == "bottom":
            restr = np.array([[ 0.0, 1.0],
                              [-1.0, 0.0]])
            # nx, ny = 0.0, -1.0
        elif relative_position == "top":
            restr = np.array([[0.0, 1.0],
                              [1.0, 0.0]])
            # nx, ny = 0.0, +1.0
        else:
            raise ValueError("Invalid value of relative_position")
        
        factor = -1.0 # minus from the def. of projection: M*q^\pm + A*u = 0, so A = -integral(...)
        gamma = (0, 0) # no scaling as it's compensated by the derivative dp/dx or dp/dy
        
        if fs_pair == (0, 0): # A_q+_u
            if relative_position in ["left", "bottom"]: # "+"-projection from a give cell is defined on the left and bottom facets
                # ! Now, q-projections are defined as the restriction of the gradient onto a facet, not the normal component. I.e. the direction of the normal vector is not taken into account. 
                # factor = factor * (-1.0) # * -1 = normal_component: nx for the left or ny for the bottom
                # Define the integrand functions
                func_k = self._fe_F[k].evaluate
                def func_ell(i, q):
                    facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                    if relative_position == "left": 
                        # Define the integrand function: projecting x-derivative of p
                        return self._fe_C[ell].evaluate_derivative(i, 0, facet_coord) # restriction of the cell-function to a facet
                    else: # relative_position == "bottom"
                        # projecting y-derivative of p
                        return self._fe_C[ell].evaluate_derivative(i, 1, facet_coord)

                # Assemble local matrix A^(k,ell)
                for i in range(ndof_k):
                    for j in range(ndof_ell):
                        A_k_ell[i, j] = 0.0
                        for q, w in zip(quad.points, quad.weights):
                            A_k_ell[i, j] += factor * w * func_k(i, q) * func_ell(j, q)

            else: # return 0 for the right and top
                A_k_ell = np.zeros((ndof_k, ndof_ell))

        elif fs_pair == (1, 0): # A_q-_u
            if relative_position in ["right", "top"]: # "-"-projection from a give cell is defined on the right and top facets
                # factor = factor * (1.0) # * 1 = normal_component: nx for the right or ny for the top
                # Define the integrand functions
                func_k = self._fe_F[k].evaluate
                def func_ell(i, q):
                    facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                    if relative_position == "right": 
                        # Define the integrand function: projecting x-derivative of p
                        return self._fe_C[ell].evaluate_derivative(i, 0, facet_coord) # restriction of the cell-function to a facet
                    else: # relative_position == "top"
                        # projecting y-derivative of p
                        return self._fe_C[ell].evaluate_derivative(i, 1, facet_coord)
                    
                # Assemble local matrix A^(k,ell)
                for i in range(ndof_k):
                    for j in range(ndof_ell):
                        A_k_ell[i, j] = 0.0
                        for q, w in zip(quad.points, quad.weights):
                            A_k_ell[i, j] += factor * w * func_k(i, q) * func_ell(j, q)

            else: # return 0 for the left and bottom
                A_k_ell = np.zeros((ndof_k, ndof_ell))

        elif fs_pair == (3, 0): # A_u+_u
            factor = -0.5 # minus from the def. of projection: M*u^\pm + A*u = 0, so A = -integral(...)
            factor = factor * (-1.0) # ! now, the plus-prujections are {(-1)*solution} and the minus-projections are {(+1)*solution}
                                     # ! This is to have the jumps (u^- - u^+) in the penalty formulation
            if relative_position in ["left", "bottom"]:
                if relative_position == "left":
                    gamma = (0, 1) # scale by hy
                else:
                    gamma = (1, 0) # scale by hx

                # Define the integrand functions
                func_k = self._fe_F[k].evaluate
                def func_ell(i, q):
                    facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                    # Define the integrand function: projecting solution u (or p)
                    return self._fe_C[ell].evaluate(i, facet_coord) # restriction of the cell-function to a facet
                
                # Assemble local matrix A^(k,ell)
                for i in range(ndof_k):
                    for j in range(ndof_ell):
                        A_k_ell[i, j] = 0.0
                        for q, w in zip(quad.points, quad.weights):
                            A_k_ell[i, j] += factor * w * func_k(i, q) * func_ell(j, q)

            else: # return 0 for the right and top
                A_k_ell = np.zeros((ndof_k, ndof_ell))

        # REMAKE THIS, PUT (3, 0) AND (4, 0) TOGETHER AS THE INTEGRANDS ARE THE SAME    
        elif fs_pair == (4, 0): # A_u-_u
            factor = -0.5 # minus from the def. of projection: M*u^\pm + A*u = 0, so A = -integral(...) 
            if relative_position in ["right", "top"]:
                if relative_position == "right":
                    gamma = (0, 1) # scale by hy
                else:
                    gamma = (1, 0) # scale by hx
                
                # Define the integrand functions
                func_k = self._fe_F[k].evaluate
                def func_ell(i, q):
                    facet_coord = np.dot(restr, np.array([1.0, q])) # gives [-1, q], [1, q], [q, -1] or [q, 1]
                    # Define the integrand function: projecting solution u (or p)
                    return self._fe_C[ell].evaluate(i, facet_coord) # restriction of the cell-function to a facet
                
                # Assemble local matrix A^(k,ell)
                for i in range(ndof_k):
                    for j in range(ndof_ell):
                        A_k_ell[i, j] = 0.0
                        for q, w in zip(quad.points, quad.weights):
                            A_k_ell[i, j] += factor * w * func_k(i, q) * func_ell(j, q)

            else: # return 0 for the left and bottom
                A_k_ell = np.zeros((ndof_k, ndof_ell))

        else:
            A_k_ell = np.zeros((ndof_k, ndof_ell))
        
        #print(k, ell, relative_position, "\n", A_k_ell, "\n")
        return {gamma: A_k_ell}
    
    ###################################
    ### Assemble weak form for A_FF ###
    ###################################
    def _assemble_A_FF(self, fs_pair, relative_position):
        """Assemble weak form for A_{FF}^{(k,ell)}
        :arg fs_pair: pair of function space indices that couple
        k = fs_pair[0]: first function space (to-space): a facet function space
        ell = fs_pair[1]: second function space (from-space): a facet function space
        :arg relative_position: relative position of the facet within a cell
        """
        k = fs_pair[0]
        ell = fs_pair[1]
        ndof_k = self._fe_F[k].ndof
        ndof_ell = self._fe_F[ell].ndof
        A_k_ell = np.empty((ndof_k, ndof_ell))
        gamma = (0, 0)

        # Quadrature rule: based on the max degree of basis in the function spaces pair
        max_degree = max(self._fe_F[k].degree, self._fe_F[ell].degree)
        quad = GaussianQuadrature(max_degree + 1)

        if fs_pair in [(0,0), (1,1), (3,3), (4,4)]: # Facet space mass matrices A_q+_q+, A_q-_q-, A_u+_u+, A_u-_u-
            func_k = self._fe_F[k].evaluate
            func_ell = self._fe_F[ell].evaluate
            factor = 0.5

            if relative_position in ["left", "right"]:
                gamma = (0, 1) # scale by hy
            else:
                gamma = (1, 0) # scale by hx

            # Assemble local matrix A^(k,ell)
            for i in range(ndof_k):
                for j in range(ndof_ell):
                    A_k_ell[i, j] = 0.0
                    for q, w in zip(quad.points, quad.weights):
                        A_k_ell[i, j] += factor * w * func_k(i, q) * func_ell(j, q)

        elif fs_pair == (2, 2): # A_qf_qf
            gamma = (0, 0)
            A_k_ell = np.eye(ndof_k)

        elif fs_pair in [(2,0), (2,1)]: # A_q+_qf, A_q-_qf
            gamma = (0, 0)
            A_k_ell = -0.5*np.eye(ndof_k)

        else:
            gamma = (0, 0)
            A_k_ell = np.zeros((ndof_k, ndof_ell))

        #print(k, ell, relative_position, "\n", A_k_ell, "\n")
        return {gamma: A_k_ell}

    ###########################
    ### Assemble RHS matrix ###
    ###########################
    # Remake this according to the note; now only works for the mass matrix for the Poisson equation with interior penalty method
    def _assemble_RHS_CC(self):
        """
        Assemble the mass matrix for the rhs.
        It is then left-multiplied by the vector of f-function evaluated at nodal points to generate the rhs vector.
        """
        
        k   = 0 # cell-cell mass matrix,
        ell = 0 # so 0-0 pair of cell function spaces
        ndof_k = self._fe_C[k].ndof
        ndof_ell = self._fe_C[ell].ndof
        MassMatrix = np.empty((ndof_k, ndof_ell))
        gamma = (1, 1) # scaling from dx dy only => hx*hy

        # Quadrature rule: based on the max degree of basis in the function spaces pair
        max_degree = max(self._fe_C[k].degree, self._fe_C[ell].degree)
        quad2d = Quadrature2d(GaussianQuadrature(max_degree + 1))

        factor = 0.25
        
        # Assemble local mass matrix
        for i in range(ndof_k):
            for j in range(ndof_ell):
                MassMatrix[i, j] = 0.0
                for q, w in zip(quad2d.points, quad2d.weights):
                    MassMatrix[i, j] += factor * w * self._fe_C[k].evaluate(i, q) * self._fe_C[ell].evaluate(j, q)

        return {gamma: MassMatrix}    


# r_press = 1
# BM = DGPoisson2dPenaltyBlockMatrix(r_press)
# pass